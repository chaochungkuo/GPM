import abc
import gc
from typing import Callable, Dict, Tuple

import lightning.pytorch as pl
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torchmetrics import MetricCollection
from torchmetrics.classification import MulticlassF1Score

from cellnet.tabnet.tab_network import TabNet


class MLP(nn.Module):

    def __init__(
        self,
        input_dim: int,
        output_dim: int,
        hidden_size: int = 128,
        n_hidden: int = 8,
        dropout: float = 0.1
    ):
        super().__init__()
        assert n_hidden >= 1

        modules = [
            nn.Linear(input_dim, hidden_size),
            nn.BatchNorm1d(hidden_size),
            nn.SiLU(),
            nn.Dropout(p=dropout)
        ]
        for _ in range(1, n_hidden):
            modules += [
                nn.Linear(hidden_size, hidden_size),
                nn.BatchNorm1d(hidden_size),
                nn.SiLU(),
                nn.Dropout(p=dropout)
            ]

        self.encoder = nn.Sequential(*modules)
        self.linear = nn.Linear(hidden_size, output_dim)

    def forward(self, x):
        return self.linear(self.encoder(x))


def augment_data(x: torch.Tensor, augmentation_vectors: torch.Tensor):
    augmentations = augmentation_vectors[
        torch.randint(0, augmentation_vectors.shape[0], (x.shape[0], ), device=x.device), :
    ]
    sign = 2. * (torch.bernoulli(.5 * torch.ones(x.shape[0], 1, device=x.device)) - .5)

    return torch.clamp(x + (sign * augmentations), min=0., max=9.)


class BaseClassifier(pl.LightningModule, abc.ABC):

    classifier: nn.Module  # classifier mapping von gene_dim to type_dim - outputs logits

    def __init__(
        self,
        # fixed params
        gene_dim: int,
        type_dim: int,
        class_weights: np.ndarray,
        child_matrix: np.ndarray,
        # params from datamodule
        train_set_size: int,
        val_set_size: int,
        batch_size: int,
        # model specific params
        learning_rate: float = 0.005,
        weight_decay: float = 0.05,
        optimizer: Callable[..., torch.optim.Optimizer] = torch.optim.AdamW,
        lr_scheduler: Callable = None,
        lr_scheduler_kwargs: Dict = None,
        gc_frequency: int = 10
    ):
        super(BaseClassifier, self).__init__()

        self.gene_dim = gene_dim
        self.type_dim = type_dim
        self.train_set_size = train_set_size
        self.val_set_size = val_set_size
        self.batch_size = batch_size
        self.gc_freq = gc_frequency
        self.lr = learning_rate
        self.weight_decay = weight_decay
        self.optim = optimizer
        self.lr_scheduler = lr_scheduler
        self.lr_scheduler_kwargs = lr_scheduler_kwargs

        metrics = MetricCollection({
            'f1_micro': MulticlassF1Score(num_classes=type_dim, average='micro'),
            'f1_macro': MulticlassF1Score(num_classes=type_dim, average='macro'),
        })
        self.train_metrics = metrics.clone(prefix='train_')
        self.val_metrics = metrics.clone(prefix='val_')
        self.test_metrics = metrics.clone(prefix='test_')
        self.register_buffer('class_weights', torch.tensor(class_weights.astype('f4')))
        self.register_buffer('child_lookup', torch.tensor(child_matrix.astype('i8')))

    @abc.abstractmethod
    def _step(self, batch, training=True) -> Tuple[torch.Tensor, torch.Tensor]:
        pass

    def hierarchy_correct(self, preds: torch.Tensor, targets: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        pred_is_child_node_or_node = torch.gt(
            torch.sum(self.child_lookup[targets, :] * F.one_hot(preds, self.type_dim), dim=1), 0
        )

        return (
            torch.where(pred_is_child_node_or_node, targets, preds),  # corrected preds
            torch.where(pred_is_child_node_or_node, preds, targets)  # corrected targets
        )

    def on_after_batch_transfer(self, batch, dataloader_idx):
        with torch.no_grad():
            batch = batch[0]
            batch['cell_type'] = torch.squeeze(batch['cell_type'])

        return batch

    def forward(self, x: torch.Tensor):
        return self.classifier(x)

    def training_step(self, batch, batch_idx):
        preds, loss = self._step(batch, training=True)
        self.log('train_loss', loss)
        f1_macro = self.train_metrics['f1_macro'](preds, batch['cell_type'])
        f1_micro = self.train_metrics['f1_micro'](preds, batch['cell_type'])
        self.log('train_f1_macro_step', f1_macro)
        self.log('train_f1_micro_step', f1_micro)
        if batch_idx % self.gc_freq == 0:
            gc.collect()

        return loss

    def validation_step(self, batch, batch_idx):
        preds, loss = self._step(batch, training=False)
        self.log('val_loss', loss)
        self.val_metrics['f1_macro'].update(preds, batch['cell_type'])
        self.val_metrics['f1_micro'].update(preds, batch['cell_type'])
        if batch_idx % self.gc_freq == 0:
            gc.collect()

    def test_step(self, batch, batch_idx):
        preds, loss = self._step(batch, training=False)
        self.log('test_loss', loss)
        self.test_metrics['f1_macro'].update(preds, batch['cell_type'])
        self.test_metrics['f1_micro'].update(preds, batch['cell_type'])
        if batch_idx % self.gc_freq == 0:
            gc.collect()

    def on_train_epoch_end(self) -> None:
        self.log('train_f1_macro_epoch', self.train_metrics['f1_macro'].compute())
        self.train_metrics['f1_macro'].reset()
        self.log('train_f1_micro_epoch', self.train_metrics['f1_micro'].compute())
        self.train_metrics['f1_micro'].reset()
        gc.collect()

    def on_validation_epoch_end(self) -> None:
        f1_macro = self.val_metrics['f1_macro'].compute()
        self.log('val_f1_macro', f1_macro)
        self.log('hp_metric', f1_macro)
        self.val_metrics['f1_macro'].reset()
        self.log('val_f1_micro', self.val_metrics['f1_micro'].compute())
        self.val_metrics['f1_micro'].reset()
        gc.collect()

    def on_test_epoch_end(self) -> None:
        self.log('test_f1_macro', self.test_metrics['f1_macro'].compute())
        self.test_metrics['f1_macro'].reset()
        self.log('test_f1_micro', self.test_metrics['f1_micro'].compute())
        self.test_metrics['f1_micro'].reset()
        gc.collect()

    def configure_optimizers(self):
        optimizer_config = {'optimizer': self.optim(self.parameters(), lr=self.lr, weight_decay=self.weight_decay)}
        if self.lr_scheduler is not None:
            lr_scheduler_kwargs = {} if self.lr_scheduler_kwargs is None else self.lr_scheduler_kwargs
            interval = lr_scheduler_kwargs.pop('interval', 'epoch')
            monitor = lr_scheduler_kwargs.pop('monitor', 'val_loss_epoch')
            frequency = lr_scheduler_kwargs.pop('frequency', 1)
            scheduler = self.lr_scheduler(optimizer_config['optimizer'], **lr_scheduler_kwargs)
            optimizer_config['lr_scheduler'] = {
                'scheduler': scheduler,
                'interval': interval,
                'monitor': monitor,
                'frequency': frequency
            }

        return optimizer_config


class LinearClassifier(BaseClassifier):

    def __init__(
        self,
        # fixed params
        gene_dim: int,
        type_dim: int,
        class_weights: np.ndarray,
        child_matrix: np.ndarray,
        # params from datamodule
        train_set_size: int,
        val_set_size: int,
        batch_size: int,
        # model specific params
        learning_rate: float = 0.005,
        weight_decay: float = 0.1,
        use_class_weights: bool = True,
        optimizer: Callable[..., torch.optim.Optimizer] = torch.optim.AdamW,
        lr_scheduler: Callable = None,
        lr_scheduler_kwargs: Dict = None,
    ):
        super(LinearClassifier, self).__init__(
            gene_dim=gene_dim,
            type_dim=type_dim,
            class_weights=class_weights,
            child_matrix=child_matrix,
            train_set_size=train_set_size,
            val_set_size=val_set_size,
            batch_size=batch_size,
            learning_rate=learning_rate,
            weight_decay=weight_decay,
            optimizer=optimizer,
            lr_scheduler=lr_scheduler,
            lr_scheduler_kwargs=lr_scheduler_kwargs
        )
        self.save_hyperparameters(ignore=['class_weights', 'parent_matrix', 'child_matrix'])

        self.use_class_weights = use_class_weights
        self.classifier = nn.Linear(gene_dim, type_dim)

    def _step(self, batch, training=True):
        x = batch['X']
        logits = self(x)
        with torch.no_grad():
            preds = torch.argmax(logits, dim=1)
            preds_corrected, targets_corrected = self.hierarchy_correct(preds, batch['cell_type'])
        if training:
            if self.use_class_weights:
                loss = F.cross_entropy(logits, batch['cell_type'], weight=self.class_weights)
            else:
                loss = F.cross_entropy(logits, batch['cell_type'])
        else:
            loss = F.cross_entropy(logits, targets_corrected)

        return preds_corrected, loss

    def predict_step(self, batch, batch_idx, dataloader_idx=None):
        x = batch['X']
        return F.softmax(self(x), dim=1)


class TabnetClassifier(BaseClassifier):

    def __init__(
        self,
        # fixed params
        gene_dim: int,
        type_dim: int,
        class_weights: np.ndarray,
        child_matrix: np.ndarray,
        augmentations: np.ndarray,
        # params from datamodule
        train_set_size: int,
        val_set_size: int,
        batch_size: int,
        # model specific params
        learning_rate: float = 0.005,
        weight_decay: float = 0.05,
        use_class_weights: bool = True,
        optimizer: Callable[..., torch.optim.Optimizer] = torch.optim.AdamW,
        lr_scheduler: Callable = None,
        lr_scheduler_kwargs: Dict = None,
        # tabnet params
        lambda_sparse: float = 1e-5,
        n_d: int = 128,
        n_a: int = 64,
        n_steps: int = 3,
        gamma: float = 1.3,
        n_independent: int = 5,
        n_shared: int = 3,
        epsilon: float = 1e-15,
        virtual_batch_size: int = 256,
        momentum: float = 0.02,
        mask_type: str = 'entmax',
        augment_training_data: bool = True,
    ):
        super(TabnetClassifier, self).__init__(
            gene_dim=gene_dim,
            type_dim=type_dim,
            class_weights=class_weights,
            child_matrix=child_matrix,
            train_set_size=train_set_size,
            val_set_size=val_set_size,
            batch_size=batch_size,
            learning_rate=learning_rate,
            weight_decay=weight_decay,
            optimizer=optimizer,
            lr_scheduler=lr_scheduler,
            lr_scheduler_kwargs=lr_scheduler_kwargs
        )
        self.save_hyperparameters(ignore=['class_weights', 'child_matrix', 'augmentations'])

        self.lambda_sparse = lambda_sparse
        classifier = TabNet(
            input_dim=gene_dim,
            output_dim=type_dim,
            n_d=n_d,
            n_a=n_a,
            n_steps=n_steps,
            gamma=gamma,
            n_independent=n_independent,
            n_shared=n_shared,
            epsilon=epsilon,
            virtual_batch_size=virtual_batch_size,
            momentum=momentum,
            mask_type=mask_type,
        )
        self.classifier = classifier

        self.use_class_weights = use_class_weights
        self.augment_training_data = augment_training_data
        if self.augment_training_data:
            self.register_buffer('augmentations', torch.tensor(augmentations.astype('f4')))

        self.predict_bottleneck = False

    def _step(self, batch, training=True):
        x = augment_data(batch['X'], self.augmentations) if (self.augment_training_data and training) else batch['X']
        logits, m_loss = self(x)
        with torch.no_grad():
            preds = torch.argmax(logits, dim=1)
            preds_corrected, targets_corrected = self.hierarchy_correct(preds, batch['cell_type'])
        if training:
            if self.use_class_weights:
                loss = F.cross_entropy(logits, batch['cell_type'],
                                       weight=self.class_weights) - self.lambda_sparse * m_loss
            else:
                loss = F.cross_entropy(logits, batch['cell_type']) - self.lambda_sparse * m_loss
        else:
            loss = F.cross_entropy(logits, targets_corrected)

        return preds_corrected, loss

    def predict_embedding(self, x: torch.Tensor):
        steps_output, _ = self.classifier.encoder(x)
        res = torch.sum(torch.stack(steps_output, dim=0), dim=0)
        return res

    def predict_cell_types(self, x: torch.Tensor):
        return F.softmax(self(x)[0], dim=1)

    def predict_step(self, batch, batch_idx, dataloader_idx=None):
        if batch_idx % self.gc_freq == 0:
            gc.collect()
        if self.predict_bottleneck:
            return self.predict_embedding(batch['X'])
        else:
            return self.predict_cell_types(batch['X'])


class MLPClassifier(BaseClassifier):

    def __init__(
        self,
        # fixed params
        gene_dim: int,
        type_dim: int,
        class_weights: np.ndarray,
        child_matrix: np.ndarray,
        augmentations: np.ndarray,
        # params from datamodule
        train_set_size: int,
        val_set_size: int,
        batch_size: int,
        # model specific params
        learning_rate: float = 0.005,
        weight_decay: float = 0.05,
        optimizer: Callable[..., torch.optim.Optimizer] = torch.optim.AdamW,
        lr_scheduler: Callable = None,
        lr_scheduler_kwargs: Dict = None,
        # MLP params
        hidden_size: int = 128,
        n_hidden: int = 8,
        dropout: float = 0.2,
        augment_training_data: bool = True,
    ):
        super(MLPClassifier, self).__init__(
            gene_dim=gene_dim,
            type_dim=type_dim,
            class_weights=class_weights,
            child_matrix=child_matrix,
            train_set_size=train_set_size,
            val_set_size=val_set_size,
            batch_size=batch_size,
            learning_rate=learning_rate,
            weight_decay=weight_decay,
            optimizer=optimizer,
            lr_scheduler=lr_scheduler,
            lr_scheduler_kwargs=lr_scheduler_kwargs
        )
        self.save_hyperparameters(ignore=['class_weights', 'child_matrix', 'augmentations'])

        self.classifier = MLP(
            input_dim=gene_dim,
            output_dim=type_dim,
            hidden_size=hidden_size,
            n_hidden=n_hidden,
            dropout=dropout
        )

        self.augment_training_data = augment_training_data
        if self.augment_training_data:
            self.register_buffer('augmentations', torch.tensor(augmentations.astype('f4')))

        self.predict_bottleneck = False

    def _step(self, batch, training=True):
        x = augment_data(batch['X'], self.augmentations) if (self.augment_training_data and training) else batch['X']
        logits = self(x)
        with torch.no_grad():
            preds = torch.argmax(logits, dim=1)
            preds_corrected, targets_corrected = self.hierarchy_correct(preds, batch['cell_type'])
        if training:
            loss = F.cross_entropy(logits, batch['cell_type'], weight=self.class_weights)
        else:
            loss = F.cross_entropy(logits, targets_corrected)

        return preds_corrected, loss

    def predict_step(self, batch, batch_idx, dataloader_idx=None):
        if batch_idx % self.gc_freq == 0:
            gc.collect()
        if self.predict_bottleneck:
            return self.classifier.encoder(batch['X'])
        else:
            return F.softmax(self(batch['X']), dim=1)
