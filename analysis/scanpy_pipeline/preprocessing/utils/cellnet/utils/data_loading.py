import numpy as np
import torch

from scipy.sparse import csc_matrix, csr_matrix, issparse
from sklearn.utils import sparsefuncs
from torch.utils.data import Dataset, DataLoader, BatchSampler, SequentialSampler, RandomSampler

"""
Data streamlining.
"""


def streamline_count_matrix(x_raw, gene_names_raw, gene_names_model):
    assert len(gene_names_raw) == len(set(gene_names_raw))
    assert len(gene_names_model) == len(set(gene_names_model))
    assert len(gene_names_raw) == x_raw.shape[1]
    assert np.isin(gene_names_raw, gene_names_model).sum() == x_raw.shape[1]
    # For fast column-wise slicing matrix has to be in csc format
    assert isinstance(x_raw, csc_matrix)
    # skip zero filling missing genes if all genes are present
    if len(gene_names_raw) == len(gene_names_model):
        return x_raw.tocsr().astype('f4')
    gene_names_raw, gene_names_model = np.array(gene_names_raw), np.array(gene_names_model)
    row, col = np.empty(x_raw.nnz, dtype='i8'), np.empty(x_raw.nnz, dtype='i8')
    data = np.empty(x_raw.nnz, dtype='f4')

    ctr = 0
    for i, gene in enumerate(gene_names_model):
        if gene in gene_names_raw:
            gene_idx = int(np.where(gene == gene_names_raw)[0])
            x_col = x_raw[:, gene_idx]
            idxs_nnz = x_col.indices.tolist()
            n_nnz = len(idxs_nnz)
            col[ctr:ctr+n_nnz] = i
            row[ctr:ctr+n_nnz] = idxs_nnz
            data[ctr:ctr+n_nnz] = x_col.data
            ctr += n_nnz

    return csr_matrix(
        (data, (row, col)),
        shape=(x_raw.shape[0], len(gene_names_model)),
        dtype='f4'
    )


def sf_normalize(x):
    """Normalize each cell to have 10000 counts. """
    x = x.copy()
    counts = np.array(x.sum(axis=1))
    # avoid zero division error
    counts += counts == 0.
    # normalize to 10000. counts
    scaling_factor = 10000. / counts

    if issparse(x):
        sparsefuncs.inplace_row_scale(x, scaling_factor)
    else:
        np.multiply(x, scaling_factor.reshape((-1, 1)), out=x)

    return x


"""
Data Loaders.
"""


class CustomDataset(Dataset):

    def __init__(self, x, obs=None):
        super(CustomDataset).__init__()
        assert any([isinstance(x, np.ndarray), isinstance(x, csr_matrix)])
        self.x = x
        self.obs = obs

    def __len__(self):
        return self.x.shape[0]

    def __getitem__(self, idx):
        if isinstance(idx, int):
            idx = [idx]
        x = self.x[idx, :]
        if isinstance(x, csr_matrix):
            x = x.toarray()

        if self.obs is not None:
            # replicate merlin dataloader output format
            out = (
                {
                    'X': torch.tensor(x.squeeze()),
                    'cell_type': torch.tensor(
                        self.obs.iloc[idx]['cell_type'].cat.codes.to_numpy().reshape((-1, 1)).astype('i8')
                    )
                }, None
            )
        else:
            out = ({'X': torch.tensor(x.squeeze())}, None)

        return out


def dataloader_factory(x, obs=None, batch_size=2048, shuffle=False):
    if shuffle:
        sampler = BatchSampler(RandomSampler(range(x.shape[0])), batch_size=batch_size, drop_last=True)
    else:
        sampler = BatchSampler(SequentialSampler(range(x.shape[0])), batch_size=batch_size, drop_last=False)

    return DataLoader(CustomDataset(x, obs), sampler=sampler, batch_size=None)
