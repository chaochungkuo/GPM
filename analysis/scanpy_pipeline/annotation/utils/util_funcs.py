import celltypist
import scanpy as sc
from anndata import AnnData
from celltypist import models as ctypist_models


def cell_typist_annotate(
    adata: AnnData, models: list[str], inplace=True, COUNTS_LAYER="counts"
) -> AnnData:
    if len(models) == 0:
        raise ValueError("The models list are empty, enter valid model names.")

    all_models = ctypist_models.models_description().model.to_list()

    for model in models:
        if model not in all_models:
            raise ValueError("{model} not found in supported cell typist models.")

    ctypist_models.download_models(force_update=True, model=models)

    adata_celltypist = adata.copy()
    adata_celltypist.X = adata.layers[COUNTS_LAYER]
    sc.pp.normalize_per_cell(adata_celltypist, counts_per_cell_after=10**4)
    sc.pp.log1p(adata_celltypist)
    adata_celltypist.X = adata_celltypist.X.toarray()

    for model in models:
        loaded_model = ctypist_models.Model.load(model=model)
        predictions = celltypist.annotate(
            adata_celltypist, model=loaded_model, majority_voting=True
        )
        predictions_adata = predictions.to_adata()
        adata.obs["celltypist_" + model + "_annotation"] = predictions_adata.obs.loc[
            adata.obs.index, "majority_voting"
        ]
        adata.obs["celltypist_" + model + "_conf_score"] = predictions_adata.obs.loc[
            adata.obs.index, "conf_score"
        ]
    if not inplace:
        return adata
