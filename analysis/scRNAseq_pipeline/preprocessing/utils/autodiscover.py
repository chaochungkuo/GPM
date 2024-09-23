"""This module provides the implementation of sample discovery from the output of multiple scRNA-seq pipelines."""

from abc import ABC, abstractmethod
from collections.abc import Callable
from functools import lru_cache, reduce
from os import path, walk

import scanpy as sc
from preprocessing_funcs import read_parsebio, splitall


class AutoDiscover(ABC):
    """Protocol for the AutoDiscover class. This class provides arbitary implmenetation to discover scRNA-seq samples in a directory.
    Implementation are free to use different heuristics to discover samples but they should implement the following methods:
    - get_samples_paths: Returns the paths to the samples in the directory.
    - get_sample_names: Returns the names of the samples.
    - get_read_function: Returns the function to read the samples.
    Optionally, the following methods can be implemented:
    - get_raw_sample_paths: Returns the paths to the raw samples in the directory.
    - get_raw_sample_read_function: Returns the function to read the raw samples.
    """

    @abstractmethod
    def __init__(self, root_path: str | None = None) -> None: ...

    @abstractmethod
    def _get_sample_paths(self) -> list[str]: ...

    @abstractmethod
    def _get_raw_sample_paths(self) -> list[str]: ...

    @abstractmethod
    def read_function(self, samples: None | dict[str, str] = None) -> Callable: ...

    @abstractmethod
    def raw_read_function(self, samples: None | dict[str, str] = None) -> Callable: ...

    @abstractmethod
    def get_samples(self) -> dict[str, str]: ...

    @abstractmethod
    def get_raw_samples(self) -> dict[str, str]: ...

    @staticmethod
    def _get_names(samples) -> dict[str, str]:

        common_prefix: str = path.commonpath(samples)
        updated_ps: list[str] = [path.relpath(p, common_prefix) for p in samples]
        split_ps = [splitall(p) for p in updated_ps]
        common: set[str] | list[str] = reduce(
            lambda x, y: set(x).intersection(set(y)), split_ps
        )

        combined_names: list[str] = list(
            map(
                lambda x: [el for el in x if el not in common],
                split_ps,
            )
        )  # type: ignore
        sample_names: list[str] = ["_".join(x) for x in combined_names]
        return dict(zip(sample_names, samples))

    @abstractmethod
    @lru_cache
    def _collect_paths(self) -> list[str]: ...


class SingeleronAutoDiscover(AutoDiscover):
    """Implementation of the AutoDiscover protocol for the Singleron scRNA-seq data."""

    def __init__(self, root_path: str | None = None) -> None:
        self.root_path: str | None = root_path
        self.components = ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]

    def _get_sample_paths(self) -> list[str]:
        paths = self._collect_paths()
        return [p for p in paths if not path.basename(p).endswith("raw")]

    def _get_raw_sample_paths(self) -> list[str]:
        paths = self._collect_paths()
        return [p for p in paths if path.basename(p).endswith("raw")]

    def get_samples(self) -> dict[str, str]:
        samples = self._get_sample_paths()
        return super()._get_names(samples)

    def get_raw_samples(self) -> dict[str, str]:
        samples = self._get_raw_sample_paths()
        return super()._get_names(samples)

    def read_function(self, samples=None) -> Callable:
        _ = samples
        return sc.read_10x_mtx

    def raw_read_function(self, samples=None) -> Callable:
        _ = samples
        return sc.read_10x_mtx

    @lru_cache
    def _collect_paths(self) -> list[str]:
        sample_paths = []
        if self.root_path is not None:
            for root, _, files in walk(self.root_path):
                if set(self.components) == set(files):
                    sample_paths.append(root)
        return sample_paths


class ParseBioAutoDiscover(AutoDiscover):

    def __init__(self, root_path: str | None = None) -> None:
        self.root_path: str | None = root_path
        self.components: list[str] = [
            "all_genes.csv",
            "cell_metadata.csv",
            "count_matrix.mtx",
        ]

    def _get_sample_paths(self) -> list[str]:
        paths: list[str] = self._collect_paths()
        filtered_samples: list[str] = [
            p for p in paths if not path.basename(p).endswith("DGE_filtered")
        ]
        return [
            p
            for p in filtered_samples
            if not path.basename(path.dirname(p)).endswith("all-sample")
        ]

    def _get_raw_sample_paths(self) -> list[str]:
        paths = self._collect_paths()
        unfiltered_samples: list[str] = [
            p for p in paths if not path.basename(p).endswith("DGE_unfiltered")
        ]
        return [
            p
            for p in unfiltered_samples
            if not path.basename(path.dirname(p)).endswith("all-sample")
        ]

    def get_samples(self) -> dict[str, str]:
        sample_paths: list[str] = self._get_sample_paths()
        return super()._get_names(sample_paths)

    def get_raw_samples(self) -> dict[str, str]:
        samples = self._get_raw_sample_paths()
        return super()._get_names(samples)

    def read_function(self, samples: dict[str, str] | None = None) -> Callable:
        _ = samples
        return read_parsebio

    def raw_read_function(self, samples: dict[str, str] | None = None) -> Callable:
        _ = samples
        return read_parsebio

    @lru_cache
    def _collect_paths(self) -> list[str]:
        sample_paths: list[str] = []
        if self.root_path is not None:
            for root, _, files in walk(self.root_path):
                # ParseBio also also contains anndata.h5ad
                if set(self.components).issubset(set(files)):
                    sample_paths.append(root)
        return sample_paths


# Documentation for CellRange output: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-overview
# https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-pipestance-structure
class CellRangerAutoDiscover(AutoDiscover):

    def __init__(self, root_path: str | None = None) -> None:
        self.root_path: str | None = root_path

    # CellRanger Multi has the prefix 'sample' in the file names in 'per_sample_outs' directory
    def _get_sample_paths(self) -> list[str]:
        sample_paths = self._collect_paths()
        final_samples = []
        for sample in sample_paths:
            if not (
                "raw_feature_bc_matrix" in path.basename(sample)
                or "molecule_info.h5" in path.basename(sample)
            ):
                final_samples.append(sample)
        return final_samples

    # For CellRanger Multi, there is no raw '.h5' file, it exists as a directory.
    def _get_raw_sample_paths(self) -> list[str]:
        all_paths = self._collect_paths()
        sample_paths = self._get_sample_paths()
        raw_paths = []
        for sample in all_paths:
            if not (
                sample in sample_paths or "molecule_info.h5" in path.basename(sample)
            ):
                raw_paths.append(sample)

        # All raw samples are in '.h5' format
        if len(raw_paths) == len(sample_paths):
            return raw_paths

        # Fallback to checking directory components
        raw_paths = []
        for sample in sample_paths:
            sample_dir = path.dirname(sample)
            potential_path = path.join(sample_dir, "sample_raw_feature_bc_matrix")
            if path.exists(potential_path):
                raw_paths.append(potential_path)

        if len(raw_paths) == len(sample_paths):
            return raw_paths

        raise LookupError("Can't find all raw counterparts for all samples.")

    def get_samples(self) -> dict[str, str]:
        sample_paths: list[str] = self._get_sample_paths()
        return super()._get_names(sample_paths)

    def get_raw_samples(self) -> dict[str, str]:
        raw_paths = self._get_raw_sample_paths()
        return super()._get_names(raw_paths)

    def read_function(self, samples: dict[str, str] | None = None) -> Callable:

        if samples is None:
            sample_paths = self._get_sample_paths()
        else:
            sample_paths = list(samples.values())

        if all([p.endswith(".h5") for p in sample_paths]):
            return sc.read_10x_h5
        return sc.read_10x_mtx

    def raw_read_function(self, samples: dict[str, str] | None = None) -> Callable:

        if samples is None:
            raw_sample_paths = self._get_sample_paths()
        else:
            raw_sample_paths = list(samples.values())

        if all([p.endswith(".h5") for p in raw_sample_paths]):
            return sc.read_10x_h5
        return sc.read_10x_mtx

    @lru_cache
    def _collect_paths(self) -> list[str]:
        exclude: list[str] = [
            "aggr",
            "multi",
            "analysis",
            # Presistance directories
            "SC_RNA_COUNTER_CS",  # CellRanger Count
            "SC_RNA_AGGREGATOR_CS",  # CellRanger Aggr
            "SC_MULTI_CS",  # CellRanger Multi
            "SC_RNA_REANALYZER_CS",  # CellRanger Reanalyze
        ]
        sample_paths: list[str] = []
        if self.root_path is not None:
            for root, dirs, files in walk(self.root_path):
                # Prune the search by execluding 'Aggr' and 'Multi' directories
                # https://stackoverflow.com/questions/19859840/excluding-directories-in-os-walk
                dirs[:] = [d for d in dirs if d not in exclude]
                for file in files:
                    if file.endswith(".h5"):
                        sample_paths.append(path.join(root, file))
        return sample_paths


discover_factory: dict[str, type[AutoDiscover]] = {
    "10x": CellRangerAutoDiscover,
    "Singleron": SingeleronAutoDiscover,
    "PraseBio": ParseBioAutoDiscover,
}
