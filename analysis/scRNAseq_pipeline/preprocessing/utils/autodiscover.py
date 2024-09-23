"""This module provides the implementation of sample discovery from the output of multiple scRNA-seq pipelines."""

from abc import ABC, abstractmethod
from functools import cache, reduce
from os import path, walk
from typing import Callable, Dict, List

import scanpy as sc
from utils.preprocessing_funcs import read_parsebio, splitall

# TODO: Refactor this module to Group similar functionality, and clean the code a bit.


class AutoDiscover(ABC):
    """Protocol for the AutoDiscover class. This class provides arbitary implmenetation to discover scRNA-seq samples in a directory.
    Implementation are free to use different heuristics to discover samples but they should implement the following methods:
    - get_samples_paths: Returns the paths to the samples in the directory.
    - get_sample_names: Returns the names of the samples.
    - get_read_function: Returns the function to read the samples.
    Optionally, the following methods can be implemented:
    - get_raw_samples_paths: Returns the paths to the raw samples in the directory.
    - get_raw_sample_read_function: Returns the function to read the raw samples.
    """

    def __init__(self, root_path: str) -> None:
        """Initialize the AutoDiscover instance with the root path of the directory containing the samples."""
        self.root_path: str = root_path

    @abstractmethod
    def get_samples_paths(self) -> List[str]:
        """Return a list paths to the samples in the directory."""
        pass

    @abstractmethod
    def get_read_function(self) -> Callable:
        """Return a read function that can be used to load the data."""
        pass
    
    @abstractmethod
    def get_raw_sample_read_function(self, samples=None) -> Callable:
        pass

    @abstractmethod
    def get_sample_names(self) -> Dict[str, str]:
        """Tries to automatically infer the sample names from the paths."""
        pass

    @abstractmethod
    def get_raw_sample_names(self) -> Dict[str, str]:
        pass

    @staticmethod
    def _get_names(samples):

        common_prefix = path.commonpath(samples)
        updated_ps = [path.relpath(p, common_prefix) for p in samples]
        split_ps = [splitall(p) for p in updated_ps]
        common = reduce(lambda x, y: set(x).intersection(set(y)), split_ps)

        sample_names = list(
            map(
                lambda x: [el for el in x if el not in common],
                split_ps,
            )
        )
        sample_names = list(map(lambda x: "_".join(x), sample_names))
        return dict(zip(sample_names, samples))

    @abstractmethod
    def _collect_paths():
        pass


class SingeleronAutoDiscover(AutoDiscover):
    """Implementation of the AutoDiscover protocol for the Singleron scRNA-seq data."""

    def __init__(self, root_path: str = None) -> None:
        self.root_path: str = root_path
        self.components = ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]

    def get_sample_paths(self) -> List[str]:
        paths = self._collect_paths()
        return [p for p in paths if not path.basename(p).endswith("raw")]

    def get_raw_samples_paths(self) -> List[str]:
        paths = self._collect_paths()
        return [p for p in paths if path.basename(p).endswith("raw")]

    def get_sample_names(self) -> Dict[str, str]:
        samples = self.get_sample_paths()
        return super()._get_names(samples)

    def get_raw_sample_names(self) -> List[str]:
        samples = self.get_raw_samples_paths()
        return super()._get_names(samples)

    def get_read_function(self, samples=None) -> Callable:
        return sc.read_10x_mtx

    def get_raw_sample_read_function(self, samples=None) -> Callable:
        return sc.read_10x_mtx

    @cache
    def _collect_paths(self) -> List[str]:
        sample_paths = []
        for root, dir, files in walk(self.root_path):
            if set(self.components) == set(files):
                sample_paths.append(root)
        return sample_paths


class ParseBioAutoDiscover(AutoDiscover):

    def __init__(self, root_path: str = None) -> None:
        self.root_path: str = root_path
        self.components: List[str] = [
            "all_genes.csv",
            "cell_metadata.csv",
            "count_matrix.mtx",
        ]

    def get_samples_paths(self) -> List[str]:
        paths: List[str] = self._collect_paths()
        filtered_samples: List[str] = [
            p for p in paths if not path.basename(p).endswith("DGE_filtered")
        ]
        return [
            p
            for p in filtered_samples
            if not path.basename(path.dirname(p)).endswith("all-sample")
        ]

    def get_raw_samples_paths(self) -> List[str]:
        paths = self._collect_paths()
        unfiltered_samples: List[str] = [
            p for p in paths if not path.basename(p).endswith("DGE_unfiltered")
        ]
        return [
            p
            for p in unfiltered_samples
            if not path.basename(path.basename(p)).endswith("all-sample")
        ]

    def get_sample_names(self) -> Dict[str, str]:
        sample_paths: List[str] = self.get_samples_paths()
        return super()._get_names(sample_paths)

    def get_raw_sample_names(self) -> Dict[str, str]:
        samples = self.get_raw_samples_paths()
        return super()._get_names(samples)

    def get_read_function(self, samples=None) -> Callable:
        return read_parsebio

    def get_raw_sample_read_function(self, samples=None) -> Callable:
        return read_parsebio

    @cache
    def _collect_paths(self) -> List[str]:
        sample_paths: List[str] = []
        for root, dir, files in walk(self.root_path):
            # ParseBio also also contains anndata.h5ad
            if set(self.components).issubset(set(files)):
                sample_paths.append(root)
        return sample_paths


# Documentation for CellRange output: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-overview
# https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-pipestance-structure
class CellRangerAutoDiscover(AutoDiscover):

    def __init__(self, root_path: str = None) -> None:
        self.root_path: str = root_path

    # CellRanger Multi has the prefix 'sample' in the file names in 'per_sample_outs' directory
    def get_samples_paths(self) -> List[str]:
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
    def get_raw_samples_paths(self) -> List[str]:
        all_paths = self._collect_paths()
        sample_paths = self.get_samples_paths()
        raw_paths = [p for p in all_paths if p not in sample_paths]

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

    def get_sample_names(self) -> Dict[str, str]:
        sample_paths: List[str] = self.get_samples_paths()
        return super()._get_names(sample_paths)

    def get_raw_sample_names(self) -> Dict[str, str]:
        raw_paths = self.get_raw_samples_paths()
        return super()._get_names(raw_paths)

    def get_read_function(self, samples: List[str] = None) -> Callable:

        if samples is None:
            sample_paths = self.get_samples_paths()
        else:
            sample_paths = samples

        if all([p.endswith(".h5") for p in sample_paths]):
            return sc.read_10x_h5
        else:
            return sc.read_10x_mtx

    def get_raw_read_function(self, samples: List[str] = None) -> Callable:

        if samples is None:
            raw_sample_paths = self.get_raw_samples_paths()
        else:
            raw_sample_paths = samples

        if all([p.endswith(".h5") for p in raw_sample_paths]):
            return sc.read_10x_h5
        else:
            return sc.read_10x_mtx

    @cache
    def _collect_paths(self) -> List[str]:
        exclude: List[str] = [
            "aggr",
            "multi",
            "analysis",
            # Presistance directories
            "SC_RNA_COUNTER_CS",  # CellRanger Count
            "SC_RNA_AGGREGATOR_CS",  # CellRanger Aggr
            "SC_MULTI_CS",  # CellRanger Multi
            "SC_RNA_REANALYZER_CS",  # CellRanger Reanalyze
        ]
        sample_paths: List[str] = []
        for root, dirs, files in walk(self.root_path):
            # Prune the search by execluding 'Aggr' and 'Multi' directories
            # https://stackoverflow.com/questions/19859840/excluding-directories-in-os-walk
            dirs[:] = [d for d in dirs if d not in exclude]
            for file in files:
                if file.endswith(".h5"):
                    sample_paths.append(path.join(root, file))
        return sample_paths


discover_factory: Dict[str, AutoDiscover] = {
    "10x": CellRangerAutoDiscover,
    "Singleron": SingeleronAutoDiscover,
    "PraseBio": ParseBioAutoDiscover,
}
