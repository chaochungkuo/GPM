"""This module provides the implementation of sample discovery from the output of multiple scRNA-seq pipelines."""

from functools import cache
from os import path, walk
from typing import Callable, List, Protocol

import scanpy as sc
from utils.preprocessing_funcs import get_sample_name


class AutoDiscover(Protocol):
    """Protocol for the AutoDiscover class. This class provides arbitary implmenetation to discover scRNA-seq samples in a directory.
    Implementation are free to use different heuristics to discover samples but they should implement the following methods:
    - get_samples_paths: Returns the paths to the samples in the directory.
    - get_sample_names: Returns the names of the samples.
    - get_read_function: Returns the function to read the samples.
    """

    def __init__(self, root_path: str) -> None:
        """Initialize the AutoDiscover instance with the root path of the directory containing the samples."""
        self.root_path: str = root_path

    def get_samples_paths(self) -> List[str]:
        """Return a list paths to the samples in the directory."""
        pass

    def get_read_function(self) -> Callable:
        """Return a read function that can be used to load the data."""
        pass

    def get_sample_names(self) -> List[str]:
        """Tries to automatically infer the sample names from the paths."""
        pass


class SingeleronAutoDiscover:
    """Implementation of the AutoDiscover protocol for the Singleron scRNA-seq data."""

    def __init__(self, root_path: str) -> None:
        self.root_path: str = root_path
        self.components = ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]

    def get_samples_paths(self) -> List[str]:
        paths = self._collect_paths()
        return [p for p in paths if not path.basename(p).endswith("raw")]

    def get_sample_names(self) -> List[str]:
        samples: List[str] = self.get_samples_paths()
        return [
            get_sample_name(s, ["raw", "filtered", "cell_calling"]) for s in samples
        ]

    def get_read_function(self) -> Callable:
        return sc.read_10x_mtx

    def get_raw_samples_paths(self) -> List[str]:
        paths = self._collect_paths()
        return [p for p in paths if path.basename(p).endswith("raw")]

    def get_raw_sample_read_function(self) -> Callable:
        return sc.read_10x_mtx

    @cache
    def _collect_paths(self) -> List[str]:
        sample_paths = []
        for root, dir, files in walk(self.root_path):
            if set(self.components) == set(files):
                sample_paths.append(root)
        return sample_paths
