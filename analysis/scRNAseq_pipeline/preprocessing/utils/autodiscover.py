from functools import cache
from os import path, walk
from typing import Callable, List, Protocol

import scanpy as sc
from utils.preprocessing_funcs import get_sample_name


class AutoDiscover(Protocol):

    def __init__(self, root_path: str) -> None:
        self.root_path: str = root_path

    def get_samples_paths(self) -> List[str]:
        pass

    def get_read_function(self) -> Callable:
        pass

    def get_sample_names(self) -> List[str]:
        pass


class SingeleronAutoDiscover:

    def __init__(self, root_path: str) -> None:
        self.root_path: str = root_path
        self.components = ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]

    @cache
    def _collect_paths(self) -> List[str]:
        sample_paths = []
        for root, dir, files in self.root_path:
            if set(self.components) == set(files):
                sample_paths.append(root)
        return sample_paths

    def get_samples_paths(self) -> List[str]:
        paths = self._collect_paths()
        return [p for p in paths if not path.basename(p).endswith("raw")]

    def get_sample_names(self) -> List[str]:
        samples = self.get_samples_paths()
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
