from os import path, walk
from typing import Callable, List, Protocol


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

    def get_samples_paths(self) -> List[str]:
        sample_paths = []
        for root, dir, files in self.root_path:
            if set(self.components) == set(files):
                sample_paths.append(root)

        # Filter raw samples paths
        for p in sample_paths:
            if path.basename(p).endswith("raw"):
                sample_paths.remove(p)

        return sample_paths

    def get_sample_names(self) -> List[str]:
        samples = self.get_samples_paths
        pass

    def get_read_function(self) -> Callable:
        pass

    def get_raw_samples_paths(self) -> List[str]:
        pass

    def get_raw_sample_read_function(self) -> Callable:
        pass
