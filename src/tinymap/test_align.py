import logging

from tinymap.index import *

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("test_align")


def take_n(n, gen):
    try:
        for i in range(n):
            yield (next(gen))
    except StopIteration:
        ...


w = 10
k = 17


def load_index():
    idx = TinymapIndex(w, k, ["chr1", "chr2", "chr3", "chr10"])
    fqs = list(find_read_files(Path.cwd()))
    reads = read_fastq(fqs[0])
    return (idx, reads)
