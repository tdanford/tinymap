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
    idx = TinymapIndex(w, k, ["chr1"])
    dir = Path.cwd() / "data"
    fqs = list(find_read_files(dir))
    if len(fqs) != 0:
        reads = read_fastq(fqs[0])
    else:
        reads = []
    return (idx, reads)
