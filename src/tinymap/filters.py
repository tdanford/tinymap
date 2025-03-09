import hashlib
import random
from typing import List, Optional
from bitarray import bitarray


class SaltedHashFunction:

    salt: bytes

    def __init__(self, salt: bytes):
        self.salt = salt

    def hash(self, raw: bytes) -> bytes:
        m = hashlib.sha256()
        m.update(self.salt)
        m.update(raw)
        return m.digest()


def random_seed(nbytes: int = 8) -> bytes:
    random.randbytes(nbytes)


def random_seeds(k: int, nbytes: int = 8) -> List[bytes]:
    return [random_seed(nbytes) for _ in range(k)]


class HashingProber:

    k: int
    n: int
    seeds: List[bytes]
    functions: List[SaltedHashFunction]

    def __init__(self, n: int, k: int, seeds: List[bytes]):
        self.n = n
        self.k = k
        self.seeds = list(seeds)
        self.functions = [SaltedHashFunction(s) for s in self.seeds]

    def probe(self, i: int, bs: bytes) -> int:
        f = self.functions[i]
        idx = int.from_bytes(f.hash(bs)) % self.n
        return idx


class CountingFilter(HashingProber):

    array: List[int]

    def __init__(
        self, n: int, k: int, seeds: List[bytes], array: Optional[List[int]] = None
    ):
        HashingProber.__init__(self, n, k, seeds)
        self.array = array or [0 for i in range(n)]

    def add_string(self, s: str):
        self.add_bytes(s.encode("utf-8"))

    def add_bytes(self, bs: bytes):
        for i in range(len(self.functions)):
            idx = self.probe(i, bs)
            self.array[idx] += 1

    def count_string(self, s: str) -> int:
        return self.count_bytes(s.encode("utf-8"))

    def count_bytes(self, bs: bytes) -> int:
        counts = [self.array[self.probe(i, bs)] for i in range(self.k)]
        return min(counts)


class BloomFilter(HashingProber):

    array: bitarray

    def __init__(
        self, n: int, k: int, seeds: List[bytes], array: Optional[bitarray] = None
    ):
        HashingProber.__init__(self, n, k, seeds)
        self.array = array or bitarray(n)

    def add_string(self, s: str):
        self.add_bytes(s.encode("utf-8"))

    def add_bytes(self, bs: bytes):
        for i in range(len(self.functions)):
            idx = self.probe(i, bs)
            self.array[idx] = 1

    def contains_string(self, s: str) -> bool:
        return self.contains_bytes(s.encode("utf-8"))

    def contains_bytes(self, bs: bytes) -> bool:
        for i in range(len(self.functions)):
            idx = self.probe(i, bs)
            if self.array[idx] == 0:
                return False
        return True
