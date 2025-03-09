import hashlib
import random
import math
from typing import Dict, List, Optional, Tuple
from bitarray import bitarray
import zlib


class SaltedHashFunction:

    salt: bytes

    def __init__(self, salt: bytes):
        self.salt = salt

    def hash(self, raw: bytes) -> bytes:
        m = hashlib.sha256()
        m.update(self.salt)
        m.update(raw)
        return m.digest()


def find_optimal_k(nbits: int, nelmts: int) -> int:
    f = nbits / nelmts
    return math.ceil(f * math.log(2.0))


def find_optimal_bits(nelmts: int, err: float) -> int:
    l2 = math.log(2.0)
    p = nelmts * math.log(err)
    return math.ceil(-p / (l2 * l2))


def optimal_filter_args(nelmts: int, err: float) -> Tuple[int, int, List[bytes]]:
    n = find_optimal_bits(nelmts, err)
    k = find_optimal_k(n, nelmts)
    seeds = random_seeds(k)
    return (n, k, seeds)


def random_seed(nbytes: int = 8) -> bytes:
    return random.randbytes(nbytes)


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

    def is_compatible(self, prober: "HashingProber") -> bool:
        return self.k == prober.k and self.n == prober.n and self.seeds == prober.seeds

    def tojson(self) -> Dict:
        return {"n": self.n, "k": self.k, "seeds": [b.hex() for b in self.seeds]}


class CountingFilter(HashingProber):

    array: List[int]

    def __init__(
        self, n: int, k: int, seeds: List[bytes], array: Optional[List[int]] = None
    ):
        HashingProber.__init__(self, n, k, seeds)
        self.array = array or [0 for i in range(n)]

    def __add__(self, other: "CountingFilter") -> "CountingFilter":
        if not self.is_compatible(other):
            raise ValueError(f"Can't add non-compatible filters")

        sums = [self.array[i] + other.array[i] for i in range(self.n)]
        return CountingFilter(self.n, self.k, self.seeds, sums)

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

    @staticmethod
    def load(d: Dict) -> "BloomFilter":
        n = d.get("n")
        k = d.get("k")
        seeds = [bytes.fromhex(h) for h in d.get("seeds")]
        array_compressed_hex = d.get("array")
        array = bitarray.frombytes(zlib.decompress(bytes.fromhex(array_compressed_hex)))
        return BloomFilter(n, k, seeds, array)

    array: bitarray

    def __init__(
        self, n: int, k: int, seeds: List[bytes], array: Optional[bitarray] = None
    ):
        HashingProber.__init__(self, n, k, seeds)
        self.array = array or bitarray(n)

    def tojson(self) -> Dict:
        d = super().tojson()
        d["array"] = zlib.compress(self.array.tobytes()).hex()
        return d

    def add_filter(self, other: "BloomFilter"):
        if not self.is_compatible(other):
            raise ValueError("Cannot union two non-compatible filters")
        self.array = self.array | other.array

    def union(self, other: "BloomFilter") -> "BloomFilter":
        if not self.is_compatible(other):
            raise ValueError("Cannot union two non-compatible filters")
        new_array = self.array | other.array
        return BloomFilter(self.n, self.k, self.seeds, new_array)

    def estimate_size(self) -> float:
        f1 = self.n / self.k
        f2 = self.array.count(1) / self.n
        return -f1 * math.log(1.0 - f2)

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
