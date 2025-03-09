from tinymap.filters import *


def test_count_strings():

    n = 100
    k = 3
    seeds = [b"001", b"010", b"100"]

    s1 = "abcdefg"
    s2 = "tuvwxyz"
    s3 = "abcdef"
    s4 = "bcdefg"
    s5 = ""

    cf = CountingFilter(n, k, seeds)

    for i in range(3):
        cf.add_string(s1)

    assert cf.count_string(s1) == 3
    assert cf.count_string(s2) == 0
    assert cf.count_string(s3) == 0

    for i in range(2):
        cf.add_string(s2)

    assert cf.count_string(s1) == 3
    assert cf.count_string(s2) == 2
    assert cf.count_string(s3) == 0

    for i in range(1):
        cf.add_string(s3)

    assert cf.count_string(s1) == 3
    assert cf.count_string(s2) == 2
    assert cf.count_string(s3) == 1


def test_add_deterministic_strings():

    n = 100
    k = 3
    seeds = [b"001", b"010", b"100"]

    s1 = "abcdefg"
    s2 = "tuvwxyz"
    s3 = "abcdef"
    s4 = "bcdefg"
    s5 = ""

    bf = BloomFilter(n, k, seeds)

    assert not bf.contains_string(s1)
    assert not bf.contains_string(s2)
    assert not bf.contains_string(s3)
    assert not bf.contains_string(s4)
    assert not bf.contains_string(s5)

    bf.add_string(s1)

    assert bf.contains_string(s1)
    assert not bf.contains_string(s2)
    assert not bf.contains_string(s3)
    assert not bf.contains_string(s4)
    assert not bf.contains_string(s5)

    bf.add_string(s2)

    assert bf.contains_string(s1)
    assert bf.contains_string(s2)
    assert not bf.contains_string(s3)
    assert not bf.contains_string(s4)
    assert not bf.contains_string(s5)

    bf.add_string(s3)

    assert bf.contains_string(s1)
    assert bf.contains_string(s2)
    assert bf.contains_string(s3)
    assert not bf.contains_string(s4)
    assert not bf.contains_string(s5)

    bf.add_string(s4)

    assert bf.contains_string(s1)
    assert bf.contains_string(s2)
    assert bf.contains_string(s3)
    assert bf.contains_string(s4)
    assert not bf.contains_string(s5)

    bf.add_string(s5)

    assert bf.contains_string(s1)
    assert bf.contains_string(s2)
    assert bf.contains_string(s3)
    assert bf.contains_string(s4)
    assert bf.contains_string(s5)

    assert not bf.contains_string(s1 + s2)
    assert not bf.contains_string(s2 + s3)
    assert not bf.contains_string(s3 + s4)


def test_union():

    n = 100
    k = 3
    seeds = [b"001", b"010", b"100"]

    s1 = "abcdefg"
    s2 = "tuvwxyz"
    s3 = "abcdef"
    s4 = "bcdefg"

    b1 = BloomFilter(n, k, seeds)
    b2 = BloomFilter(n, k, seeds)

    b1.add_string(s1)
    b1.add_string(s2)
    b2.add_string(s3)
    b2.add_string(s4)

    assert b1.contains_string(s1)
    assert b1.contains_string(s2)
    assert not b1.contains_string(s3)
    assert not b1.contains_string(s4)

    assert int(b1.estimate_size()) == 2

    assert not b2.contains_string(s1)
    assert not b2.contains_string(s2)
    assert b2.contains_string(s3)
    assert b2.contains_string(s4)

    assert int(b2.estimate_size()) == 2

    b3 = b1.union(b2)
    assert b3.contains_string(s1)
    assert b3.contains_string(s2)
    assert b3.contains_string(s3)
    assert b3.contains_string(s4)
    assert int(b3.estimate_size()) == 4
