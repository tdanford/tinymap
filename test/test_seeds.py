from tinymap.seeds import *


def test_minimizers():
    seq = "GTCATGC"
    w = 3
    k = 4
    kitr = kmers(k, seq)
    mins = list(minimizers(w, k, kitr))
    assert mins == [("CATG", 2), ("ATGC", 3)]
