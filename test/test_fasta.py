
from pathlib import Path 
from tinymap.fasta import read_fasta, revcomp, seq_filter, kmers, revcomp_kmers

def test_kmers(): 
    assert list(kmers(3, 'aaattt')) == ['aaa', 'aat', 'att', 'ttt']

def test_revcomp_kmers(): 
    assert list(revcomp_kmers(3, 'aaaggg')) == ['ccc', 'cct', 'ctt', 'ttt']

def test_read_fasta(simple_fasta: Path): 
    strs = read_fasta(simple_fasta)
    assert len(strs) == 2 
    assert strs['seq1'] == 'ATGCA'
    assert strs['seq2 something'] == 'TTTAAACCCGGG'

def test_revcomp(): 
    assert revcomp("ATGC") == "GCAT"
    assert revcomp("") == ""

def test_seq_filter(): 
    f1 = seq_filter("abc") 
    assert f1("abc")
    assert not f1("cde") 
    assert not f1("") 
    assert not f1("ab")
    assert not f1("abc de")