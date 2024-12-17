
from tinymap.fasta import kmers 

def test_simple_kmers(): 
    seq = "ATGCACG"
    assert list(kmers(4, seq)) == ['ATGC', 'TGCA', 'GCAC', 'CACG']