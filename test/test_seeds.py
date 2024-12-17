
from tinymap.seeds import * 

def test_minimizers(): 
    seq = "GTCATGC" 
    w = 3 
    k = 4 
    mins = list(minimizers(w, k, seq))
    assert mins == [('CATG', 2), ('ATGC', 3)]