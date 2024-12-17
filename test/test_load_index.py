from pytest import mark

from tinymap.test_align import load_index

@mark.skip(reason="Just here for stepping into and debugging load_index")
def test_load(): 
    (idx, reads) = load_index()
