
from hypothesis import given, strategies as st 
from tinymap.chains import Anchor

@given(st.integers(min_value=1))
def test_anchor_ordering(k):
    a1 = Anchor(10, 100, k) 
    a1a = Anchor(10, 100, k) 
    a2 = Anchor(10, 200, k) 
    a3 = Anchor(5, 100, k) 
    a4 = Anchor(5, 200, k) 

    assert a1 == a1a 
    assert a1 < a2 
    assert a3 < a1 
    assert a1 < a4

@given(st.integers(min_value=1))
def test_anchor_diffs(k): 
    a1 = Anchor(10, 100, k) 
    a2 = Anchor(20, 200, k) 

    assert a2.read_diff(a1) == 10 
    assert a2.reference_diff(a1) == 100 
    assert a1.read_diff(a1) == 0
    assert a1.reference_diff(a1) == 0
    assert a1.read_diff(a2) == -10 
    assert a1.reference_diff(a2) == -100