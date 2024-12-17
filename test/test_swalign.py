
from tinymap.swalign import AlignOp, AlignmentParams, Aligner

def test_align_ops_right_advance(): 
    assert AlignOp.ORIGIN.right_advance() == 0 
    assert AlignOp.MATCH.right_advance() == 1
    assert AlignOp.MISMATCH.right_advance() == 1
    assert AlignOp.LEFT_OPEN.right_advance() == 1
    assert AlignOp.LEFT_EXTEND.right_advance() == 1
    assert AlignOp.RIGHT_OPEN.right_advance() == 0
    assert AlignOp.RIGHT_EXTEND.right_advance() == 0

def test_align_ops_left_advance(): 
    assert AlignOp.ORIGIN.left_advance() == 0 
    assert AlignOp.MATCH.left_advance() == 1
    assert AlignOp.MISMATCH.left_advance() == 1
    assert AlignOp.LEFT_OPEN.left_advance() == 0
    assert AlignOp.LEFT_EXTEND.left_advance() == 0
    assert AlignOp.RIGHT_OPEN.left_advance() == 1
    assert AlignOp.RIGHT_EXTEND.left_advance() == 1

def test_align_ops_is_gap(): 
    assert not AlignOp.ORIGIN.is_gap()
    assert not AlignOp.MATCH.is_gap()
    assert not AlignOp.MISMATCH.is_gap()
    assert AlignOp.LEFT_OPEN.is_gap()
    assert AlignOp.LEFT_EXTEND.is_gap()
    assert AlignOp.RIGHT_OPEN.is_gap()
    assert AlignOp.RIGHT_EXTEND.is_gap()

def test_align_ops_prev_score(): 
    p = AlignmentParams(mismatch=3.0, gap_open=5.0, gap_extend=7.0)
    assert AlignOp.MATCH.next_score(p, 1.0, 0) == 1.0 
    assert AlignOp.MISMATCH.next_score(p, 1.0, 0) == 4.0
    assert AlignOp.LEFT_OPEN.next_score(p, 1.0, 1) == 6.0
    assert AlignOp.LEFT_EXTEND.next_score(p, 1.0, 2) == 13.0
    assert AlignOp.RIGHT_OPEN.next_score(p, 1.0, 1) == 6.0
    assert AlignOp.RIGHT_EXTEND.next_score(p, 1.0, 2) == 13.0

def test_no_gap_alignment(): 
    left = "tim" 
    right = "tim"
    p = AlignmentParams(mismatch=1.0, gap_open=1.25, gap_extend=0.5)
    aligner = Aligner(params=p, left=left, right=right)
    aligner.align()
    assert aligner.gapped_left_string == left
    assert aligner.gapped_right_string == right 

def test_gapped_alignment(): 
    left = "tim" 
    right = "tifoom"
    p = AlignmentParams(mismatch=1.0, gap_open=1.25, gap_extend=0.5)
    aligner = Aligner(params=p, left=left, right=right)
    aligner.align()
    assert aligner.gapped_left_string == "ti---m"
    assert aligner.gapped_right_string == "tifoom" 
    assert aligner.cigar == [AlignOp.MATCH, AlignOp.MATCH, AlignOp.LEFT_OPEN, AlignOp.LEFT_EXTEND, AlignOp.LEFT_EXTEND, AlignOp.MATCH]

def test_one_string_alignment(): 
    left = ""
    right = "abc"
    p = AlignmentParams(mismatch=1.0, gap_open=1.25, gap_extend=0.5)
    aligner = Aligner(params=p, left=left, right=right)
    aligner.align()
    assert aligner.gapped_left_string == '---'
    assert aligner.gapped_right_string == 'abc'

    aligner = Aligner(params=p, left=right, right=left)
    aligner.align()
    assert aligner.gapped_left_string == 'abc'
    assert aligner.gapped_right_string == '---'

def test_full_left_alignment(): 
    p = AlignmentParams(3.0, 2.0, 0.5)

    left = "ccc"
    right = "cccatgatgctc" 
    aligner = Aligner(params=p, left=left, right=right)
    aligner.align()
    assert aligner.gapped_left_string == 'ccc---------'
    assert aligner.gapped_right_string == 'cccatgatgctc'

def test_sw_left_alignment(): 
    p = AlignmentParams(3.0, 2.0, 0.5)
    left = "GAA"
    right = "GAAATTTTAGTTGTCGTAGTAGG"
    aligner = Aligner(params=p, left=left, right=right)
    aligner.align(left_only=True)

    assert aligner.gapped_left_string ==  'GAA                    '
    assert aligner.gapped_right_string == 'GAAATTTTAGTTGTCGTAGTAGG'



    