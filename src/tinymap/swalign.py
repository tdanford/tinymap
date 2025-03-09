import logging
from enum import Enum
from dataclasses import dataclass, field
from typing import Tuple, List
from math import ceil, log10


@dataclass
class AlignmentParams:

    mismatch: float
    gap_open: float
    gap_extend: float
    match: float = field(default=0.0)


class AlignOp(Enum):
    ORIGIN = 1
    MATCH = 2
    MISMATCH = 3
    LEFT_OPEN = 4
    LEFT_EXTEND = 5
    RIGHT_OPEN = 8
    RIGHT_EXTEND = 9

    def is_gap(self) -> bool:
        return (self.value & 12) != 0

    def is_right_gap(self) -> bool:
        return self.value == 8 or self.value == 9

    def is_left_gap(self) -> bool:
        return self.value == 4 or self.value == 5

    def left_advance(self) -> int:
        return 1 if self.advances_left() else 0

    def right_advance(self) -> int:
        return 1 if self.advances_right() else 0

    def advances_left(self) -> bool:
        return (self.value & 10) != 0

    def advances_right(self) -> bool:
        return (self.value & 6) != 0

    def arrow(self) -> str:
        (left, right) = (self.left_advance(), self.right_advance())
        if left > 0 and right > 0:
            if self == AlignOp.MATCH:
                return LEFT_DOWN_DOUBLE_ARROW
            else:
                return LEFT_DOWN_ARROW
        elif left > 0:
            return RIGHT_ARROW
        elif right > 0:
            return LEFT_ARROW
        else:
            return "\u21bb"

    def __lt__(self, other: "AlignOp") -> bool:
        return self.value < other.value

    def prev_idx(self, idx: Tuple[int, int]) -> Tuple[int, int]:
        (i, j) = idx
        return (i - self.left_advance(), j - self.right_advance())

    def next_score(
        self, params: AlignmentParams, prev_score: float, new_gap_length: int
    ) -> float:
        if self == AlignOp.MATCH:
            return prev_score + params.match
        elif self == AlignOp.MISMATCH:
            return prev_score + params.mismatch
        elif self == AlignOp.LEFT_OPEN or self == AlignOp.RIGHT_OPEN:
            return prev_score + params.gap_open
        elif self == AlignOp.LEFT_EXTEND or self == AlignOp.RIGHT_EXTEND:
            return (
                prev_score + params.gap_open + (new_gap_length - 1) * params.gap_extend
            )
        else:
            raise ValueError()


LEFT_ARROW = "\u2190"
RIGHT_ARROW = "\u2193"
LEFT_DOWN_ARROW = "\u2199"
LEFT_DOWN_DOUBLE_ARROW = "\u21d9"


class Aligner:

    params: AlignmentParams

    M: int
    left_seq: str

    N: int
    right_seq: str

    scores: List[List[float]]
    ops: List[List[AlignOp]]
    gaps: List[List[int]]
    cigar: List[AlignOp]

    def __init__(self, params: AlignmentParams, left: str, right: str):
        self.M = len(left) + 1
        self.left_seq = left
        self.N = len(right) + 1
        self.right_seq = right
        self.params = params
        self.scores = [[0.0 for j in range(self.N)] for i in range(self.M)]
        self.ops = [[AlignOp.MISMATCH for j in range(self.N)] for i in range(self.M)]
        self.cigar = []
        self.gaps = []

    @property
    def score(self) -> float:
        return self.scores[self.M - 1][self.N - 1]

    @property
    def gapped_left_string(self) -> str:
        letters = []
        i = -1
        for op in self.cigar:
            if op is None:
                letters.append(" ")
            elif op.advances_left():
                i += 1
                letters.append(self.left_seq[i])
            else:
                letters.append("-")
        return "".join(letters)

    @property
    def gapped_right_string(self) -> str:
        letters = []
        i = -1
        for op in self.cigar:
            if op is None or op.advances_right():
                i += 1
                letters.append(self.right_seq[i])
            else:
                letters.append("-")
        return "".join(letters)

    @property
    def alignment(self) -> str:
        lstr = self.gapped_left_string
        rstr = self.gapped_right_string
        n = len(lstr)
        coord_rows = int(ceil(log10(max(1, n))))
        rows = [rstr, lstr]
        d = 1
        for i in range(coord_rows):
            offset = 0
            row = [" " for _ in range(n)]
            for k in range(0, n, d):
                row[k] = str(offset)
                offset = (offset + 1) % 10
            d *= 10
            rows.append("".join(row))
        return "\n".join(reversed(rows))

    @property
    def score_matrix(self) -> str:
        lines = list(
            reversed([" ".join(f"{s:0.1f}" for s in row) for row in self.scores])
        )
        return "\n".join(lines)

    @property
    def gap_matrix(self) -> str:
        lines = list(reversed([" ".join(f"{s:02}" for s in row) for row in self.gaps]))
        return "\n".join(lines)

    @property
    def backtrack_matrix(self) -> str:
        lines = list(reversed(["".join(op.arrow() for op in row) for row in self.ops]))
        return "\n".join(lines)

    def align(self, left_only: bool = False) -> float:
        logger = logging.getLogger("Aligner")
        self.ops[0][0] = AlignOp.ORIGIN
        self.scores[0][0] = 0.0
        gap_lengths = [[0 for j in range(self.N)] for i in range(self.M)]

        if self.M > 1:
            self.ops[1][0] = AlignOp.RIGHT_OPEN
            self.scores[1][0] = AlignOp.RIGHT_OPEN.next_score(
                self.params, self.scores[0][0], 1
            )
            gap_lengths[1][0] = 1

        if self.N > 1:
            self.ops[0][1] = AlignOp.LEFT_OPEN
            self.scores[0][1] = AlignOp.LEFT_OPEN.next_score(
                self.params, self.scores[0][0], 1
            )
            gap_lengths[0][1] = 1

        for i in range(2, self.M):
            self.ops[i][0] = AlignOp.RIGHT_EXTEND
            self.scores[i][0] = AlignOp.RIGHT_EXTEND.next_score(
                self.params, self.scores[0][0], i
            )
            gap_lengths[i][0] = i

        for i in range(2, self.N):
            self.ops[0][i] = AlignOp.LEFT_EXTEND
            self.scores[0][i] = AlignOp.LEFT_EXTEND.next_score(
                self.params, self.scores[0][0], i
            )
            gap_lengths[0][i] = i

        for i in range(1, self.M):
            for j in range(1, self.N):
                left_gap_op = (
                    AlignOp.RIGHT_EXTEND
                    if self.ops[i - 1][j].is_right_gap()
                    else AlignOp.RIGHT_OPEN
                )
                left_gap_length = (
                    gap_lengths[i - 1][j] + 1
                    if left_gap_op == AlignOp.RIGHT_EXTEND
                    else 1
                )
                prev_left_score = self.scores[i - left_gap_length][j]
                left_gap_score = left_gap_op.next_score(
                    self.params, prev_left_score, left_gap_length
                )

                right_gap_op = (
                    AlignOp.LEFT_EXTEND
                    if self.ops[i][j - 1].is_left_gap()
                    else AlignOp.LEFT_OPEN
                )
                right_gap_length = (
                    gap_lengths[i][j - 1] + 1
                    if right_gap_op == AlignOp.LEFT_EXTEND
                    else 1
                )
                prev_right_score = self.scores[i][j - right_gap_length]
                right_gap_score = right_gap_op.next_score(
                    self.params, prev_right_score, right_gap_length
                )

                mm_op = (
                    AlignOp.MATCH
                    if self.left_seq[i - 1] == self.right_seq[j - 1]
                    else AlignOp.MISMATCH
                )
                mm_score = mm_op.next_score(self.params, self.scores[i - 1][j - 1], 0)

                min_op_scores = sorted(
                    [
                        (right_gap_score, right_gap_op, right_gap_length),
                        (left_gap_score, left_gap_op, left_gap_length),
                        (mm_score, mm_op, 0),
                    ]
                )
                logger.debug(f"{i=} {j=} {min_op_scores=}")
                (next_score, next_op, next_gap_length) = min_op_scores[0]
                self.ops[i][j] = next_op
                self.scores[i][j] = next_score
                gap_lengths[i][j] = next_gap_length

        self.gaps = gap_lengths

        if left_only:
            min_j = 0
            top_i = self.M - 1
            for j in range(1, self.N):
                if self.scores[top_i][j] < self.scores[top_i][min_j]:
                    min_j = j
            (i, j) = (top_i, min_j)
        else:
            (i, j) = (self.M - 1, self.N - 1)

        self.cigar = self.backtrack_from(i, j)

    def backtrack_from(self, i: int, j: int) -> List[AlignOp]:
        first_j = j
        backtrack = []

        while i > 0 or j > 0:
            op = self.ops[i][j]
            backtrack.append(op)
            (i, j) = op.prev_idx((i, j))

        core_cigar = list(reversed(backtrack))
        return core_cigar + [None for _ in range(first_j + 1, self.N)]
