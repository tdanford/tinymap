from dataclasses import dataclass, field
from typing import List, Dict
from math import log2, log
from heapq import heapify, heappush
from functools import reduce
import logging

def lambda_c(l: int, w: int) -> float: 
    if l <= 0: return 0.0 
    else: return 0.1 * w * l + 0.5 * log2(l)

INF = float('inf')
NEG_INF = float('-inf')

@dataclass
class Anchor: 
    """An Anchor matches a region on the reference sequence to a region on the query (read) sequence. 

    From the minimap2 paper: "An anchor is a 3-tuple (x, y, w), indicating interval [x- w +1, x] on 
    the reference matching interval [y-w+1, y] on the query"

    Unlike in the minimap2 paper, we index the relative locations of the Anchor by the starting point 
    of the match, not the ending -- this has some subtle differences with some of the formulae, later.
    """
    read_offset: int 
    reference_offset: int
    k: int = field(compare=False)

    def __lt__(self, other: 'Anchor') -> bool: 
        return ( 
            self.reference_offset < other.reference_offset or 
            (self.reference_offset == other.reference_offset and self.read_offset < other.read_offset)
        )
    
    @property 
    def read_end(self) -> int: return self.read_offset + self.k

    @property 
    def reference_end(self) -> int: return self.reference_offset + self.k

    def read_region(self) -> 'Region': 
        return Region(self.read_offset, self.read_end) 
    
    def ref_region(self) -> 'Region': 
        return Region(self.reference_offset, self.reference_end)

    def read_diff(self, other: 'Anchor') -> int: 
        return self.read_offset - other.read_offset 

    def reference_diff(self, other: 'Anchor') -> int: 
        return self.reference_offset - other.reference_offset
    
    def match_score(self, w: int, other: 'Anchor') -> float: 
        return min(other.read_diff(self), other.reference_diff(self), w)
    
    def gap_score(self, w: int, max_distance: int, other: 'Anchor') -> float | None: 
        # self is j, other is i 
        if other.read_offset <= self.read_offset or max(other.reference_diff(self), other.read_diff(self)) > max_distance: 
            return INF
        else: 
            return lambda_c(other.read_diff(self) - other.reference_diff(self), w)

@dataclass 
class Region: 
    """A region is an interval [left, right) on a sequence.  

    Here left <= right (where equality means 'a zero-width interval') and we note the *exclusive* nature 
    of the right bound.
    """
    left: int 
    right: int 

    def extract(self, v: str) -> str: return v[self.left:self.right]

    def contains(self, pt: int) -> bool: return self.left <= pt and pt < self.right

    def overlaps(self, other: 'Region') -> bool: return self.contains(other.left) or other.contains(self.left) 

    @property 
    def width(self) -> int: return self.right - self.left 

    def __len__(self) -> int: return self.width 

    def count_overlap(self, other: 'Region') -> int: 
        max_left = max(self.left, other.left) 
        min_right = min(self.right, other.right) 
        return min_right - max_left

@dataclass
class Chain: 
    """A Chain is a linear sequence of Anchors between a query (read) and reference sequence.
    """

    score: float
    hits: List[Anchor] = field(compare=False)
    secondary: List['Chain'] = field(default_factory=list, compare=False, repr=False)

    def __getitem__(self, key): 
        if isinstance(key, slice): 
            return Chain(self.score, self.hits[key], [])
        else: 
            return self.hits[key]

    def __len__(self) -> int: return len(self.hits) 

    def __lt__(self, other: 'Chain') -> bool: 
        return self.score < other.score

    def add_secondary_chain(self, chain: 'Chain'): 
        heappush(self.secondary, chain)
    
    def chunks(self) -> List[List[int]]: 
        def combiner(chunks: List[List[int]], i: int) -> List[List[int]]: 
            if len(chunks) == 0: 
                return [[i]]
            elif self.hits[chunks[-1][-1]].read_end >= self.hits[i].read_offset: 
                chunks[-1].append(i) 
                return chunks 
            else: 
                return chunks + [[i]]
        return reduce(
            combiner, 
            range(len(self.hits)), 
            []
        )
    
    @property 
    def best_secondary_score(self) -> float: 
        if len(self.secondary) == 0: 
            return 0.0 
        else: 
            return max(self.secondary).score 
    
    @property 
    def mapq(self) -> float: 
        # mapq = 40 * (1 - f2 / f1) * min{ 1, m / 10 } * log f1
        # where m = # of anchors
        # f1 is this score, and f2 the best secondary score
        # log is natural logarithm 
        m = len(self.hits) 
        f21 = self.best_secondary_score / self.score if self.score != 0.0 else 1.0
        return 40 * (1 - f21) * min(1, m / 10) * log(self.score) 

    @property 
    def reference_start(self) -> int: 
        return self.hits[0].reference_offset
    
    @property
    def reference_end(self) -> int: 
        return self.hits[-1].reference_end

    @property 
    def read_start(self) -> int: 
        return self.hits[0].read_offset
    
    @property
    def read_end(self) -> int: 
        return self.hits[-1].read_end

    @property 
    def reference_region(self) -> Region:
        return Region(left=self.reference_start, right=self.reference_end)

    @property 
    def read_region(self) -> Region:
        return Region(left=self.read_start, right=self.read_end)
    
    @property
    def reference_width(self) -> int: 
        return self.reference_end - self.reference_start 
    
    def contains_reference(self, i: int) -> bool: 
        return self.reference_start <= i and i < self.reference_end
    
    def overlaps(self, other: 'Chain') -> bool: 
        return ( 
            self.contains_reference(other.reference_start) or 
            other.contains_reference(self.reference_start)
        )
    
    def overlap_pct(self, other: 'Chain'): 
        srr: Region = self.reference_region 
        orr: Region = other.reference_region 
        if self.reference_width < other.reference_width: 
            return (srr.count_overlap(orr)) / max(1, srr.width)
        else: 
            return (srr.count_overlap(orr)) / max(1, orr.width)



@dataclass
class ChainParams: 
    w : int
    max_distance: int = field(default=1000)
    max_iter: int = field(default=50)

class Chainer: 

    params: ChainParams
    logger: logging.Logger
    
    hits: List[Anchor]
    scores: List[float] 
    backtrack: List[int] 
    chains: List[Chain]

    def __init__(self, params: ChainParams, hits: List[Anchor]): 
        self.params = params
        self.hits = hits 
        self.scores = [0.0 for i in range(len(self.hits))] 
        self.backtrack = [-1 for i in range(len(self.hits))]
        self.hit_chain = None
        self.logger = logging.getLogger("Chainer")
    
    def compute_score(self, j: int, i: int) -> float: 
        assert self.hits[j] < self.hits[i], f"{j=} {self.hits[j]} and {i=} {self.hits[i]} were specified in the wrong order"
        msc = self.hits[j].match_score(self.params.w, self.hits[i])
        gsc = self.hits[j].gap_score(self.params.w, self.params.max_distance, self.hits[i])
        sc = msc - gsc 

        assert self.scores[j] is not None, f"Score for hit {j=} was None, {self.scores=}"
        return self.scores[j] + sc 

    def recover_chains(self): 
        self.logger.log(logging.DEBUG, f"recover_chains()")

        chain_membership = [None for _ in range(len(self.hits))]
        total_chains = []

        chain_idx = 0
        for tail_i in range(len(self.hits)-1, -1, -1): 
            if chain_membership[tail_i] is None: 
                c = [] 
                ii = tail_i
        
                while ii >= 0: 
                    c.append(self.hits[ii]) 
                    chain_membership[ii] = chain_idx
                    ii = self.backtrack[ii] 
                
                chain_idx += 1
                hit_chain = Chain(score=self.scores[tail_i], hits=list(reversed(c)))
                total_chains.append(hit_chain)

        self.cluster_chains(total_chains)
    
    def cluster_chains(self, raw_chains: List[Chain]): 
        clusters = [] 
        raw_chains.sort(reverse=True)  # sort in descending 'chain score' order 

        for chain in raw_chains: 
            primary_chain: Chain = None 
            for c in clusters: 
                if chain.overlap_pct(c) >= 0.5: 
                    primary_chain = c 
                    break 
            if primary_chain is None: 
                clusters.append(chain) 
            else: 
                primary_chain.add_secondary_chain(chain)

        self.chains = clusters 

    def chain_align(self): 
        n = len(self.hits)
        self.backtrack = [-1 for _ in range(n)]
        self.scores = [0.0 for _ in range(n)]

        for i in range(n): 
            self.logger.log(logging.DEBUG, f"Hit {i=} {self.hits[i]}")
            max_score = 0.0
            max_j = None

            distant = i-1
            while distant >= 0 and self.hits[i].reference_diff(self.hits[distant]) <= self.params.max_distance: 
                distant -= 1
            st = max(0, i-self.params.max_iter, distant)
            self.logger.log(logging.DEBUG, f"{st=}")

            for j in range(i-1, st-1, -1): 
                score = self.compute_score(j, i) 
                self.logger.log(logging.DEBUG, f"{j=} {score=}")
                if score > NEG_INF and (max_j is None or score > max_score):
                    self.logger.log(logging.DEBUG, f"new MAX {j=} {score=} -> {max_score=} {max_j=}")
                    max_score = score 
                    max_j = j 
                else: 
                    self.logger.log(logging.DEBUG, f"Iter {j=} {score=} vs. {max_j=} {max_score=}")

            if max_j is not None: 
                self.logger.log(logging.DEBUG, f"MAX found {max_j=} {max_score=}")
                self.scores[i] = max_score 
                self.backtrack[i] = max_j
        
        self.recover_chains()
        
