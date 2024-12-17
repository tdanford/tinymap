import json 
import gzip 
import pandas as pd

from dataclasses import dataclass, field
from typing import Iterable, List, Tuple, Generator, Dict
from pathlib import Path

from .fasta import kmers, open_with_gz
from .chains import Anchor, Chain


def argmin_kmer(window: List[str]) -> int: 
    """Takes a list of kmers, and returns the index of the minimum kmer within the list

    Arbitrarily returns an index if there is tie.
    """
    i = 0 
    for j in range(1, len(window)): 
        if window[j] < window[i]: 
            i = j 
    return i 

def minimizers(w: int, k: int, itrkmers: Iterable[kmers]) -> Generator[Tuple[str, int], None, None]: 
    """Generates a stream of minimizer/index tuples for a given sequence 

    A 'minimizer' is a k-mer which is the minimum (by lexicographic order) in a window of w successive k-mers. 

    This method generates a tuples, (s, i) where each s is a minimizer, and i is the offset of the start 
    of the minimizer s in the 'seq' string argument.  
    """
    window_offset = 0
    min_kmer = 0
    window = [] 
    for kmer in itrkmers: 
        if len(window) < w-1: 
            window.append(kmer) 
        else: 
            window.append(kmer) 
            while len(window) > w: 
                del(window[0])
            min_kmer -= 1 
            if min_kmer < 0: 
                min_kmer = argmin_kmer(window) 
                yield (window[min_kmer], window_offset+min_kmer)
            elif kmer < window[min_kmer]: 
                min_kmer = w-1
                yield (kmer, window_offset + w - 1)
            window_offset += 1

def create_seeds_df(w: int, k: int, seq: str) -> pd.DataFrame: 
    seed_dict = { 'kmer': [], 'offset': [] }
    for (kmer, offset) in minimizers(w, k, kmers(k, seq)): 
        seed_dict['kmer'].append(kmer) 
        seed_dict['offset'].append(offset) 
    return pd.DataFrame.from_dict(seed_dict)

def create_seeds(w: int, k: int, seq: str) -> Dict[str, List[int]]: 
    """Creates a dictionary whose keys are the minimizers of a sequence 'seq', and whose values 
    are a list of all the offsets in the 'seq' string where that key was found to be a minimizer.
    """
    seeds: Dict[str, List[int]] = {} 
    for (k, i) in minimizers(w, k, seq): 
        if k in seeds: 
            seeds[k].append(i) 
        else: 
            seeds[k] = [i] 
    return seeds


def binary_search(offsets: List[int], value: int) -> Tuple[int, int]: 
    if len(offsets) == 0: return (0, 0)
    if len(offsets) == 1: return (0, 1)
    left = 0
    right = len(offsets) 
    while right - left > 1: 
        middle = left + int((right - left) / 2.0)
        if offsets[middle] <= value: 
            left = middle 
        else: 
            right = middle 
    return (left, right) 

def filter_offsets(offsets: List[int], left: int, right: int) -> List[int]: 
    (lower, _) = binary_search(offsets, left) 
    (_, upper) = binary_search(offsets, right) 
    return offsets[lower:upper]
