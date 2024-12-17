from dataclasses import dataclass, field 
from typing import Dict, List, Callable, Generator, Tuple, Iterable
from pathlib import Path 
from contextlib import contextmanager
from typing import TextIO
from itertools import groupby
import logging

import csv 
import gzip 
import random 
import re

def seq_filter(*seqs: str) -> Callable[[str], bool]: 
    if len(seqs) == 0: return lambda s: True
    def accepter(name: str) -> bool:
        for s in seqs: 
            if name == s: 
                return True
        return False
    return accepter

REV = {
    'A': 'T', 
    'T': 'A', 
    'G': 'C', 
    'C': 'G', 
    'N': 'N',
    'a': 't', 
    't': 'a', 
    'g': 'c', 
    'c': 'g', 
    'n': 'n'
}

def revcomp(seq: str) -> str: 
    return ''.join(REV[seq[i]] for i in range(len(seq)-1, -1, -1))

@dataclass
class Read: 
    name: str 
    seq: str 
    qual: str 

FASTQ_SEQ_HEADER = re.compile("@([^ ]+) (.*)")
FASTQ_QUAL_HEADER = re.compile("\\+(.*)")

@contextmanager
def open_with_gz(p: Path): 
    inf = None
    try: 
        if is_gzip(p): 
            inf = gzip.open(p, 'rt') 
        else: 
            inf = p.open('rt') 
        yield inf 
    finally: 
        if inf is not None: 
            inf.close()

def find_read_files(dir: Path) -> Generator[Path, None, None]: 
    if dir.is_dir(): 
        for p in dir.iterdir(): 
            if p.is_file() and is_fastq(p): 
                yield p 

def find_ref_files(dir: Path) -> Generator[Path, None, None]: 
    if dir.is_dir(): 
        for p in dir.iterdir(): 
            if p.is_file() and is_fasta(p): 
                yield p 

def is_fasta(p: Path) -> bool: 
    nn = p.name.lower()
    ff = nn.find('.fasta') != -1 or nn.find('.fa.') != -1 or nn.endswith('.fa') 
    return ff 

def is_fastq(p: Path) -> bool: 
    return p.name.lower().find(".fastq") != -1

def is_gzip(p: Path) -> bool: 
    return p.name.lower().endswith(".gz")

def parse_fastq_stream(inf: TextIO) -> Generator[Read, None, None]: 
    try: 
        while True: 
            seq_header = next(inf).strip()
            seq = next(inf).strip()
            qual_header = next(inf).strip() 
            qual = next(inf).strip()
            sheader_matcher = FASTQ_SEQ_HEADER.match(seq_header) 
            qheader_matcher = FASTQ_QUAL_HEADER.match(qual_header) 
            if not sheader_matcher or not qheader_matcher: 
                raise ValueError(f"{seq_header} didn't match format")
            name = sheader_matcher.groups()[0] 
            yield Read(name, seq, qual)
    except StopIteration: 
        ...

def index_fasta_files(index_file: Path, *files: Path): 
    """Uses fasta_seq_offsets to generate sequence offsets in a fasta file, and then writes it to a file
    """
    with index_file.open('wt') as outf: 
        writer = csv.writer(outf, dialect=csv.excel_tab) 
        writer.writerow(["filename", "seq_name", "offset"])
        for fasta in files: 
            for (seq_name, offset) in fasta_seq_offsets(fasta): 
                writer.writerow([fasta.name, seq_name, offset])
                print(f"{seq_name=}, {offset=}")

@dataclass
class FastaIndex: 

    offsets: Dict[str, List[Tuple[Path, int]]] = field(repr=False)

    def load_seq(self, name: str) -> str: 
        logger = logging.getLogger("FastaIndex")
        logger.info(f"load_seq({name}) = {self.offsets.get(name, [])}")
        name_seq_tuple: List[Tuple[str, str]] = [
            next(read_fasta_stream(file, offset=offset)) for (file, offset) in self.offsets.get(name, [])
        ]
        if len(name_seq_tuple) < 1: 
            raise ValueError(f"Couldn't find seq {name} in FastaIndex {self.offsets.keys()}")
        return name_seq_tuple[0][1]

def load_fasta_index(index_file: Path) -> FastaIndex:
    logger = logging.getLogger("fasta.py") 
    logger.info(f"Loading FASTA index from {index_file.name}")
    with index_file.open('rt') as inf: 
        reader = csv.reader(inf, dialect=csv.excel_tab) 
        _ = next(reader) # header
        triples = list(reader) 
        seqs = groupby(triples, lambda t: t[1]) 
        return FastaIndex(offsets={
            seq_name: [(Path(filename), int(offset)) for (filename, _, offset) in triples]
            for (seq_name, triples) in seqs
        })

def fasta_seq_offsets(p: Path) -> Generator[Tuple[str, int], None, None]: 
    """Generates (seq_name, offset) tuples for distinct sequences in a FASTA file 
    """
    with open_with_gz(p) as inf: 
        try: 
            while True: 
                offset = inf.tell() 
                line = next(inf).rstrip('\n') 
                if line.startswith('>'): 
                    name = line[1:]
                    yield (name, offset) 
        except StopIteration: 
            ...

def read_fastq(p: Path) -> Generator[Read, None, None]: 
    with open_with_gz(p) as inf: 
        yield from parse_fastq_stream(inf) 

def parse_fasta_stream(inf: TextIO, offset: int = None) -> Generator[Tuple[str, str], None, None]: 
    def process_string(seq: str) -> str: 
        return seq.upper()
    current_name = None 
    current_sequence = []
    if offset is not None: 
        inf.seek(offset)
    try: 
        while True:
            line = next(inf).strip()
            if len(line) > 0: 
                if line.startswith('>'): 
                    if current_name is not None: 
                        yield (current_name, process_string(''.join(current_sequence)))
                    name = line[1:] 
                    current_name = name 
                    current_sequence = [] 
                else: 
                    if current_name is not None: 
                        current_sequence.append(line)
    except StopIteration: 
        ...
    if current_name is not None: 
        yield (current_name, process_string(''.join(current_sequence)))

def read_fasta_stream(p: Path, offset: int = None) -> Generator[Tuple[str, str], None, None]: 
    with open_with_gz(p) as inf: 
        yield from parse_fasta_stream(inf, offset=offset)

def read_fasta(p: Path, offset: int = None) -> Dict[str, str]: 
    with open_with_gz(p) as inf: 
        return {
            name: seq for (name, seq) in parse_fasta_stream(inf, offset=offset)
        }
        
def kmers(k: int, seq: str) -> Generator[str, None, None]: 
    for i in range(len(seq) - k + 1): 
        yield seq[i:i+k]

def revcomp_kmers(k: int, seq: str) -> Generator[str, None, None]: 
    for i in range(len(seq) - k, -1, -1): 
        yield revcomp(seq[i:i+k])

def random_read(seq: str, read_length: int) -> Tuple[str, int]: 
    offset = random.randint(0, len(seq) - read_length) 
    return seq[offset:offset+read_length], offset
