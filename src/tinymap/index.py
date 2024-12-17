
import json
import gzip 
import logging
import pandas as pd 

from dataclasses import dataclass, field 
from typing import Generator, List, Dict, Tuple, Iterable
from pathlib import Path
from rich.progress import track
from itertools import groupby 

from .fasta import FastaIndex, load_fasta_index, read_fasta, read_fastq, find_read_files, find_ref_files, kmers, Read
from .chains import Anchor, Chain, Chainer, ChainParams, Region
from .seeds import create_seeds, create_seeds_df, minimizers
from .swalign import Aligner, AlignmentParams

@dataclass 
class Alignment: 

    read_name: str 
    read_seq: str 
    ref_name: str
    chain: Chain 

    @property 
    def mapq(self) -> float: return self.chain.mapq

    @property 
    def chain_score(self) -> float: return self.chain.score

    def alignable_chunks(self, buffer: int = 20) -> List[Tuple[Region, Region]]: 
        """Break up the regions of the read into separately alignable chunks

        Each chunk is a pair of Regions that need to be aligned together; one from the read 
        and the second from the reference itself. 
        """
        chunks: List[List[int]] = self.chain.chunks() 
        read_regions = [] 
        ref_regions = [] 
        logger = logging.getLogger("Alignment") 
        logger.debug(f"{chunks=}")

        for i in range(len(chunks)): 
            left_anchor = self.chain[chunks[i-1][-1]] if i > 0 else None
            right_anchor = self.chain[chunks[i][0]]
            read_regions.append(self.read_region(left_anchor, right_anchor)) 
            ref_regions.append(self.reference_region(left_anchor, right_anchor, buffer=buffer))

            left_anchor: Anchor = self.chain[chunks[i][0]] 
            right_anchor: Anchor = self.chain[chunks[i][-1]]
            read_regions.append(Region(left_anchor.read_offset, right_anchor.read_end))
            ref_regions.append(Region(left_anchor.reference_offset, right_anchor.reference_end))
        
        left_anchor = self.chain[chunks[-1][-1]] 
        right_anchor = None
        read_regions.append(self.read_region(left_anchor, right_anchor)) 
        ref_regions.append(self.reference_region(left_anchor, right_anchor, buffer=buffer))

        logger.debug(f"{read_regions=}")
        logger.debug(f"{ref_regions=}")

        return [
            (rer, rfr) for (rer, rfr) in zip(read_regions, ref_regions) 
            if rer.width > 0 or rfr.width > 0
        ]

    def alignable_regions(self, buffer: int = 20) -> List[Tuple[Region, Region]]: 
        anchor_pairs = list(zip([None] + self.chain.hits, self.chain.hits + [None]))
        read_regions = [self.read_region(*r) for r in anchor_pairs]
        ref_regions = [self.reference_region(r[0], r[1], buffer=buffer) for r in anchor_pairs]
        return [
            (rer, rfr) for (rer, rfr) in zip(read_regions, ref_regions) 
            if rer.width > 0 and rfr.width > 0
        ]

    def read_region(self, left_anchor: Anchor, right_anchor: Anchor) -> Region: 
        if left_anchor is None and right_anchor is None: 
            raise ValueError(f"Both left and right anchors can't be None, when defining a region") 
        if left_anchor is None: 
            return Region(0, right_anchor.read_offset) 
        elif right_anchor is None: 
            return Region(left_anchor.read_end, len(self.read_seq)) 
        else: 
            return Region(left_anchor.read_end, right_anchor.read_offset) 
    
    def reference_region(self, left_anchor: Anchor, right_anchor: Anchor, buffer: int = 20) -> Region: 
        if left_anchor is None and right_anchor is None: 
            raise ValueError(f"Both left and right anchors can't be None, when defining a region") 
        if left_anchor is None: 
            left_edge = right_anchor.read_offset + buffer
            return Region(max(0, right_anchor.reference_offset - left_edge), right_anchor.reference_offset) 
        elif right_anchor is None: 
            right_edge = (len(self.read_seq) - left_anchor.read_end) + buffer
            return Region(left_anchor.reference_end, left_anchor.reference_end + right_edge) 
        else: 
            return Region(left_anchor.reference_end, right_anchor.reference_offset) 

    def __lt__(self, other: 'Alignment'): 
        return self.chain < other.chain

    def __repr__(self) -> str: 
        read_region = self.chain.read_region
        ref_region = self.chain.reference_region
        return f"Align({self.read_name}, {self.chain_score} score, {self.mapq:0.02} mapq, {len(self.read_seq)} bases, {len(self.chain)} anchors, [{read_region.left},{read_region.right}) -> {self.ref_name} [{ref_region.left}, {ref_region.right}), {read_region.width}:{ref_region.width})"


class TinymapIndex: 

    w: int 
    k: int
    seq_indices: Dict[str, 'TinymapSeqIndex']
    logger: logging.Logger 
    fasta_index: FastaIndex

    def __init__(self, w: int, k: int, ref_seqs: Iterable[str]): 
        self.logger = logging.getLogger("TinymapIndex")
        self.w = w 
        self.k = k 
        self.logger.info(f"TinymapIndex({w=}, {k=}, {ref_seqs=})")
        self.seq_indices = {}
        self.fasta_index = load_fasta_index(Path("index.csv"))
        for name in ref_seqs: 
            self.add_seq(name)
    
    def add_seq(self, seq_name: str): 
        if not self.load_seq(seq_name): 
            seq = self.fasta_index.load_seq(seq_name)
            self.index_seq(seq_name, seq)
            self.logger.info(f"Saving {seq_name}")
            self.seq_indices[seq_name].save()

    def load_seq(self, seq_name: str) -> bool: 
        p = Path(f"index_{seq_name}.json.gz")
        if p.exists(): 
            self.logger.info(f"Loading {seq_name} from {p.name}")
            self.seq_indices[seq_name] = TinymapSeqIndex.load(p) 
            return True
        else: 
            self.logger.info(f"{p.name} doesn't exist")
            return False
    
    def save(self): 
        for (seq_name, idx) in self.seq_indices.items(): 
            self.logger.info(f"Saving {seq_name}")
            idx.save()
    
    def index_seq(self, seq_name: str, seq: str): 
        self.logger.info(f"Indexing {seq_name[0:100]}")
        self.seq_indices[seq_name] = ( 
            TinymapSeqIndex(w=self.w, k=self.k, seq_name=seq_name, seq=seq, seeds=create_seeds_df(self.w, self.k, seq))
        )
    
    def align_next(self, reads: Generator[Read, None, None]) -> Tuple[Read, List[Alignment]]: 
        while True: 
            read = next(reads) 
            aligns = self.align(read) 
            if len(aligns) > 0: 
                return (read, aligns) 
    
    def aligns(self, reads: List[Read]) -> Dict[Read, List[Alignment]]: 
        d = {} 
        for read in track(reads, description=f"Aligning {len(reads)} reads"): 
            d[read] = self.align(read)
        return d
    
    def base_align(self, align: Alignment) -> Tuple[str, str]: 
        params = AlignmentParams(2.0, 3.0, 0.75) 
        def align_middle_chunk(chunk: Tuple[Region, Region]) -> Tuple[str, str]:
            read_str = chunk[0].extract(align.read_seq)
            ref_str = chunk[1].extract(self.seq_indices[align.ref_name].seq)
            aligner = Aligner(params, read_str, ref_str)
            aligner.align()
            return (aligner.gapped_left_string, aligner.gapped_right_string) 
        def align_right_chunk(chunk: Tuple[Region, Region]) -> Tuple[str, str]:
            read_str = chunk[0].extract(align.read_seq)
            ref_str = chunk[1].extract(self.seq_indices[align.ref_name].seq)
            aligner = Aligner(params, read_str, ref_str)
            aligner.align(left_only=True)
            gls = aligner.gapped_left_string.rstrip(' ')
            grs = aligner.gapped_right_string[:len(gls)]
            return (gls, grs)
        def align_left_chunk(chunk: Tuple[Region, Region]) -> Tuple[str, str]:
            read_str = chunk[0].extract(align.read_seq)[::-1]
            ref_str = chunk[1].extract(self.seq_indices[align.ref_name].seq)[::-1]
            aligner = Aligner(params, read_str, ref_str)
            aligner.align(left_only=True)
            gls = aligner.gapped_left_string.rstrip(' ')
            grs = aligner.gapped_right_string[:len(gls)]
            return (gls[::-1], grs[::-1])

        chunks = align.alignable_chunks() 
        [left_chunk, *middle_chunks, right_chunk] = chunks
        left_align = align_left_chunk(left_chunk)
        middle_aligns = [align_middle_chunk(c) for c in middle_chunks]
        right_align = align_right_chunk(right_chunk) 
        aligns = [left_align, *middle_aligns, right_align] 
        return (
            ''.join(x[0] for x in aligns), 
            ''.join(x[1] for x in aligns) 
        )
    
    def align(self, read: Read) -> List[Alignment]: 
        self.logger.debug(f"Aligning {read}")
        hits= [
            align 
            for (_, idx) in self.seq_indices.items()
            for align in idx.align(read)
        ]
        hits.sort()
        self.logger.debug(f"# Hits={len(hits)}")
        return hits
    

@dataclass
class TinymapSeqIndex: 

    @staticmethod 
    def load(p: Path) -> 'TinymapSeqIndex': 
        logger = logging.getLogger("TinymapSeqIndex")
        with p.open('rt') as inf: 
            d = json.loads(inf.read().encode('UTF-8'))
            logger.info(f"Loaded index manifest {d}")
            w = d.get('w') 
            k = d.get('k') 
            seq_name = d.get('seq_name') 
            seq_path = d.get('seq') 
            seed_path = d.get('seeds')
        seed_df = pd.read_parquet(Path(seed_path)) 
        logger.info(f"Seed count: {len(seed_df)}")
        seqs = read_fasta(Path(seq_path))
        seq = seqs[seq_name]
        logger.info(f"Sequence length: {len(seq)}")
        return TinymapSeqIndex(
            w, k, seq_name, seq, seed_df
        )

    w: int 
    k: int 
    seq_name: str 
    seq: str = field(repr=False, compare=False, hash=False)
    seeds: pd.DataFrame = field(repr=False, compare=False, hash=False)
    search_width: int = field(default=1000)

    def save(self): 
        p = Path(f'index_{self.seq_name}.json')
        seed_path = p.with_suffix('.parquet') 
        fasta_path = p.with_suffix('.fasta.gz')
        with gzip.open(p, 'wt') as gout: 
            gout.write(json.dumps({
                "w": self.w, 
                "k": self.k, 
                "seq_name": self.seq_name, 
                "seq": fasta_path.name, 
                "seeds": seed_path.name
            }))
        with gzip.open(fasta_path, 'wt') as outf: 
            outf.write(f">{self.seq_name}\n")
            for offset in range(0, len(self.seq), 100): 
                outf.write(f"{self.seq[offset:offset+100]}\n")
        self.seeds.to_parquet(seed_path)
    
    def align(self, read: Read) -> List[Alignment]: 
        chains = self.align_to_seeds(read.seq) 
        return [
            Alignment(read_name=read.name, read_seq=read.seq, ref_name=self.seq_name, chain=chain) 
            for chain in chains
            if chain.score > 0.0
        ]

    def align_to_seeds(self, read: str) -> List[Chain]: 
        logger = logging.getLogger(f"TinymapIndex_{self.seq_name}")
        logger.log(logging.DEBUG, f"align_to_seeds({read=})")
        read_minimizers: List[Tuple[str, int]] = ( 
            list(minimizers(self.w, self.k, read))
        )
        min_strs: List[str] = [x[0] for x in read_minimizers]
        hit_df = self.seeds[self.seeds['kmer'].isin(min_strs)]

        hit_kmers = hit_df['kmer'].to_list()
        hit_offsets = hit_df['offset'].to_list()
        kmer_lists = { k: [x[1] for x in v] for (k, v) in groupby(zip(hit_kmers, hit_offsets), key=lambda x: x[0]) }
        logger.log(logging.DEBUG, f"{kmer_lists=}")
        hits: List[Anchor] = [
            Anchor(read_offset=read_offset, reference_offset=ref_offset, k=self.k) 
            for (kmer, read_offset) in read_minimizers
            for ref_offset in kmer_lists.get(kmer, [])
        ]
        logger.log(logging.DEBUG, f"# hits: {len(hits)}")

        hits.sort()
        params = ChainParams(w = self.w, max_distance=self.search_width)
        self.chainer = Chainer(params, hits)
        self.chainer.chain_align()
        return self.chainer.chains

def index_seq(w: int, k: int, name: str, seq: str) -> TinymapIndex:
    tt = TinymapIndex(w, k, name, seq, create_seeds_df(w, k, seq)) 
    return tt

