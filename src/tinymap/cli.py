import click
from .swalign import Aligner, AlignmentParams


@click.group()
def main(): ...


@main.command("align")
@click.argument("seq1")
@click.argument("seq2")
def align_sequences(seq1: str, seq2: str):
    params = AlignmentParams(mismatch=1.0, gap_open=1.25, gap_extend=0.5)
    aligner = Aligner(params=params, left=seq1, right=seq2)
    aligner.align()
    click.echo(aligner.alignment)
