import typer

from .basic import combine_trim_clip_pe as _combine_trim_clip_pe
from .camii_isolates import isolate_150_preprocess as _isolate_150_preprocess
from ..trim import resolve_input_path


app = typer.Typer()


@app.command()
def combine_trim_clip_pe(
    input_fastq_dir: str = typer.Option(
        ..., "-i", "--input-dir", help="Input interleaved FASTQ path"
    ),
    output_fastq: str = typer.Option(
        ..., "-o1", "--output1", help="Output FASTQ for read 1"
    ),
    output_fastq_r2: str = typer.Option(
        ..., "-o2", "--output2", help="Output FASTQ for read 2"
    ),
    primer_set: str = typer.Option(..., help="Primer set name"),
    first_k: int | None = typer.Option(
        None, "-k1", "--first-k1", help="First k bases to keep from read 1"
    ),
    first_k_r2: int | None = typer.Option(
        None, "-k2", "--first-k2", help="First k bases to keep from read 2"
    ),
    min_length: int = typer.Option(
        0, "-l1", "--min-length1", help="Minimum length for read 1"
    ),
    min_length_r2: int = typer.Option(
        0, "-l2", "--min-length2", help="Minimum length for read 2"
    ),
    quality_trimming: bool = typer.Option(
        False, help="Enable quality trimming with fastp"
    ),
    cores: int = typer.Option(16, help="Number of cores to use"),
) -> None:
    return _combine_trim_clip_pe(
        input_fastq_dir=input_fastq_dir,
        output_fastq=output_fastq,
        output_fastq_r2=output_fastq_r2,
        primer_set=primer_set,
        first_k=first_k,
        first_k_r2=first_k_r2,
        min_length=min_length,
        min_length_r2=min_length_r2,
        quality_trimming=quality_trimming,
        cores=cores,
    )


@app.command()
def isolate_150_preprocess(
    fastq_dir: str = typer.Option(..., "-i", "--input-dir", help="Input FASTQ directory"),
    barcode_fwd_fasta: str = typer.Option(..., "-fb", "--barcode-fwd", help="Forward barcode FASTA"),
    barcode_rev_fasta: str = typer.Option(..., "-rb", "--barcode-rev", help="Reverse barcode FASTA"),
    rename_pattern: str = typer.Option(..., "-pt", "--rename-pattern", help="Renaming pattern for output reads"),
    output_fastq: str = typer.Option(..., "-o1", "--output1", help="Output FASTQ for read 1"),
    output_fastq_r2: str = typer.Option(..., "-o2", "--output2", help="Output FASTQ for read 2"),
    # args for combine_trim_clip_pe
    primer_set: str = typer.Option(..., help="Primer set name"),
    first_k: int | None = typer.Option(
        None, "-k1", "--first-k1", help="First k bases to keep from read 1"
    ),
    first_k_r2: int | None = typer.Option(
        None, "-k2", "--first-k2", help="First k bases to keep from read 2"
    ),
    min_length: int = typer.Option(
        0, "-l1", "--min-length1", help="Minimum length for read 1"
    ),
    min_length_r2: int = typer.Option(
        0, "-l2", "--min-length2", help="Minimum length for read 2"
    ),
    # early_stop: bool = False,
    quality_trimming: bool = typer.Option(
        False, help="Enable quality trimming with fastp"
    ),
    cores: int = typer.Option(16, help="Number of cores to use"),
) -> None:
    barcode_fwd_fasta = resolve_input_path(barcode_fwd_fasta)
    barcode_rev_fasta = resolve_input_path(barcode_rev_fasta)
    rename_pattern = resolve_input_path(rename_pattern)
    return _isolate_150_preprocess(
        fastq_dir=fastq_dir,
        barcode_fwd_fasta=barcode_fwd_fasta,
        barcode_rev_fasta=barcode_rev_fasta,
        rename_pattern=rename_pattern,
        output_fastq=output_fastq,
        output_fastq_r2=output_fastq_r2,
        primer_set=primer_set,
        first_k=first_k,
        first_k_r2=first_k_r2,
        min_length=min_length,
        min_length_r2=min_length_r2,
        quality_trimming=quality_trimming,
        cores=cores,
    )