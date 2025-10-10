from typing import Annotated

from cyclopts import App, Parameter

from .basic import (
    combine_trim_clip_pe as _combine_trim_clip_pe,
    combine_trim_merge_pe as _combine_trim_merge_pe,
)
from .camii_isolates import isolate_150_preprocess as _isolate_150_preprocess
from .rnatag_seq import rnatagseq_150_preprocess as _rnatagseq_150_preprocess
from ..trim import resolve_input_path


app = App()


@app.command
def combine_trim_clip_pe(
    input_fastq_dir: Annotated[
        str, Parameter(name=["-i", "--input-dir"], help="Input interleaved FASTQ path")
    ],
    output_fastq: Annotated[
        str, Parameter(name=["-o1", "--output1"], help="Output FASTQ for read 1")
    ],
    output_fastq_r2: Annotated[
        str, Parameter(name=["-o2", "--output2"], help="Output FASTQ for read 2")
    ],
    primer_set: Annotated[str | None, Parameter(help="Primer set name")] = None,
    first_k: Annotated[
        int | None,
        Parameter(name=["-k1", "--first-k1"], help="First k bases to keep from read 1"),
    ] = None,
    first_k_r2: Annotated[
        int | None,
        Parameter(name=["-k2", "--first-k2"], help="First k bases to keep from read 2"),
    ] = None,
    min_length: Annotated[
        int, Parameter(name=["-l1", "--min-length1"], help="Minimum length for read 1")
    ] = 0,
    min_length_r2: Annotated[
        int, Parameter(name=["-l2", "--min-length2"], help="Minimum length for read 2")
    ] = 0,
    quality_trimming: Annotated[
        bool, Parameter(help="Enable quality trimming with fastp")
    ] = False,
    cores: Annotated[int, Parameter(help="Number of cores to use")] = 16,
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


@app.command
def isolate_150_preprocess(
    fastq_dir: Annotated[
        str, Parameter(name=["-i", "--input-dir"], help="Input FASTQ directory")
    ],
    barcode_fwd_fasta: Annotated[
        str, Parameter(name=["-fb", "--barcode-fwd"], help="Forward barcode FASTA")
    ],
    barcode_rev_fasta: Annotated[
        str, Parameter(name=["-rb", "--barcode-rev"], help="Reverse barcode FASTA")
    ],
    rename_pattern: Annotated[
        str,
        Parameter(
            name=["-pt", "--rename-pattern"], help="Renaming pattern for output reads"
        ),
    ],
    output_fastq: Annotated[
        str, Parameter(name=["-o1", "--output1"], help="Output FASTQ for read 1")
    ],
    output_fastq_r2: Annotated[
        str, Parameter(name=["-o2", "--output2"], help="Output FASTQ for read 2")
    ],
    # args for combine_trim_clip_pe
    primer_set: Annotated[str, Parameter(help="Primer set name")],
    first_k: Annotated[
        int | None,
        Parameter(name=["-k1", "--first-k1"], help="First k bases to keep from read 1"),
    ] = None,
    first_k_r2: Annotated[
        int | None,
        Parameter(name=["-k2", "--first-k2"], help="First k bases to keep from read 2"),
    ] = None,
    min_length: Annotated[
        int, Parameter(name=["-l1", "--min-length1"], help="Minimum length for read 1")
    ] = 0,
    min_length_r2: Annotated[
        int, Parameter(name=["-l2", "--min-length2"], help="Minimum length for read 2")
    ] = 0,
    # early_stop: bool = False,
    quality_trimming: Annotated[
        bool, Parameter(help="Enable quality trimming with fastp")
    ] = False,
    cores: Annotated[int, Parameter(help="Number of cores to use")] = 16,
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


@app.command
def combine_trim_merge_pe(
    input_dir: Annotated[
        str, Parameter(name=["-i", "--input-dir"], help="Input FASTQ directory")
    ],
    output_fastq: Annotated[
        str, Parameter(name=["-o", "--output"], help="Output path for merged FASTQ")
    ],
    primer_set: Annotated[str | None, Parameter(help="Primer set name")] = None,
    overlap_len_require: Annotated[
        int | None,
        Parameter(
            name="--overlap-len-require",
            help="Required overlap length for merging reads",
        ),
    ] = None,
    overlap_diff_limit: Annotated[
        int | None,
        Parameter(
            name="--overlap-diff-limit",
            help="Maximum mismatches allowed in the overlap",
        ),
    ] = None,
    overlap_diff_percent_limit: Annotated[
        float | None,
        Parameter(
            name="--overlap-diff-percent-limit",
            help="Maximum mismatch percentage allowed in the overlap",
        ),
    ] = None,
    min_length: Annotated[
        int | None,
        Parameter(name="--min-length", help="Minimum required length after merging"),
    ] = None,
    cores: Annotated[
        int,
        Parameter(name="--cores", help="Number of cores (threads) for fastp to use"),
    ] = 8,
) -> None:
    """
    Combines, quality trims, and merges paired-end reads from a directory using fastp.
    """
    return _combine_trim_merge_pe(
        input_dir=input_dir,
        output_fastq=output_fastq,
        primer_set=primer_set,
        overlap_len_require=overlap_len_require,
        overlap_diff_limit=overlap_diff_limit,
        overlap_diff_percent_limit=overlap_diff_percent_limit,
        min_length=min_length,
        cores=cores,
    )


@app.command
def rnatagseq_150_preprocess(
    input_: Annotated[
        list[str],
        Parameter(
            name=["-i", "--input"],
            help="Input directory, glob, or list of FASTQ files.",
            consume_multiple=True,
        ),
    ],
    barcode_fasta: Annotated[
        str,
        Parameter(
            name=["-b", "--barcode-fasta"], help="Barcode FASTA file for demultiplexing"
        ),
    ],
    output_dir: Annotated[
        str,
        Parameter(
            name=["-o", "--output-dir"], help="Directory to store processed outputs"
        ),
    ],
    cores: Annotated[int, Parameter(help="Number of CPU threads to use")] = 16,
) -> None:
    """
    Preprocess RNAtag-Seq 150bp reads by trimming, demultiplexing, and writing outputs.

    Args:
        input_: Input FASTQ files or directory/glob pattern.
        barcode_fasta: Path to barcode FASTA file.
        output_dir: Directory where processed FASTQs will be written.
        cores: Number of threads to use.
    """
    barcode_fasta = resolve_input_path(barcode_fasta)
    if len(input_) == 1:
        input_resolved = input_[0]
    else:
        input_resolved = input_
    _rnatagseq_150_preprocess(
        input_=input_resolved,
        barcode_fasta=barcode_fasta,
        output_dir=output_dir,
        cores=cores,
    )
