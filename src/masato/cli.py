import os
import sys
from typing import Annotated

import cyclopts
from cyclopts import App, Parameter

from .utils import cat_fastq as _cat_fastq

app = App(
    help="Concatenate paired-end FASTQ files.",
    help_on_error=True,
    version_flags=("--version", "-v"),
)


@app.command
def cat_fastq(
    input_pos: Annotated[
        list[str] | None,
        Parameter(help="Positional input: directory, glob, or list of FASTQ files."),
    ] = None,
    input_kw: Annotated[
        list[str] | None,
        Parameter(
            name=["--input", "-i"],
            help="Keyword-style input: directory, glob, or list of FASTQ files.",
            consume_multiple=True,
        ),
    ] = None,
    output1: Annotated[
        str | None,
        Parameter(
            name=["-o1", "--output1"], help="Output file or stream for R1 reads."
        ),
    ] = None,
    output2: Annotated[
        str | None,
        Parameter(
            name=["-o2", "--output2"], help="Output file or stream for R2 reads."
        ),
    ] = None,
    sample_list: Annotated[
        str | None,
        Parameter(name="--sample-list", help="Optional sample metadata file."),
    ] = None,
    remove_undet: Annotated[
        bool,
        Parameter(name="--remove-undet", help="Skip samples named 'Undetermined'."),
    ] = True,
    have_sample_name: Annotated[
        bool,
        Parameter(
            name="--have-sample-name", help="Include sample name in renamed read ID."
        ),
    ] = False,
):
    """
    Concatenate paired-end FASTQ files from keyword or positional input,
    with optional renaming and filtering.
    """
    if not input_kw and not input_pos:
        print("Error: No input provided.", file=sys.stderr)
        raise SystemExit(1)
    if input_kw and input_pos:
        print(
            "Error: Provide input either via --input/-i or positional, not both.",
            file=sys.stderr,
        )
        raise SystemExit(1)

    input_ = input_kw or input_pos
    if not input_:
        print("Error: No input provided.", file=sys.stderr)
        raise SystemExit(1)

    if len(input_) == 1 and os.path.isdir(input_[0]):
        resolved_input = input_[0]
    else:
        resolved_input = input_

    _cat_fastq(
        directory=resolved_input,
        output_fp_r1=output1,
        output_fp_r2=output2,
        metadata=sample_list,
        _remove_undet=remove_undet,
        _have_sample_name=have_sample_name,
    )


def cat_fastq_cli():
    cyclopts.run(cat_fastq)
