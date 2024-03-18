cutadapt \
    -e 0.15 \
    -a "GTGTGYCAGCMGCCGCGGTAA;required...ATTAGAWACCCBNGTAGTCCGG;optional" \
    -A "CCGGACTACNVGGGTWTCTAAT;required...TTACCGCGGCKGCTGRCACAC;optional" \
    --minimum-length 110 \
    --pair-filter any \
    -o \
    ./temp/merged_isolate_R1.fq.gz \
    -p \
    ./temp/merged_isolate_R2.fq.gz \
    --untrimmed-output \
    ./temp/untrimmed_1.fq.gz \
    --untrimmed-paired-output \
    ./temp/untrimmed_2.fq.gz \
    --too-short-output \
    ./temp/too_short_1.fq.gz \
    --too-short-paired-output \
    ./temp/too_short_2.fq.gz \
    --cores 4 \
    -l 110 \
    /mnt/c/aws_data/20240315_arl_boneyard_bulk_isolate_i_redo/fastq_isolate/ABYPD4A_S4_L001_R1_001.fastq.gz \
    /mnt/c/aws_data/20240315_arl_boneyard_bulk_isolate_i_redo/fastq_isolate/ABYPD4A_S4_L001_R2_001.fastq.gz


cutadapt \
    -e 0.15 \
    -g ^file:data/isolate/barcodes_fwd.fasta \
    -G ^file:data/isolate/barcodes_rev.fasta \
    -o ./temp/{name1}-{name2}_R1.fq.gz \
    -p ./temp/{name1}-{name2}_R2.fq.gz \
    /mnt/c/aws_data/20240315_arl_boneyard_bulk_isolate_i_redo/fastq_isolate/ABYPD4A_S4_L001_R1_001.fastq.gz \
    /mnt/c/aws_data/20240315_arl_boneyard_bulk_isolate_i_redo/fastq_isolate/ABYPD4A_S4_L001_R2_001.fastq.gz
    # --no-indels \

    cutadapt \
        -a \
        GTGTGYCAGCMGCCGCGGTAA...ATTAGAWACCCBNGTAGTCCGG \
        -A \
        CCGGACTACNVGGGTWTCTAAT...TTACCGCGGCKGCTGRCACAC \
        --minimum-length \
        110 \
        --pair-filter \
        first \
        -o \
        ./temp/merged_isolate_R1.fq.gz \
        -p \
        ./temp/merged_isolate_R2.fq.gz \
        --cores \
        4 \
        /mnt/c/aws_data/20240306_arl_boneyard_bulk_isolate_i/fastq_isolate/ABYPD2A_S2_L001_R1_001.fastq.gz \
        /mnt/c/aws_data/20240306_arl_boneyard_bulk_isolate_i/fastq_isolate/ABYPD2A_S2_L001_R2_001.fastq.gz