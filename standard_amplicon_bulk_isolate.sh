
set amplicon_type its
set working_dir /mnt/c/aws_data/20240315_arl_boneyard_bulk_isolate_i_redo
set fastq_dir fastq
set output_unoise3_dir output_unoise3_$amplicon_type
set metadata_file metadata/bulk.tsv
set db_16s_local /mnt/c/aws_data/data/blast/db/16s_ribosomal_rna/16S_ribosomal_RNA
set db_16s_online 16S_ribosomal_RNA
set db_its_local /mnt/c/aws_data/data/blast/db/its_refseq_fungi/ITS_RefSeq_Fungi
set db_its_online ITS_RefSeq_Fungi
set db_its_unite_local /mnt/c/aws_data/data/unite/general/sh_general_release_dynamic_25.07.2023_dev_simple
set db_nt_online nt
if test $amplicon_type = 16s
    set rdp_db 16srrna
    set prefix 16S-O
    set unknown_name 16S-O_UNKNOWN
    set level genus
    set fig_dir figs_16s
    set domain bacteria
else
    set rdp_db fungalits_unite
    # set rdp_db fungalits_unite_v9
    set prefix ITS-O
    set level genus species
    set unknown_name ITS-O_UNKNOWN
    set fig_dir figs_its
    set domain fungi
end


mkdir -p $working_dir/$output_unoise3_dir
./trim.py \
    --mode isolate_150 \
    -p $amplicon_type \
    -i $working_dir/fastq_isolate \
    -o $working_dir/$output_unoise3_dir/merged_isolate.fq.gz \
    -fb data/isolate/barcodes_fwd.fasta \
    -rb data/isolate/barcodes_rev.fasta \
    -pt data/isolate/patterns.txt \
    -k 110 \
    -l 110

./trim.py \
    --mode r1 \
    -i $working_dir/fastq_bulk_$amplicon_type \
    -o $working_dir/$output_unoise3_dir/merged_bulk.fq.gz \
    -p $amplicon_type \
    -k 110 \
    -l 110

zcat $working_dir/$output_unoise3_dir/merged_isolate.fq.gz $working_dir/$output_unoise3_dir/merged_bulk.fq.gz | gzip > $working_dir/$output_unoise3_dir/merged.fq.gz

./usearch_workflow.py unoise3 \
    -i $working_dir/$output_unoise3_dir/merged.fq.gz \
    -o $working_dir/$output_unoise3_dir/unoise3_zotu.fa \
    -l $prefix

rdp_classifier \
    $working_dir/$output_unoise3_dir/unoise3_zotu.fa \
    -d $rdp_db \
    -o $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb.tsv \
    -h $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb.hier.tsv

./process_rrndb.py \
    -i $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb.tsv \
    -o $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_processed.tsv

# with 16S database
rdp_classifier \
    $working_dir/$output_unoise3_dir/unoise3_zotu.fa \
    -d 16srrna \
    -o $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_16s.tsv \
    -h $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_16s.hier.tsv
./process_rrndb.py \
    -i $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_16s.tsv \
    -o $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_16s_processed.tsv

# with ITS database
rdp_classifier \
    $working_dir/$output_unoise3_dir/unoise3_zotu.fa \
    -d fungalits_unite_v9 \
    -o $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_its.tsv \
    -h $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_its.hier.tsv
./process_rrndb.py \
    -i $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_its.tsv \
    -o $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_its_processed.tsv

./usearch_workflow.py search_global \
    -i $working_dir/$output_unoise3_dir/merged.fq.gz \
    -d $working_dir/$output_unoise3_dir/unoise3_zotu.fa \
    --unknown_name $unknown_name

mkdir -p $working_dir/$fig_dir
./plot.py abundance_group heatmap_log10 \
    -i $working_dir/$output_unoise3_dir/unoise3_zotu.tsv \
    -m $working_dir/$metadata_file \
    -t $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_processed.tsv \
    -f $working_dir/$fig_dir/abundance_isolate \
    -l $level \
    -fo taxonomy_tree \
    -a 0.001 \
    -ii $working_dir/interaction/count_otu_$domain.tsv \
    -im $working_dir/interaction/isolate_metadata.tsv \
    -ir sample_type