
set amplicon_type 16s
# set working_dir /mnt/c/aws_data/20240318_arl_boneyard_bulk_isolate_ii
set working_dir /mnt/c/aws_data/20240315_arl_boneyard_bulk_isolate_i_redo
set fastq_dir fastq
set output_unoise3_dir output_unoise3_$amplicon_type
set metadata_file metadata/bulk_$amplicon_type.tsv
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

./get_tree.py -i $working_dir/$output_unoise3_dir/unoise3_zotu.fa

set data_collection_dir /mnt/c/Users/quym/Dropbox/Hyperspectral_imaging/data_collection
set data_name 202401_darpa_arcadia_arl_boneyard_b1b2_1
# set data_name 202401_darpa_arcadia_arl_boneyard_b3b4_1

mkdir -p $working_dir/interaction
./isolate_utils.py get_metadata \
    -i $data_collection_dir/camii_pick/$data_name \
    -pm $data_collection_dir/plate_metadata/$data_name.csv \
    -o $working_dir/interaction/isolate_metadata.tsv

./isolate_utils.py combine_count \
    -m $working_dir/interaction/isolate_metadata.tsv \
    -cb $working_dir/output_unoise3_16s/unoise3_zotu.tsv \
    -tb $working_dir/output_unoise3_16s/unoise3_zotu_rrndb_processed.tsv \
    -cf $working_dir/output_unoise3_its/unoise3_zotu.tsv \
    -tf $working_dir/output_unoise3_its/unoise3_zotu_rrndb_processed.tsv \
    -o $working_dir/interaction

./isolate_utils.py combine_count \
    -m $working_dir/interaction/isolate_metadata.tsv \
    -cb $working_dir/output_unoise3_16s/unoise3_zotu.tsv \
    -tb $working_dir/output_unoise3_16s/unoise3_zotu_rrndb_processed.tsv \
    -cf $working_dir/output_unoise3_its/unoise3_zotu.tsv \
    -tf $working_dir/output_unoise3_its/unoise3_zotu_rrndb_processed.tsv \
    -o $working_dir/interaction \
    -pl otu domain \
    -cp 0.5

./plot.py abundance_group heatmap_log10 \
    -i $working_dir/$output_unoise3_dir/unoise3_zotu.tsv \
    -m $working_dir/$metadata_file \
    -t $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_processed.tsv \
    -f $working_dir/figs/abundance_isolate_$amplicon_type \
    -a 0 \
    -l $level \
    -fo taxonomy_tree \
    -ii $working_dir/interaction/count_otu_$domain.tsv \
    -im $working_dir/interaction/isolate_metadata.tsv \
    -ir sample_type \
    -isp spike_in_$amplicon_type

./plot.py abundance_group heatmap_log10 \
    -i $working_dir/$output_unoise3_dir/unoise3_zotu.tsv \
    -m $working_dir/$metadata_file \
    -t $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_processed.tsv \
    -f $working_dir/figs/abundance_isolate_top_$amplicon_type \
    -a 5 \
    -l $level \
    -fo taxonomy_tree \
    -ii $working_dir/interaction/count_otu_$domain.tsv \
    -im $working_dir/interaction/isolate_metadata.tsv \
    -ir sample_type \
    -isp spike_in_$amplicon_type


./plot.py abundance_group heatmap_raw \
    -i $working_dir/interaction/count_otu_$domain.tsv \
    -m $working_dir/interaction/isolate_metadata.tsv \
    -t $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_processed.tsv \
    -f $working_dir/figs/isolate_count_$amplicon_type \
    -l $level \
    -a 0.01 \
    -r plate_name \
    -s sample_group \
    -sp spike_in_$amplicon_type \
    -fo taxonomy_tree

# ./plot.py abundance_group heatmap_raw \
#         -i $working_dir/interaction/count_otu_fungi.tsv \
#         -m $working_dir/interaction/isolate_metadata.tsv \
#         -t $working_dir/output_unoise3_its/unoise3_zotu_rrndb_processed.tsv \
#         -f $working_dir/figs/isolate_count_its \
#         -l species \
#         -r plate_name \
#         -s sample_group
