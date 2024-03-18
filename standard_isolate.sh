set data_collection_dir /mnt/c/Users/quym/Dropbox/Hyperspectral_imaging/data_collection
set data_name 202401_darpa_arcadia_arl_boneyard_b1b2_1
set amplicon_dir /mnt/c/aws_data/20240315_arl_boneyard_bulk_isolate_i_redo

mkdir -p $amplicon_dir/interaction
./isolate_utils.py get_metadata \
    -i $data_collection_dir/camii_pick/$data_name \
    -pm $data_collection_dir/plate_metadata/$data_name.csv \
    -o $amplicon_dir/interaction/isolate_metadata.tsv

./isolate_utils.py combine_count \
    -m $amplicon_dir/interaction/isolate_metadata.tsv \
    -cb $amplicon_dir/output_unoise3_16s/unoise3_zotu.tsv \
    -tb $amplicon_dir/output_unoise3_16s/unoise3_zotu_rrndb_processed.tsv \
    -cf $amplicon_dir/output_unoise3_its/unoise3_zotu.tsv \
    -tf $amplicon_dir/output_unoise3_its/unoise3_zotu_rrndb_processed.tsv \
    -o $amplicon_dir/interaction

./isolate_utils.py combine_count \
    -m $amplicon_dir/interaction/isolate_metadata.tsv \
    -cb $amplicon_dir/output_unoise3_16s/unoise3_zotu.tsv \
    -tb $amplicon_dir/output_unoise3_16s/unoise3_zotu_rrndb_processed.tsv \
    -cf $amplicon_dir/output_unoise3_its/unoise3_zotu.tsv \
    -tf $amplicon_dir/output_unoise3_its/unoise3_zotu_rrndb_processed.tsv \
    -o $amplicon_dir/interaction \
    -pl otu domain \
    -cp 0.5

mkdir -p $working_dir/$fig_dir
./plot.py abundance_group heatmap_log10 \
    -i $working_dir/$output_unoise3_dir/unoise3_zotu.tsv \
    -m $working_dir/$metadata_file \
    -t $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_processed.tsv \
    -f $working_dir/$fig_dir/abundance_isolate \
    -l genus \
    -fo taxonomy_tree \
    -a 0.001 \
    -ii $working_dir/interaction/count_otu_bacteria.tsv \
    -im $working_dir/interaction/isolate_metadata.tsv \
    -ir sample_type

./plot.py abundance_group heatmap_raw \
    -i $amplicon_dir/interaction/count_otu_bacteria.tsv \
    -m $amplicon_dir/interaction/isolate_metadata.tsv \
    -t $amplicon_dir/output_unoise3_16s/unoise3_zotu_rrndb_processed.tsv \
    -f $amplicon_dir/figs/isolate_count_16s \
    -l genus \
    -r plate_name \
    -s sample_group

./plot.py abundance_group heatmap_raw \
    -i $amplicon_dir/interaction/count_otu_fungi.tsv \
    -m $amplicon_dir/interaction/isolate_metadata.tsv \
    -t $amplicon_dir/output_unoise3_its/unoise3_zotu_rrndb_processed.tsv \
    -f $amplicon_dir/figs/isolate_count_its \
    -l species \
    -r plate_name \
    -s sample_group
