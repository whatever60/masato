set data_collection_dir /mnt/c/Users/quym/Dropbox/Hyperspectral_imaging/data_collection
set data_name 202304_darpa_arcadia_soil
set amplicon_dir /mnt/c/aws_data/20240307_arcadia_soil_interaction_reanalysis

./isolate_utils.py get_metadata \
    -i $data_collection_dir/camii_pick/$data_name \
    -pm $data_collection_dir/plate_metadata/$data_name.csv \
    -o $amplicon_dir/interaction/isolate_metadata.tsv

./isolate_utils.py combine_count \
    -m $amplicon_dir/interaction/isolate_metadata.tsv \
    -cb $amplicon_dir/output_unoise3/unoise3_zotu.tsv \
    -tb $amplicon_dir/output_unoise3/unoise3_zotu_rrndb_processed.tsv \
    -o $amplicon_dir/interaction
    # -cf $amplicon_dir/output_unoise3_its/unoise3_zotu.tsv \
    # -tf $amplicon_dir/output_unoise3_its/unoise3_zotu_rrndb_processed.tsv \

./isolate_utils.py combine_count \
    -m $amplicon_dir/interaction/isolate_metadata.tsv \
    -cb $amplicon_dir/output_unoise3/unoise3_zotu.tsv \
    -tb $amplicon_dir/output_unoise3/unoise3_zotu_rrndb_processed.tsv \
    -o $amplicon_dir/interaction \
    -pl otu domain \
    -cp 0.5
    # -cf $amplicon_dir/output_unoise3_its/unoise3_zotu.tsv \
    # -tf $amplicon_dir/output_unoise3_its/unoise3_zotu_rrndb_processed.tsv \

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
