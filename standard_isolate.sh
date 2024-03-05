set data_collection_dir /mnt/c/Users/quym/Dropbox/Hyperspectral_imaging/data_collection
set data_name 202401_darpa_arcadia_arl_boneyard_b1b2_1
set amplicon_dir /mnt/c/aws_data/20240224_logan_tyndall_boneyard_interaction

./isolate_utils.py get_metadata \
    -i $data_collection_dir/camii_pick/$data_name \
    -pm $data_collection_dir/plate_metadata/$data_name.csv \
    -o $amplicon_dir/interaction/isolate_metadata.tsv

./isolate_utils.py combine_count \
    -m $amplicon_dir/interaction/isolate_metadata.tsv \
    -cb $amplicon_dir/16s/output_unoise3/unoise3_zotu.tsv \
    -cf $amplicon_dir/its/output_unoise3/unoise3_zotu.tsv \
    -tb $amplicon_dir/16s/output_unoise3/unoise3_zotu_rrndb_processed.tsv \
    -tf $amplicon_dir/its/output_unoise3/unoise3_zotu_rrndb_processed.tsv \
    -o $amplicon_dir/interaction


./isolate_utils.py combine_count \
    -m $amplicon_dir/interaction/isolate_metadata.tsv \
    -cb $amplicon_dir/16s/output_unoise3/unoise3_zotu.tsv \
    -cf $amplicon_dir/its/output_unoise3/unoise3_zotu.tsv \
    -tb $amplicon_dir/16s/output_unoise3/unoise3_zotu_rrndb_processed.tsv \
    -tf $amplicon_dir/its/output_unoise3/unoise3_zotu_rrndb_processed.tsv \
    -o $amplicon_dir/interaction \
    -pl otu domain