
set amplicon_type 16s
set working_dir /mnt/c/aws_data/20240306_arl_boneyard_bulk_isolate_i
set fastq_dir fastq
set output_unoise3_dir output_unoise3_$amplicon_type
set metadata_file metadata/bulk.tsv
if test $amplicon_type = 16s
    set rdp_db 16srrna
    set prefix 16S-ZOTU
    set unknown_name 16S-ZOTU_UNKNOWN
    set level genus 
else
    set rdp_db fungalits_unite
    set rdp_db fungalits_unite_v9
    set prefix ITS-ZOTU
    set level genus species
    set unknown_name ITS-ZOTU_UNKNOWN
end

set fig_dir figs

# ./trim.py --mode simple -i $working_dir/$fastq_dir -o $working_dir/$output_unoise3_dir/merged.fq.gz
./usearch_workflow.py unoise3 \
    -i $working_dir/$output_unoise3_dir/merged.fq.gz \
    -o $working_dir/$output_unoise3_dir/unoise3_zotu.fa \
    -l $prefix
./usearch_workflow.py search_global \
    -i $working_dir/$output_unoise3_dir/merged.fq.gz \
    -d $working_dir/$output_unoise3_dir/unoise3_zotu.fa \
    --unknown_name $unknown_name

rdp_classifier \
    $working_dir/$output_unoise3_dir/unoise3_zotu.fa \
    -d $rdp_db \
    -o $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb.tsv \
    -h $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb.hier.tsv

./process_rrndb.py -i $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb.tsv -o  $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_processed.tsv

./plot.py abundance_group all \
    -i $working_dir/$output_unoise3_dir/unoise3_zotu.tsv \
    -m $working_dir/$metadata_file \
    -t $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_processed.tsv \
    -f $working_dir/$fig_dir/abundance \
    -l genus family order \
    -a 0.01

./plot.py abundance_group all \
    -i $working_dir/$output_unoise3_dir/unoise3_zotu.tsv \
    -m $working_dir/$metadata_file \
    -t $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_processed.tsv \
    -f $working_dir/$fig_dir/abundance_top \
    -l $level \
    -a 1

./plot.py abundance_group all \
    -i $working_dir/$output_unoise3_dir/unoise3_zotu.tsv \
    -m $working_dir/$metadata_file \
    -t $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_processed.tsv \
    -f $working_dir/$fig_dir/abundance_clust \
    -l genus family order \
    -fc \
    -a 0.01

./plot.py abundance_group all \
    -i $working_dir/$output_unoise3_dir/unoise3_zotu.tsv \
    -m $working_dir/$metadata_file \
    -t $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_processed.tsv \
    -f $working_dir/$fig_dir/abundance_top_clust \
    -l $level \
    -fc \
    -a 1

./plot.py stats_sample_count \
    -i $working_dir/$output_unoise3_dir/unoise3_zotu.tsv \
    -m $working_dir/$metadata_file \
    -t $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_processed.tsv \
    -f $working_dir/$fig_dir/alpha_diversity_rarefying \
    -l otu genus family order \
    -rv 500 \
    -rn 50

./plot.py stats_sample_count \
    -i $working_dir/$output_unoise3_dir/unoise3_zotu.tsv \
    -m $working_dir/$metadata_file \
    -t $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_processed.tsv \
    -f $working_dir/$fig_dir/alpha_diversity \
    -l otu genus family order

./plot_dm.py \
    -i $working_dir/$output_unoise3_dir/unoise3_zotu.tsv \
    -m $working_dir/$metadata_file \
    -t $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_processed.tsv \
    -f $working_dir/$fig_dir/beta_diversity \
    -l otu genus \
    -d euclid --log10 -tr relative \
    --hue location_short

./plot_dm.py \
    -i $working_dir/$output_unoise3_dir/unoise3_zotu.tsv \
    -m $working_dir/$metadata_file \
    -t $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_processed.tsv \
    -f $working_dir/$fig_dir/beta_diversity \
    -l otu genus \
    -d braycurtis \
    --hue location_short

./plot_dm.py \
    -i $working_dir/$output_unoise3_dir/unoise3_zotu.tsv \
    -m $working_dir/$metadata_file \
    -t $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_processed.tsv \
    -f $working_dir/$fig_dir/beta_diversity_annot \
    -l otu genus \
    -d euclid --log10 -tr relative \
    --hue location_short \
    -a

./plot_dm.py \
    -i $working_dir/$output_unoise3_dir/unoise3_zotu.tsv \
    -m $working_dir/$metadata_file \
    -t $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_processed.tsv \
    -f $working_dir/$fig_dir/beta_diversity_annot \
    -l otu genus \
    -d braycurtis \
    --hue location_short \
    -a
