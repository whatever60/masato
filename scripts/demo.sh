#!/bin/bash

# Set the data directory from the environment variable
data_dir=/home/whatever60/dev/easy_amplicon_restructure/tests/data
data_dir=~/data/easy_amplicon_camii_test_data
data_dir=/mnt/c/aws_data/data/camii_demo/demo_20241105_output

# Define barcode files
barcodes_fwd="isolate/barcodes_fwd.fasta"
barcodes_rev="isolate/barcodes_rev.fasta"
patterns="isolate/patterns.txt"

# Create output directory
mkdir -p "$data_dir/output_unoise3_16s"

# Run trim.py for isolate_150 mode with barcode arguments
trim.py --mode isolate_150 \
  --input_dir "$data_dir/fastq_isolate" \
  --output "$data_dir/output_unoise3_16s/merged_isolate.fq.gz" \
  --primer_set "16s" \
  --first_k "110" \
  --min_length "110" \
  --barcode_fwd "$barcodes_fwd" \
  --barcode_rev "$barcodes_rev" \
  --pattern "$patterns"

# Run trim.py for r1 mode without barcode arguments
trim.py --mode r1 \
  --input_dir "$data_dir/fastq_bulk_16s" \
  --output "$data_dir/output_unoise3_16s/merged_bulk.fq.gz" \
  --primer_set "16s" \
  --first_k "110" \
  --min_length "110"

# Run trim.py for simple mode without barcode arguments
trim.py --mode simple \
  --input_dir "$data_dir/fastq_bulk_16s" \
  --output "$data_dir/output_unoise3_16s/merged.fq.gz"

# Command 1: Run unoise3
usearch_workflow.py unoise3 \
  --input_fastq "$data_dir/output_unoise3_16s/merged.fq.gz" \
  --output_fasta "$data_dir/output_unoise3_16s/unoise3_zotu.fa" \
  --relabel_prefix "16S-U"

# Before the next command let's subset merged.fq.gz to 100,000 reads for faster processing. In real life, you would not do this.
seqtk sample -s 100 "$data_dir/output_unoise3_16s/merged.fq.gz" 100000 | pigz > "$data_dir/output_unoise3_16s/merged_1m.fq.gz"

# Command 2: Run search_global
usearch_workflow.py search_global \
  --input_fastq "$data_dir/output_unoise3_16s/merged_1m.fq.gz" \
  --zotu_fasta "$data_dir/output_unoise3_16s/unoise3_zotu.fa" \
  --output_tsv "$data_dir/output_unoise3_16s/unoise3_zotu.biom" \
  --unknown_name "16S-U_UNKNOWN"

# Command 3: Run workflow_per_sample
usearch_workflow.py workflow_per_sample \
  --input_fastq "$data_dir/output_unoise3_16s/merged_isolate.fq.gz" \
  --output_path "$data_dir/output_unoise3_16s/unoise3_zotu_isolate.json" \
  --search

# Command 4: Run aggregate_samples
usearch_workflow.py aggregate_samples \
  --input_json "$data_dir/output_unoise3_16s/unoise3_zotu_isolate.json" \
  --output_fasta "$data_dir/output_unoise3_16s/unoise3_zotu_isolate.fa" \
  --output_count "$data_dir/output_unoise3_16s/unoise3_zotu_isolate.biom" \
  --prefix "16S-U"

# Run RDP classifier
run_rdp_classifier.py --input "$data_dir/output_unoise3_16s/unoise3_zotu.fa" --database "16srrna"

# Generate phylogenetic tree
get_tree.py --input_fasta "$data_dir/output_unoise3_16s/unoise3_zotu.fa"

fig_dir_order_pairs=("tax taxonomy_tree" "abc alphabetical" "clust hierarchical")
subcommands=("heatmap" "heatmap_log10" "stacked_bar")
threshold_suffix_pairs=("1e-2 1e-2" "3 top3")

# Loop through each fig_dir and order_option pair
for pair in "${fig_dir_order_pairs[@]}"; do
    # Extract the figure directory and ordering option from the pair
    fig_dir_name=$(echo $pair | cut -d ' ' -f1)
    order_option=$(echo $pair | cut -d ' ' -f2)

    # Loop through each threshold_suffix pair
    for threshold_pair in "${threshold_suffix_pairs[@]}"; do
        # Extract the abundance threshold and directory suffix from the pair
        abundance_threshold=$(echo $threshold_pair | cut -d ' ' -f1)
        suffix=$(echo $threshold_pair | cut -d ' ' -f2)

        # Adjust the figure directory name based on the suffix
        adjusted_fig_dir="${fig_dir_name}_${suffix}"

        # Loop through each subcommand
        for subcommand in "${subcommands[@]}"; do
            # Construct and execute the command
            echo plot.py abundance_group "$subcommand" \
                --otu_count_tsv "$data_dir/output_unoise3_16s/unoise3_zotu.biom" \
                --metadata "$data_dir/metadata/bulk_16s.tsv" \
                --otu_taxonomy_tsv "$data_dir/output_unoise3_16s/unoise3_zotu_rrndb_processed.tsv" \
                --fig_dir "$data_dir/figs/abundance_${adjusted_fig_dir}_16s" \
                --tax_levels "genus" \
                --feature_ordering "$order_option" \
                --rel_ab_thresholds "$abundance_threshold"
            plot.py abundance_group "$subcommand" \
                --otu_count_tsv "$data_dir/output_unoise3_16s/unoise3_zotu.biom" \
                --metadata "$data_dir/metadata/bulk_16s.tsv" \
                --otu_taxonomy_tsv "$data_dir/output_unoise3_16s/unoise3_zotu_rrndb_processed.tsv" \
                --fig_dir "$data_dir/figs/abundance_${adjusted_fig_dir}_16s" \
                --tax_levels "genus" \
                --feature_ordering "$order_option" \
                --rel_ab_thresholds "$abundance_threshold"
        done
    done
done

# Plot alpha diversity
plot.py stats_sample_count \
    --otu_count_tsv "$data_dir/output_unoise3_16s/unoise3_zotu.biom" \
    --metadata "$data_dir/metadata/bulk_16s.tsv" \
    --otu_taxonomy_tsv "$data_dir/output_unoise3_16s/unoise3_zotu_rrndb_processed.tsv" \
    --fig_dir "$data_dir/figs/alpha_diversity_16s" \
    --tax_levels "otu" "genus"

plot.py stats_sample_count \
    --otu_count_tsv "$data_dir/output_unoise3_16s/unoise3_zotu.biom" \
    --metadata "$data_dir/metadata/bulk_16s.tsv" \
    --otu_taxonomy_tsv "$data_dir/output_unoise3_16s/unoise3_zotu_rrndb_processed.tsv" \
    --fig_dir "$data_dir/figs/alpha_diversity_rarefying_16s" \
    --tax_levels "otu" "genus" \
    --rarefying_value 500 \
    --rarefying_repeat 50

# Command 1: Euclidean distance with log10, no annotation
plot_dm.py \
    --otu_count_tsv "$data_dir/output_unoise3_16s/unoise3_zotu.biom" \
    --metadata "$data_dir/metadata/bulk_16s.tsv" \
    --otu_taxonomy_tsv "$data_dir/output_unoise3_16s/unoise3_zotu_rrndb_processed.tsv" \
    --fig_dir "$data_dir/figs/beta_diversity_16s" \
    --tax_levels "otu" "genus" \
    --distance "euclid" --log10 --transform "relative" \
    --hue "rep_group"

# Command 2: Bray-Curtis distance, no annotation
plot_dm.py \
    --otu_count_tsv "$data_dir/output_unoise3_16s/unoise3_zotu.biom" \
    --metadata "$data_dir/metadata/bulk_16s.tsv" \
    --otu_taxonomy_tsv "$data_dir/output_unoise3_16s/unoise3_zotu_rrndb_processed.tsv" \
    --fig_dir "$data_dir/figs/beta_diversity_16s" \
    --tax_levels "otu" "genus" \
    --distance "braycurtis" \
    --hue "rep_group"

# Command 3: Euclidean distance with log10 and annotation
plot_dm.py \
    --otu_count_tsv "$data_dir/output_unoise3_16s/unoise3_zotu.biom" \
    --metadata "$data_dir/metadata/bulk_16s.tsv" \
    --otu_taxonomy_tsv "$data_dir/output_unoise3_16s/unoise3_zotu_rrndb_processed.tsv" \
    --fig_dir "$data_dir/figs/beta_diversity_annot_16s" \
    --tax_levels "otu" "genus" \
    --distance "euclid" --log10 --transform "relative" \
    --hue "rep_group" \
    --annotate_dots

# Command 4: Bray-Curtis distance with annotation
plot_dm.py \
    --otu_count_tsv "$data_dir/output_unoise3_16s/unoise3_zotu.biom" \
    --metadata "$data_dir/metadata/bulk_16s.tsv" \
    --otu_taxonomy_tsv "$data_dir/output_unoise3_16s/unoise3_zotu_rrndb_processed.tsv" \
    --fig_dir "$data_dir/figs/beta_diversity_annot_16s" \
    --tax_levels "otu" "genus" \
    --distance "braycurtis" \
    --hue "rep_group" \
    --annotate_dots
