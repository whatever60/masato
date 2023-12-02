# Bacteria 16S rRNA-seq analysis pipeline

## Requirement

Successful execution of this pipeline requires the following software:

- [vsearch](https://github.com/torognes/vsearch)
- [RDP classifier](https://github.com/rdpstaff/classifier)
- [seqtk](https://github.com/lh3/seqtk)
- [seqkit](https://github.com/shenwei356/seqkit)

If you read the code, you'll notice `usearch` could also be used by changing `backend` argument. However, I recommend using `vsearch` because it is open source, though their implementations and CLI commands may vary and `vsearch` could be slower and miss some features of `usearch`

## Usage 1: UNOISE3 denoising and RDP-based taxonomy

### Step 1: Prepare demultiplexed fastq files and a metadata file.

The input directory should contains all fastq files of interest. In general, sequencing results of samples from the same study could all be included in the same directory, and samples across study should be put into different directory and analyzed separately. File names should look like `<sample_name>_R1.fastq` and `<sample_name>_R2.fastq`.

Theaa metadata file should be a tsv text file with the these columns:

- `sample`: it contains sample names. All entries of this column should have corresponding fastq files in the input directory.
- `rep_group`: it specifies replication that should be eventually averaged. Intuitively, the entries in this column will appear in figures you want to present or include in a paper.
- `sample_group`: it is just for the convenience of plotting. In the final plotting step, samples with the same `sample_group` will be plotted in the same panel.
- `reference`: it is for subsequent group-wise subsampling (see Step 5). It points to the sample name whose sequencing depth will be used for subsampling.
- `spike_in`: optional, in the format of `<taxon_level>;<taxon_name>` such as `genus;Sporosarcina`. It specifies spike-in taxa that should be removed when getting the abundance and plotting.
- `sample_weight`: optional. It is used for normaliation when getting absolute abundance. Note that a normalization factor will be calculated based on the product of spike in count and sample weight, so depending on study design, one normalization might suffice and you don't need both `sample_weight` and `spike_in`.

Note that the order of samples in the metadata file matters since it determines the order of samples in the figures we will plot. 

### Step 2: Populate metadata file with sequencing depth.

```shell
usearch_workflow.py add_depth_to_metadata --fastq_dir <path_to_fastq_dir> --metadata <path_to_metadata>
```

The metadata is modified in-place by adding an integer column `read_count`.

### Step 3: Merge paired-end reads.

```shell
usearch_workflow.py merge_pairs -i <path_to_fastq_dir> -o <output_dir> -t <threads>
```

Pool all reads across all samples and merge paired-end reads. A merged fastq file and fasta file will be generated. Sequences are relabeld to `sample=<sample_name> <read_index> <original_name>`.

### Step 4: Unique high quality non-singleton sequence database construction.

```shell
usearch_workflow.py db_construct -i <path_to_merged_fastq> -t <threads>
```

Output files will be in the same directory as the input fastq file. The file that contains the database ends with `.uniq.fa`.

Though pooling reads across all samples is recommended to boost abundance, I still perform subsampling after database construction. This is because subsampling could cause spurious abundance increase and thus add more noise to the database.

### Step 5: Denoising.

Denoise database with `UNOISE3` and assign all reads (from which the database is derived) to a ZOTU. `UPARSE` serves similar purpose but should not be used anymore.

The database is the `.uniq.fa` file generated in `Step 4`, or, if you didn't follow the previous step, any fasta file with unique high quality sequences will do. Note that in the sequence labels there must be a `size` field. There is no need to exclude low abundance (i.e., low `size` value) sequences in this field, since subsequent denoising will take care of that.

```shell
usearch_workflow.py cluster_unoise3 -i <path_to_fastq> -d <path_to_db> -o <output_dir> -t <threads>
```

### Step 6: Taxonomy assignment for OTUs.

Use `rdp_classifier` (RDP Naive Bayesian Classifier) to assign taxonomy to each OTU.

```shell
# run RDP classifier
java -Xmx1g -jar <path_to_classifier.jar> classify -t <path_to_16srrna_rRNAClassifier.properties> -c 0.8 -o <output_path> <path_to_fasta>
# parse RDP classifier output
process_rrndb.py -i <rdp_classifier_output> -o <processed_output>
```

The output of `RDP classifier` is a tab-delimited text file with no header. Its columns are OTU label, taxonomy (start with `Root`), taxonomy rank (start with `rootrank`), confidence score, taxonomy of next level, taxonomy rank of next level, confidence score of next level, and so on.

We process the output into a tsv file with columns `otu`, `domain`, `domain_p`, ..., `species`, `species_p`.


### Step 7: Generate OTU and taxonomy abundance table.

Aggregate OTU count table into sample-level OTU abundance (both relative and absolute) and group-level taxonomy relative abundance, i.e., samples from the same `rep_group` are averaged. Spike-in taxa are removed in this step. We can certainly group OTUs and samples arbitrarily, but I think interesting ones are:
- OTU x sample absolute abundance
- OTU x sample relative abundance, 
- OTU/genus/family/... x group relative abundance

```shell
get_abundance.py \
    -i <path_to_otu_count_table> \
    -m <path_to_metadata> \
    -t <path_to_processed_rrndb_taxonomy> \
    -o <output_directory> \
    -l <a list of taxonomy levels of interest> \
    -s <key_in_metadata_that_specify_spikein> \
    -r <key_in_metadata_that_specify_replication_group> \
    -w <key_in_metadata_that_specify_sample_weight>
    -a <threshold_to_consider_a_taxon_as_rare> \
    --keep_others \
    --keep_unknown
```

### Step 8 (last step, well done): Make figures.

Plot figures that show relative abundance, one for each taxonomy level you specified.

Three choices are available: stacked bar plot, heatmap, and heatmap with `y=log10(x+1e-4)` transformation. The plots will have several axes, and each axes will show the relative abundance of samples from the same `sample_group`. Samples that are also in the same `rep_group` are averaged and not shown individually in these figures. The axis tick labels will be the `rep_group` names.

```shell
plot.py abundance_group \
    [stacked_bar|heatmap|heatmap_log10|all] \
    -i <directory_of_sample_level_otu_relative_abundance> \
    -m <path_to_metadata> \
    -f <directory_to_store_figures> \
    -t <a_list_of_taxonomy_levels_of_interest>
    -s <key_in_metadata_that_specify_sample_group> \
    -r <key_in_metadata_that_specify_replication_group> \
```

Plot barplot figures that show per `rep_group` statistics, one for each statistics and each taxonomy level you specified. In most cases 

```shell
plot.py stats_sample \
    -i <directory_of_smaple_level_otu_relative_abundance> \
    -m <path_to_metadata> \
    -f <directory_to_store_figures> \
    -t <a_list_of_taxonomy_levels_of_interest> \
    -s <key_in_metadata_that_specify_sample_group> \
    -r <key_in_metadata_that_specify_replication_group> \
```

## Usage 2: SINTAX-based taxonomy

Use `SINTAX` (simple non-Bayesian taxonomy classifier) to assign taxonomy to each sequence.

```shell
# run SINTAX algorithm
usearch_workflow.py tax_sintax -i <path_to_fasta> -d <path_to_db> -o <output_dir> -t <threads>
# parse SINTAX output
analysis_sintax
```

The output of `SINTAX` is a tab-delimited text file with no header. Its columns are read label, taxonomy with confidence score, strand, and taxonomy without confidence score.

## Usage 3: Subsampling (aka. rarefying)

```shell
usearch_workflow.py subsample -i <path_to_fastq_dir> -m <path_to_metadata> -o <fastq_output_dir> -om <subsample_metadata_path> -n <subsample_n_times> --mode [reference|min]
```

Subsample a few times for each pairs of fastq files in the input directory. For a fastq file named `<sample_name>_R1/2.fastq`, the resulting subsampled files would be named as `<original_sample_name>-subsampled-<k>_R1/2.fastq`. Thus the sample names of the new fastq files are changed to `<original_sample_name>-subsampled-<k>`.

A output metadata will be generated for the new samples. `sample` now has the new sample names, and a new column `original_sample` is added to indicate the original sample name. The `read_count` column now have the read count of the new samples. All other columns in the original metadata are copied to all newly subsampled samples from their original samples.

There are two modes of subsampling, `reference` and `min`. In `reference` mode, the subsampling depth is the read count of the sample specified by the `reference` column in the original metadata. In case where a sample have lower sequencing depth than its reference, these samples will only be subsampled once by essentially copying over, no matter of the `-n` in the command. In `min` mode, the subsampling depth is the minimum read count of all samples in the original metadata minus 1.

The `reference` mode might be confusing at first. It is implemented because in one study we have paired samples obtained from environment direcly and plate scraping of culture plates inoculated with biomass from the same environment. Theoretically the former is more likely to have lower sequencing depth, since there is less material to start with, so we think the environment samples do not need subsampling and all other samples should be subsampled to the depth of the corresponding environment sample.

To combine subsampling with Usage 1, after subsampling:
1. Perform paired-end merging for subsampled reads as specified in Usage 1 Step 3.
2. `cat` the two merged big `.fastq` into an even bigger one.
3. Do denoising as specified in Usage 1 Step 5, but on the big `.fastq` file. Since each sample should have a different name, they could be directly distinguishable by the following scripts.
4. Get abundance tables and plot figures separately for original samples and subsampled ones with their own metadata files. This is the common use case, but it is also possible that you want to combine the abundance tables or collectively plot things. Write your some custom scripts for that.

QUESTION: should we remove spike-in reads before subsampling?

## Usage 3: QIIME2-based taxonomy

TODO

## A note on taxonomy assignment algorithms

Two classic algorithms are widely used for taxonomy assignment: RDP and SINTAX. RDP is a naive Bayesian classifier, while SINTAX is a k-mer based algorithm. The former is more accurate and the latter is faster.

RDP is implemented in multiple softwares, including the original RDP classifier, USEARCH, rRDP R package, and mothur, as well as online tools like [rrnDB](https://rrndb.umms.med.umich.edu/). Unfortunately, it is not implemented in `vsearch`, so one more dependency :/

SINTAX is available in USEARCH and its open source alternative, VSEARCH.
