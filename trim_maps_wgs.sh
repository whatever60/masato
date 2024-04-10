input_fastq_dir=/mnt/c/aws_data/20240328_liyuan_maps_wgs_random_hexamer/fastq
# get parent directory
data_dir=$(dirname $input_fastq_dir)
output_dir=$data_dir/fastq_qc
output_bam_dir=$data_dir/map
mkdir -p $output_dir
mkdir -p $output_bam_dir

# maps random hexamer test
mkdir -p $output_dir/round1
./trim.py \
  --mode maps_rand_hex_test \
  -i $input_fastq_dir \
  -o $output_dir/round1/merged.fq.gz \
  -fb data/maps_seq/BC1.fasta

# first round
mkdir -p $output_dir/round1
./trim.py \
  --mode maps_round1 \
  -i $input_fastq_dir \
  -o $output_dir/round1/merged.fq.gz \
  -fb data/maps_seq/BC1.fasta

# second round
mkdir -p $output_dir/round2
./trim.py \
  --mode maps_round2 \
  -i $output_dir/round1/merged.fq.gz \
  -o $output_dir/round2/merged.fq.gz \
  -fb data/maps_seq/BC2.fasta

# third round
mkdir -p $output_dir/round3
./trim.py \
  --mode maps_round3 \
  -i $output_dir/round2/merged.fq.gz \
  -o $output_dir/round3/merged.fq.gz \
  -fb data/maps_seq/BC3.fasta

zcat $output_dir/round3/merged.fq.gz | awk ' 
BEGIN {
    OFS="\n" # Set the Output Field Separator to a newline character
}
{
    if (NR % 4 == 1) { # Process every 4th line starting with the first
        split($0, a, " "); # Split the line by spaces into array a
        split(a[1], b, "="); # Split the first element of a by "=" into array b
        print "@"b[2]"_"a[2]"_"a[3] # Print the modified identifier
    } else {
        print $0 # Print lines that are not processed by the if statement
    }
}
' | gzip > $output_dir/round3/merged_rename.fq.gz

bwa mem \
    -t 16 $data_dir/ref_genomes/combined.fasta \
    $output_dir/round3/merged_rename.fq.gz \
    | samtools sort -o $output_bam_dir/default.sorted.bam -
sambamba index -t 16 $output_bam_dir/default.sorted.bam
cmd='String[] parts = record.getReadName().split(\"_\"); record.setAttribute(\"CB\", parts[0] + \"_\" + parts[1] + \"_\" + parts[2] + \"_\" + parts[3]); return record;'
jvarkit samjdk \
    -e "$cmd" \
    $output_bam_dir/default.sorted.bam \
    --samoutputformat BAM \
    -o $output_bam_dir/default.tagged.sorted.bam
sambamba index -t 16 $output_bam_dir/default.tagged.sorted.bam
umi_tools dedup \
    --ignore-umi \
    --cell-tag CB \
    --extract-umi-method tag \
    --assigned-status-tag XS \
    -I $output_bam_dir/default.tagged.sorted.bam \
    -S $output_bam_dir/default.tagged.dedup.sorted.bam
sambamba index -t 16 $output_bam_dir/default.tagged.dedup.sorted.bam
