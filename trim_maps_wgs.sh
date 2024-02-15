input_fastq_dir=/mnt/c/aws_data/20240214_liyuan_maps_wgs/fastq_raw
# get parent directory
data_dir=$(dirname $input_fastq_dir)
output_dir=$data_dir/fastq_qc

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
