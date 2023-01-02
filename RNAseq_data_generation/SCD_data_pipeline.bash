#!/bin/bash
# The interpreter used to execute the script

# your job commands below
#for Q6
LP='XXX'
module load Bioinformatics
module load hisat2
module load samtools
module load bedtools2/2.29.2
module load sratoolkit/2.9.6
module load homer/4.10
module load stringtie


name_list="SRR14479828 SRR14479829 SRR14479830"

for input_name in $name_list
do 
	mkdir -p $LP/Analysis_${input_name}
	# Step_1 download data
	fastq-dump --outdir $LP/Analysis_${input_name} --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip ${input_name}
	mkdir $LP/Analysis_${input_name}/QC

	# Step_2 QC
	fastqc -o $LP/Analysis_${input_name}/QC $LP/Analysis_${input_name}/${input_name}_pass.fastq.gz 

	# Step_3 remove barcode and adapter 
	TrimmomaticPE $LP/Analysis_${input_name}/${input_name}_pass_1.fastq.gz  $LP/Analysis_${input_name}/${input_name}_pass_2.fastq.gz\
				  $LP/Analysis_${input_name}/${input_name}_pass_1.trimmed.fastq $LP/Analysis_${input_name}/${input_name}_pass_1un.trimmed.fastq \
				  $LP/Analysis_${input_name}/${input_name}_pass_2.trimmed.fastq $LP/Analysis_${input_name}/${input_name}_pass_2un.trimmed.fastq \
				  -phred33 ILLUMINACLIP:$LP/Analysis_${input_name}/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36



	# Step_4 map to reference genome
	hisat2 -q -x $LP/index_R64/r64_tran/genome_tran -1 $LP/Analysis_${input_name}/${input_name}_pass_1.trimmed.fastq -2 $LP/Analysis_${input_name}/${input_name}_pass_2.trimmed.fastq --max-intronlen 2500 --rna-strandness FR -S $LP/Analysis_${input_name}/mapped_result.sam



	# Step_5 sort, filter and index sam output, I remove reads with mapping quality lower than 20

	samtools view -bS -F 4 -q 20 $LP/Analysis_${input_name}/mapped_result.sam | samtools sort -o $LP/Analysis_${input_name}/mapped_result_sorted.bam

	samtools index $LP/Analysis_${input_name}/mapped_result_sorted.bam

	# Step_6 convert bam into bed file 

	bedtools bamtobed -i $LP/Analysis_${input_name}/mapped_result_sorted.bam > $LP/Analysis_${input_name}/mapped_result_raw.bed

	# get rid of reads mapped, two ends of which are longer than 5kb and remove reads mapped to mitochondria
	awk '{if($3-$2 <= 5000) print}' $LP/Analysis_${input_name}/mapped_result_raw.bed | grep -v 'Mito' > $LP/Analysis_${input_name}/mapped_result.bed



	mkdir -p $LP/Analysis_${input_name}/expression
	########################################################################################################################################################

	# Sliding window method
	# genic expression
	bedtools coverage -a $LP/Final_bed_file/SW_genic_w20_s20.bed -b $LP/Analysis_${input_name}/mapped_result.bed -s -sorted> $LP/Analysis_${input_name}/expression/Expression_SW_stranded_genic_w20_s20
	# intergenic expression
	bedtools coverage -a $LP/Final_bed_file/SW_intergenic_N1_w20_s20.bed -b $LP/Analysis_${input_name}/mapped_result.bed -s -sorted> $LP//Analysis_${input_name}/expression/Expression_SW_stranded_intergenic_N1_w20_s20
	bedtools coverage -a $LP/Final_bed_file/SW_intergenic_N2_w20_s20.bed -b $LP/Analysis_${input_name}/mapped_result.bed -s -sorted> $LP//Analysis_${input_name}/expression/Expression_SW_stranded_intergenic_N2_w20_s20
	bedtools coverage -a $LP/Final_bed_file/SW_intergenic_N3_w20_s20.bed -b $LP/Analysis_${input_name}/mapped_result.bed -s -sorted> $LP//Analysis_${input_name}/expression/Expression_SW_stranded_intergenic_N3_w20_s20
	bedtools coverage -a $LP/Final_bed_file/SW_intergenic_N4_w20_s20.bed -b $LP/Analysis_${input_name}/mapped_result.bed -s -sorted> $LP//Analysis_${input_name}/expression/Expression_SW_stranded_intergenic_N4_w20_s20
	bedtools coverage -a $LP/Final_bed_file/SW_intergenic_N5_w20_s20.bed -b $LP/Analysis_${input_name}/mapped_result.bed -s -sorted> $LP//Analysis_${input_name}/expression/Expression_SW_stranded_intergenic_N5_w20_s20
	bedtools coverage -a $LP/Final_bed_file/SW_intergenic_N6_w20_s20.bed -b $LP/Analysis_${input_name}/mapped_result.bed -s -sorted> $LP//Analysis_${input_name}/expression/Expression_SW_stranded_intergenic_N6_w20_s20
	bedtools coverage -a $LP/Final_bed_file/SW_intergenic_N7_w20_s20.bed -b $LP/Analysis_${input_name}/mapped_result.bed -s -sorted> $LP//Analysis_${input_name}/expression/Expression_SW_stranded_intergenic_N7_w20_s20
	bedtools coverage -a $LP/Final_bed_file/SW_intergenic_N8_w20_s20.bed -b $LP/Analysis_${input_name}/mapped_result.bed -s -sorted> $LP//Analysis_${input_name}/expression/Expression_SW_stranded_intergenic_N8_w20_s20
	bedtools coverage -a $LP/Final_bed_file/SW_intergenic_N9_w20_s20.bed -b $LP/Analysis_${input_name}/mapped_result.bed -s -sorted> $LP//Analysis_${input_name}/expression/Expression_SW_stranded_intergenic_N9_w20_s20
	bedtools coverage -a $LP/Final_bed_file/SW_intergenic_N10_w20_s20.bed -b $LP/Analysis_${input_name}/mapped_result.bed -s -sorted> $LP//Analysis_${input_name}/expression/Expression_SW_stranded_intergenic_N10_w20_s20
	bedtools coverage -a $LP/Final_bed_file/SW_intergenic_N11_w20_s20.bed -b $LP/Analysis_${input_name}/mapped_result.bed -s -sorted> $LP//Analysis_${input_name}/expression/Expression_SW_stranded_intergenic_N11_w20_s20
	bedtools coverage -a $LP/Final_bed_file/SW_intergenic_N12_w20_s20.bed -b $LP/Analysis_${input_name}/mapped_result.bed -s -sorted> $LP//Analysis_${input_name}/expression/Expression_SW_stranded_intergenic_N12_w20_s20
	bedtools coverage -a $LP/Final_bed_file/SW_intergenic_N13_w20_s20.bed -b $LP/Analysis_${input_name}/mapped_result.bed -s -sorted> $LP//Analysis_${input_name}/expression/Expression_SW_stranded_intergenic_N13_w20_s20
	bedtools coverage -a $LP/Final_bed_file/SW_intergenic_N14_w20_s20.bed -b $LP/Analysis_${input_name}/mapped_result.bed -s -sorted> $LP//Analysis_${input_name}/expression/Expression_SW_stranded_intergenic_N14_w20_s20


	# do not use sliding windows

	# genic expression
	stringtie $LP/Analysis_${input_name}/mapped_result_sorted.bam -G $LP/Saccharomyces_cerevisiae.R64-1-1.104.gtf -m 100 -t -j 3 -o $LP/Analysis_${input_name}/expression/With_annotation.gtf -A $LP/Analysis_${input_name}/expression/Expression_with_annotation.tab
	# Step_10 check reads number 
	wc -l $LP/Analysis_${input_name}/mapped_result.bed
########################################################################################################################################################
done