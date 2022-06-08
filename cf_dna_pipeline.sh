########### FASTQ FILE NAME INFO ############
## for the fastq files, you will need to split them up into a 'prefix' and 'suffix' (prefix + suffix == full file name)
#### e.g. fastq file name: BYU-B1_S1_L001_R1_001.fastq.gz

#### so an example 'prefix' could be 'BYU-B1_S1'
####### prefixes are used as identifiers for each sample and will be used to help name downstream analysis files for each sample

#### and an example 'suffix' could be '_L001_R1_001.fastq.gz'
####### suffixes are used to identify the forward and reverse reads (in the suffix example 'R1' refers to forward read and 'R2' refers to reverse read)
#############################################


################################################
#### ARGUMENTS - PLEASE MODIFY AS NECESSARY ####
################################################
fastq_dir="/Users/ryan/Desktop/cf_dna_data/first_run/fastqs" # full path to the directory containing fastq files; e.g. "/Users/ryan/Desktop/cf_dna_data/first_run/fastqs"
# sample_prefixes=(BYU-B1_S1 BYU-B2_S3 BYU-B3_S5 BYU-B4_S7 BYU-B5_S9 BYU-B6_S10 BYU-G1_S11 BYU-G2_S12 BYU-G3_S2 BYU-G4_S4 BYU-G5_S6 BYU-G6_S8)
sample_prefixes=(BYU-B1_S1) # unique, identifiable prefixes for each sample you want to analyze e.g. (BYU-B1_S1 BYU-B2_S3 BYU-B3_S5)
forward_read_suffix="_L001_R1_001.fastq.gz" # e.g. "_L001_R1_001.fastq.gz"
reverse_read_suffix="_L001_R2_001.fastq.gz" # e.g. "_L001_R2_001.fastq.gz"

bismark_cores=3 # number of cores to use when running bismark, 3 works well on my 1.1 GHz Quad-Core Intel Core i5, 16GB RAM, MacBook Air. If you have a beefier machine, feel free to increase this argument, but I would suggest increasing by increments of 1 at a time
bismark_reference_dir="/Users/ryan/Desktop/cf_dna_data/reference" # full path to directory containing reference genome and Bismark index files; e.g. "/Users/ryan/Desktop/cf_dna_data/reference"

python_script_path="/Users/ryan/Desktop/cf-seq/cfseq.py" # full path to cfseq.py script; e.g. "/Users/ryan/Desktop/cf-seq/cfseq.py"
targets_file_path="/Users/ryan/Desktop/cf-seq/targets.bed" # full path to target CpG bed file; e.g. "/Users/ryan/Desktop/cf-seq/targets.bed"

#################################################################
#### END OF ARGUMENTS; NO NEED TO CHANGE ANYTHING BELOW THIS ####
#################################################################

# making sure last char of fastq_dir is '/'
# if [ ${fastq_dir: -1} != "/" ]
# then
#     fastq_dir="${fastq_dir}/"
# fi

cd $fastq_dir

for i in ${sample_prefixes[*]}
do 
	forward_read=${i}${forward_read_suffix}
	reverse_read=${i}${reverse_read_suffix}

    ## run fastqc
	echo "# running fastQC on ${forward_read} and ${reverse_read}"
	echo ""
    fastqc ${forward_read} ${reverse_read} -o ${fastq_dir}
    # remove fastqc .zip files
    rm *_fastqc.zip


	## align with bismark
	echo "# aligning ${forward_read} and ${reverse_read}"
	echo ""
	bismark --bowtie2 --multicore $bismark_cores --genome_folder $bismark_reference_dir -1 $forward_read -2 $reverse_read


	## convert .bam to .sam
	echo "# converting .bam output file to .sam file for further analysis"
	echo ""
	samtools view -h ${i}*.bam -o ${i}.sam


	## deleting .bam file
	rm ${i}*.bam


	## running molecule seq
	echo "# analyzing molecule-specific methylation on ${i}.sam"
	echo ""
	python $python_script_path --sam ${i}.sam --targets $targets_file_path


    ## separate output into separate files based on chromosome
    ## WARNING - this is very much 'hard-coded' at this point and will need to be changed to accomodate changes to assay
    awk -F"\t" '{print $1"\t"$2"\t"$3}' "${i}.sam-methylation_results.tsv" | grep -ve 'NA\tNA' > "${i}_chr6_results.tsv"
    awk -F"\t" '{print $1"\t"$4"\t"$5}' "${i}.sam-methylation_results.tsv" | grep -ve 'NA\tNA' > "${i}_chr7_results.tsv"
    awk -F"\t" '{print $1"\t"$6"\t"$7}' "${i}.sam-methylation_results.tsv" | grep -ve 'NA\tNA' > "${i}_chr10_results.tsv"
    awk -F"\t" '{print $1"\t"$8"\t"$9"\t"$10}' "${i}.sam-methylation_results.tsv" | grep -ve 'NA\tNA\tNA' > "${i}_chr17_results.tsv"

	#break
done
