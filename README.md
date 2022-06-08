\**These instructions have only been tested on a Mac, although directions for Linux machines should be very similar (you will just need to install Conda for Linux instead of Mac). Let me know if you are having issues setting everything up.

# Conda
Before anything, you will need to get Conda running on your machine. Conda is a great tool that allows you to create environments where you only install certain packages and dendencies to help eliminate package dependency and versioning nightmares on your entire machine
- A great Conda [cheat sheet](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf "cheat sheet") I refer to often
- If Conda becomes too difficult to use, we can migrate to other options such as Docker

## Install Conda
- Navigate to [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html "https://docs.conda.io/en/latest/miniconda.html")
- Click 'Miniconda3 MaxOSX 64-bit bash' and download the installer file to an easily accessible location on your computer (e.g. desktop)
- Open Terminal and navigate to the directory containing the installer file you just downloaded
	- `cd ~/Desktop`
- Run the installer
	- `bash Miniconda3-latest-MacOSX-x86_64.sh`
- If after installing Conda you would like to remove '(base)' that now appears on the command line, you can run `conda config --set auto_activate_base false`
- If you hate Conda, here is a [link](https://stackoverflow.com/questions/42182706/how-to-uninstall-anaconda-completely-from-macos "link") explaining how to uninstall Conda (you will use `miniconda3` instead of `anaconda3` which is used in the examples)

## Create environment to run `cfseq`
- On command line, run `conda create --name cfseq python=3.7` to create a new environment to run `cfseq`
- Type `y` and press 'Return' to proceed installing listed packages
- Run `conda activate cfseq` to enter new environment
- Run `conda install -c bioconda bismark=0.23.1 fastqc=0.11.9 bowtie2=2.2.5 bedtools=2.30.0 simplesam=0.1.3.1 samtools=1.9` to install all the necessary packages/versions to run this pipeline
- Type `y` and press 'Return' to proceed
- Now you should see `(cfseq)` on the command line meaning you are in the `cfseq` conda environment and are ready to proceed
	- Enter a Conda environment `conda activate environment_name`
	- Exit a Conda environment `conda deactivate`

# Reference genome
Now, you will need to download the necessary reference DNA sequences and create the necessary reference index files so you can align your methylation sequencing reads.

## Grab reference sequences

### Chromosome-by-chromosome approach
Here, I only download the reference sequences for the chromsomes which contain target CpGs (to save space on my computer). If hard drive space is not an issue for you and/or if you expect to add many more target CpGs on many chromosomes later on to your panel, you may consider downloading and indexing the entire genome to prevent the need to repeat this step when expand your panel to include targets on new chromosomes (see **'Entire genome approach'**).
- Download the reference sequence (FASTA file, e.g. chr7.fa.gz) for each chromosome containing a target CpG from UCSC ([https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/ "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/"))
	- Click on desired file and download to your computer OR right-click desired file and select 'Save Link As...' and save to your computer
- Combine all separate chromosome sequence files into one
	- On the command line, navigate to the directory containing your FASTA reference files and concatenate all your chromosome FASTA files into one file `cat chr6.fa.gz chr7.fa.gz chr10.fa.gz chr17.fa.gz > reference_chromosomes.fa.gz` and delete the individual `.fa.gz` files (if desired).

### Entire genome approach
- Download the reference sequence for the entire genome from UCSC ([https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/ "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/"))
	- Click on `hg19.fa.gz` and download to your computer OR right-click file and select 'Save Link As...' and save to your computer

## Create Bismark methylation reference genome index files
- On command line, run: `bismark_genome_preparation  /path/to/reference/genome/directory`
- Bismark will recognize any \*.fasta or \*.fa file here and make a reference index from it
- This step can take > 30 minutes depending on the size of your reference genome

# Download pipeline from Github
The code for this pipeline is hosted on [Github](https://github.com/Inherent-Biosciences/cf-seq "Github") which will allow us to stay up to date on the latest version of the pipeline

## Download via command line
- Make sure you have 'git' installed on your computer ([https://git-scm.com/download/mac](https://git-scm.com/download/mac "https://git-scm.com/download/mac"))
- Login to your Github account:
	- `git config user.name "your_username"`
	- `git config user.email "your_email_address@example.com"`
- Navigate to directory where you want to download the pipeline and run:
	- `git clone https://github.com/Inherent-Biosciences/cf-seq`

## Download from internet browser
- Go to [https://github.com/Inherent-Biosciences/cf-seq](https://github.com/Inherent-Biosciences/cf-seq "https://github.com/Inherent-Biosciences/cf-seq")
- 'Code' > 'Download ZIP'
- Double-click zipped file to unzip code directory

# Run pipeline
For now, the pipeline is deployed as a bash script which runs all the different parts of the analysis. In the future, we can encapsulate this into an all-in-one package with argument flags like any other typical command line program, but for the moment, I think this does the job. 

To run the pipeline, all you will need are:
- reference files you made previously
- FASTQ files to analyze
- .bed file containing the information about each CpG you want to analyze.

## .bed targets file
Necessary characteristics of your targets .bed file
- header line
- tab-separated
- the columns need to be
	- chromosome (in `chr6` format)
	- start position
	- end position (start position or (start position + 1))
	- name of CpG (name as you see fit)

Example: `targets.bed`
```
Chrom	Start	Stop	Name
chr6	25962987	25962988	cg08501292-PCR7
chr6	25963022	25963023	cg11571304-PCR7
chr17	80255457	80255458	cg21435684-PCR11
chr17	80255510	80255511	cg12593223-PCR11
chr17	80255734	80255735	cg26663796-PCR11
chr10	45719880	45719881	cg20447730-PCR12
chr10	45719956	45719957	cg01782066-PCR12
chr7	45018757	45018758	cg12416569-PCR31
chr7	45018849	45018850	cg10673833-PCR31
```

## Modify arguments as needed in `cf_dna_pipeline.sh`
Open `cf_dna_pipeline.sh` in a text editor, make the necessary changes, and save the file.
```
########### FASTQ FILE NAME INFO ############
## for the fastq files, you will need to split them up into a 'prefix' and 'suffix' (prefix + suffix == full file name)
#### e.g. fastq file name: BYU-B1_S1_L001_R1_001.fastq.gz

#### so an example 'prefix' could be 'BYU-B1_S1'
####### prefixes are used as identifiers for each sample and will be used to help name downstream analysis files for each sample

#### and an example 'suffix' could be '_L001_R1_001.fastq.gz'
####### suffixes are used to identify the forward and reverse reads (in the suffix example 'R1' refers to forward read and 'R2' refers to reverse read)
#############################################

fastq_dir="/Users/ryan/Desktop/cf_dna_data/first_run/fastqs" # full path to the directory containing fastq files; e.g. "/Users/ryan/Desktop/cf_dna_data/first_run/fastqs"
sample_prefixes=(BYU-B1_S1 BYU-B2_S3 BYU-B3_S5) # unique, identifiable prefixes for each sample you want to analyze e.g. (BYU-B1_S1 BYU-B2_S3 BYU-B3_S5)
forward_read_suffix="_L001_R1_001.fastq.gz" # e.g. "_L001_R1_001.fastq.gz"
reverse_read_suffix="_L001_R2_001.fastq.gz" # e.g. "_L001_R2_001.fastq.gz"

bismark_cores=3 # number of cores to use when running bismark, 3 works well on my 1.1 GHz Quad-Core Intel Core i5, 16GB RAM, MacBook Air. If you have a beefier machine, feel free to increase this argument, but I would suggest increasing by increments of 1 at a time
bismark_reference_dir="/Users/ryan/Desktop/cf_dna_data/reference" # full path to directory containing reference genome and Bismark index files; e.g. "/Users/ryan/Desktop/cf_dna_data/reference"

python_script_path="/Users/ryan/Desktop/cf-seq/cfseq.py" # full path to cfseq.py script; e.g. "/Users/ryan/Desktop/cf-seq/cfseq.py"
targets_file_path="/Users/ryan/Desktop/cf-seq/targets.bed" # full path to target CpG bed file; e.g. "/Users/ryan/Desktop/cf-seq/targets.bed"
```

## Execute pipeline
In the command line, navigate to directory containing file 'cf_dna_pipeline.sh' and run: `bash cf_dna_pipeline.sh`

## Pipeline Outputs
For each sample:
- FastQC .html report file for forward and reverse read (e.g. BYU-B1_S1_L001_R1_001_fastqc.html and BYU-B1_S1_L001_R2_001_fastqc.html)
- Bismark report file (e.g. BYU-B1_S1_L001_R1_001_bismark_bt2_PE_report.txt)
- .sam file containing alignment info for each paired-end read (e.g. BYU-B1_S1.sam)
- Output file with methylation state of each read (e.g. BYU-B1_S1.sam-methylation_results.tsv)
- For each chromosome, an output file with methylation state of reads covering a target CpG on that chromosome (e.g. BYU-B1_S1_chr6_results.tsv, etc.)

