#!/bin/bash -l

#################################################
## Shell script to run steps 2-5 together      ##
#################################################

####Array script for fastp trimming#####
####Modified from Laxman Adhikari and Sandesh Shrestha#####
##https://docs.rc.fas.harvard.edu/kb/submitting-large-numbers-of-jobs/####

#SBATCH --job-name=bio_process
#SBATCH --time=00-23:59:00 #Use the from DD-HH:MM:SS
#SBATCH --mem-per-cpu=36G #Memory per core, use --mem= form memory per node
#SBATCH --mail-type=ALL   # same as =BEGIN,FAIL,END 
#SBATCH --output="Step_2_5_%A_%a.out"
#SBATCH --error="Step_2_5_%A_%a.err"

ulimit -u 8192
module load HISAT2
module load SAMtools/1.10-GCC-8.3.0
module list

INPUTDIR=$1 #get working directory for files

mkdir -p ${INPUTDIR}/02_Trimmed ${INPUTDIR}/out_error/02_Trimmed
mkdir -p ${INPUTDIR}/03_Aligned ${INPUTDIR}/out_error/03_Aligned 
mkdir -p ${INPUTDIR}/04_Filtered ${INPUTDIR}/out_error/04_Filtered
mkdir -p ${INPUTDIR}/05_Indexed ${INPUTDIR}/out_error/05_Indexed



FILES=($(ls ${INPUTDIR}/01_Samples/*R1.fastq.gz))

FILENAME=${FILES[$SLURM_ARRAY_TASK_ID]} #this gets file name from outside array
echo ${FILENAME}
FILE_PROCESS=$(basename ${FILENAME} "_R1.fastq.gz") #remove _R1.fq.gz ending

echo "Current working directory: `pwd`"

echo "Doing JOB  on |" $(date) " |  $INPUTDIR "
echo "" 
echo $FILE_PROCESS

######### TRIMMING #########################
/homes/jcrain/Software/fastp/fastp -i ${INPUTDIR}/01_Samples/${FILE_PROCESS}_R1.fastq.gz -I ${INPUTDIR}/01_Samples/${FILE_PROCESS}_R2.fastq.gz -o ${INPUTDIR}/02_Trimmed/${FILE_PROCESS}_R1.fastq.gz -O ${INPUTDIR}/02_Trimmed/${FILE_PROCESS}_R2.fastq.gz --thread=1 --html=${INPUTDIR}/out_error/02_Trimmed/${FILE_PROCESS}_file.html --json=${INPUTDIR}/out_error/02_Trimmed/${FILE_PROCESS}_file.json --detect_adapter_for_pe --qualified_quality_phred=20 --length_required=75

############ HiSAT2 Alignment ###############################
hisat2 -p 1 -x /bulk/jpoland/genome/intermedium/C4-5353T1/index/v3/Thinopyrum_intermedium.mainGenome -1 ${INPUTDIR}/02_Trimmed/${FILE_PROCESS}_R1.fastq.gz -2 ${INPUTDIR}/02_Trimmed/${FILE_PROCESS}_R2.fastq.gz -S ${INPUTDIR}/03_Aligned/${FILE_PROCESS}.sam --no-spliced-alignment --no-unal &> ${INPUTDIR}/out_error/03_Aligned/${FILE_PROCESS}.log

########## Filtering Unique Reads ############################
samtools view -h -f 0x2 ${INPUTDIR}/03_Aligned/${FILE_PROCESS}.sam | awk 'substr($1, 0, 1)=="@" || $5 == 60' | samtools view -h -b > ${INPUTDIR}/04_Filtered/${FILE_PROCESS}_filtered.bam

######### index bam file for SNP calling ####################
samtools sort ${INPUTDIR}/04_Filtered/${FILE_PROCESS}_filtered.bam -o ${INPUTDIR}/05_Indexed/${FILE_PROCESS}_sorted.bam -O BAM
samtools index -c ${INPUTDIR}/05_indexed/${FILE_PROCESS}_sorted.bam


echo "Finished JOB on `date`";




#########To RUN###########
#get an array list of all R1 reads in bash before sbatch as 
#% ls $INPUTDIR/01_samples/$R1* | wc -l  #output total number of files

#sbatch --array=0-(output-1) Step_2_5_Pipeline.sh


