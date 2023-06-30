#!/bin/bash -l

#################################################
## Shell script to run FastP trimming         ##
#################################################

####Array script for fastp trimming#####
####Modified from Laxman Adhikari and Sandesh Shrestha#####
##https://docs.rc.fas.harvard.edu/kb/submitting-large-numbers-of-jobs/####

#SBATCH --job-name=trimfastp
#SBATCH --time=00-23:59:00 #Use the from DD-HH:MM:SS
#SBATCH --mem-per-cpu=2G #Memory per core, use --mem= form memory per node
#SBATCH --mail-type=ALL   # same as =BEGIN,FAIL,END 
#SBATCH --output="FastP_%A_%a.out"
#SBATCH --error="FastP_%A_%a.err"

ulimit -u 8192
INPUTDIR=$1 #get working directory for files

mkdir -p ${INPUTDIR}/02_Trimmed ${INPUTDIR}/out_error/02_Trimmed
FILES=($(ls ${INPUTDIR}/01_Samples/*R1.fastq.gz))

FILENAME=${FILES[$SLURM_ARRAY_TASK_ID]} #this gets file name from outside array
echo ${FILENAME}
FILE_PROCESS=$(basename ${FILENAME} "_R1.fastq.gz") #remove _R1.fq.gz ending

echo "Current working directory: `pwd`"

echo "Doing JOB  on |" $(date) " |  $INPUTDIR "
echo "" 
echo $FILE_PROCESS


/homes/jcrain/Software/fastp/fastp -i ${INPUTDIR}/01_Samples/${FILE_PROCESS}_R1.fastq.gz -I ${INPUTDIR}/01_Samples/${FILE_PROCESS}_R2.fastq.gz -o ${INPUTDIR}/02_Trimmed/${FILE_PROCESS}_R1.fastq.gz -O ${INPUTDIR}/02_Trimmed/${FILE_PROCESS}_R2.fastq.gz --thread=1 --html=${INPUTDIR}/out_error/02_Trimmed/${FILE_PROCESS}_file.html --json=${INPUTDIR}/out_error/02_Trimmed/${FILE_PROCESS}_file.json --detect_adapter_for_pe --qualified_quality_phred=20 --length_required=150


echo "Finished JOB on `date`";


#########To RUN###########
#get an array list of all R1 reads in bash before sbatch as 
#% ls $INPUTDIR/01_samples/$R1* | wc -l  #output total number of files

#sbatch --array=0-(output-1) FASTQ_Trimming.sh


