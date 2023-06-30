#!/bin/bash -l

#################################################
## Shell script to run HiSat2 Alginment       ##
#################################################

#SBATCH --job-name=hisat2
#SBATCH --time=00-23:59:00 #Use the from DD-HH:MM:SS
#SBATCH --mail-type=ALL   # same as =BEGIN,FAIL,END
#SBATCH --mem-per-cpu=36G #Memory per core, use --mem= form memory per node
#SBATCH --output="Hisat2_%A_%a.out"
#SBATCH --error="Hisat2_%A_%a.err"


ulimit -u 8192
mkdir -p ${INPUTDIR}/03_Aligned ${INPUTDIR}/out_error/03_Aligned 
INPUTDIR=$1 #get working directory for files

FILES=($(ls ${INPUTDIR}/02_Trimmed/*R1.fastq.gz))
module load HISAT2
module list

FILENAME=${FILES[$SLURM_ARRAY_TASK_ID]} #this gets file name from outside array
echo ${FILENAME}
FILE_PROCESS=$(basename ${FILENAME} "_R1.fastq.gz") #remove _R1.fq.gz ending

echo "Current working directory: `pwd`"

echo "Doing JOB  on |" $(date) " |  $INPUTDIR "
echo "" 
echo $FILE_PROCESS


hisat2 -p 1 -x /bulk/jpoland/genome/intermedium/C4-5353T1/index/v3/Thinopyrum_intermedium.mainGenome -1 ${INPUTDIR}/02_Trimmed/${FILE_PROCESS}_R1.fastq.gz -2 ${INPUTDIR}/02_Trimmed/${FILE_PROCESS}_R2.fastq.gz -S ${INPUTDIR}/03_Aligned/${FILE_PROCESS}.sam --no-spliced-alignment --no-unal &> ${INPUTDIR}/out_error/03_Trimmed/${FILE_PROCESS}.log


#########To RUN###########
#get an array list of all R1 reads in bash before sbatch as 
#% ls $DIR/02_Trimmed/$R1* | wc -l  #output total number of files

#sbatch --array=0-(output number -1) HiSat2_Alignment.sh


echo "Finished JOB on `date`";




