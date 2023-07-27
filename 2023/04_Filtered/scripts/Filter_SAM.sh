#!/bin/bash -l

#################################################
## Shell script to filter SAMs for unique align ##
#################################################

#SBATCH --job-name=SAMFilter
#SBATCH --time=00-23:00:00 #Use the from DD-HH:MM:SS
#SBATCH --mail-type=ALL   # same as =BEGIN,FAIL,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1G   # Memory per core, use --mem= for memory per node
#SBATCH --output="samtools_%A_%a.out"
#SBATCH --error="samtools_%A_%a.err"


ulimit -u 8192
module load SAMtools/1.10-GCC-8.3.0
module list


INPUTDIR=$1 #get working directory for files
echo ${INPUTDIR}
mkdir -p ${INPUTDIR}/04_Filtered ${INPUTDIR}/out_error/04_Filtered 

FILES=($(ls ${INPUTDIR}/03_Aligned/*.sam))

FILENAME=${FILES[$SLURM_ARRAY_TASK_ID]} #this gets file name from outside array
echo ${FILENAME}
FILE_PROCESS=$(basename ${FILENAME} ".sam") #remove _R1.fq.gz ending

echo "Doing JOB  on |" $(date) " |  $INPUTDIR "

echo $FILE_PROCESS

DIR=`pwd` #get working directory


samtools view -h -f 0x2 ${INPUTDIR}/03_Aligned/${FILE_PROCESS}.sam | awk 'substr($1, 0, 1)=="@" || $5 == 60' | samtools view -h -b > ${INPUTDIR}/04_Filtered/${FILE_PROCESS}_filtered.bam

echo "Finished JOB on `date`";


#########To RUN###########
#get an array list of all sam files in bash before sbatch as 
#% ls $DIR/03_Aligned/*.sam | wc -l  #output total number of files

#sbatch --array=0-(output number -1) Filter_SAM.sh
