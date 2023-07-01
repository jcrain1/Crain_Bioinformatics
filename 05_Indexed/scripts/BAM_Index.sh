#!/bin/bash -l

#################################################
## Shell script to index BAMs for SNP calling ##
#################################################


#SBATCH --job-name=Index
#SBATCH --mail-type=ALL   # same as =BEGIN,FAIL,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00-23:59:00   # Use the form DD-HH:MM:SS
#SBATCH --mem-per-cpu=2G   # Memory per core, use --mem= for memory per node
#SBATCH --output="index_%A_%a.out"
#SBATCH --error="index_%A_%a.err"

ulimit -u 8192
module load SAMtools/1.10-GCC-8.3.0
module list

INPUTDIR=$1 #get working directory for files
echo ${INPUTDIR}
mkdir -p ${INPUTDIR}/05_Indexed ${INPUTDIR}/out_error/05_Indexed


FILES=($(ls ${INPUTDIR}/04_Filtered/*.bam))


FILENAME=${FILES[$SLURM_ARRAY_TASK_ID]} #this gets file name from outside array
echo ${FILENAME}
FILE_PROCESS=$(basename ${FILENAME} "_filtered.bam") #remove _R1.fq.gz ending

echo "Doing JOB  on |" $(date) " |  $INPUTDIR "

echo $FILE_PROCESS


## index bam file
samtools sort ${INPUTDIR}/04_Filtered/${FILE_PROCESS}_filtered.bam -o ${INPUTDIR}/05_Indexed/${FILE_PROCESS}_sorted.bam -O BAM
samtools index -c ${INPUTDIR}/05_indexed/${FILE_PROCESS}_sorted.bam


echo "Finished JOB on `date`";



#########To RUN###########
#get an array list of all R1 reads in bash before sbatch as 
#export FILES=($(ls *.ba))
# get size of array
#NUMFASTQ=${#FILES[@]}
# now subtract 1 as we have to use zero-based indexing (first cell is 0)
#ARRAYNUMFASTQ=$(($NUMFASTQ - 1))
# now submit to SLURM

#sbatch --array=0-$ARRAYNUMFASTQ BAM_Index.sh




