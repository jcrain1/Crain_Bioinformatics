#!/bin/bash -l

#################################################
## Shell script to concatenate fastq files     ##
#################################################

#SBATCH --job-name=Concat
#SBATCH --mail-type=ALL   # same as =BEGIN,FAIL,END --mail-user=jcrain@ksu.edu
#SBATCH --time=0-23:00:00   # Use the form DD-HH:MM:SS
#SBATCH --mem-per-cpu=4G   # Memory per core, use --mem= for memory per node
#SBATCH --partition=ksu-plantpath-jpoland.q,batch.q,killable.q

ulimit -u 8192

inputfile=$1
inputdir=$2

mkdir -p $outputdir/01_Samples
cd $inputdir/00_Demultiplex
echo "Current working directory: `pwd`"

echo "Doing JOB  on |" $(date) " | $inputfile | $inputdir "
echo "" 

cat $inputfile | while read LINE; do
    
    echo "Working on $LINE"
    cat_files_R1=$(find . -type f -name "${LINE}*R1_001.fastq.gz")
    echo "$cat_files_R1"
    echo "$cat_files" | xargs cat > ${inputdir}/01_Samples/${LINE}_R1.fastq.gz
   #cat ${cat_files_R1} > ${inputdir}/01_Samples/${LINE}_R1.fastq.gz
   #cat_files_R2=$(find . -type f -name "${LINE}*R2_001.fastq.gz")
  # echo "$cat_files_R2"
  # cat ${cat_files_R2} > ${inputdir}/01_Samples/${LINE}_R2.fastq.gz


done


echo "Finished JOB on `date`";

