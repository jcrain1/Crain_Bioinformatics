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
outputdir=$3

mkdir -p $outputdir/01_Samples

echo "Doing JOB  on |" $(date) " | $inputfile | $inputdir "
echo "" 

cat $inputfile | while read LINE; do
    
    echo "Working on $LINE"
   cat "$inputdir/00_Demultiplex/$LINE*R1_001.fastq.gz" > $inputdir/01_Samples/$LINE_R1.fastq.gz
   cat "$inputdir/00_Demultiplex/$LINE*R2_001.fastq.gz" > $inputdir/01_Samples/$LINE_R2.fastq.gz


done


echo "Finished JOB on `date`";

