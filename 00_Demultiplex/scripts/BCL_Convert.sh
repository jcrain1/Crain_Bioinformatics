#!/bin/bash 
## SLURM Resource requirement:  
#SBATCH --nodes=1 
#SBATCH --mem=64GB 
#SBATCH --cpus-per-task=16
#SBATCH --job-name=NEXUT
#SBATCH --output=bcl_demux_%A_%a.out  
#SBATCH --error=bcl_demux_%A_%a.err 
#SBATCH --time=00-23:59:00   # Use the form DD-HH:MM:SS
#SBATCH --mail-type=ALL   # same as =BEGIN,FAIL,END

#BCL-convert requires high ulimits. -n is the number of files, -u is the number of processes.
ulimit -n 65535
ulimit -u 32768

INPUTDIR=$1

OUTPUTDIR=$2

SAMPLESHEET=$3

echo "Doing JOB  on |" $(date) " | $SAMPLESHEET | $INPUTDIR | $OUTPUTDIR "
echo ""   

#use bcl-convert
/homes/jcrain/Software/bcl_convert/usr/bin/bcl-convert   --sample-sheet $SAMPLESHEET  --force --bcl-input-directory $INPUTDIR --output-directory $OUTPUTDIR 








echo "Finished JOB on `date`";

