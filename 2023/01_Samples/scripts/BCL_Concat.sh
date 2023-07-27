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

####input variables full path to where the sample names are and output directory will be
inputfile=$1
inputdir=$2


#make output directories and change
mkdir -p $outputdir/01_Samples
cd $inputdir/00_Demultiplex
echo "Current working directory: `pwd`"

echo "Doing JOB  on |" $(date) " | $inputfile | $inputdir "
echo "" 

#reads sample sheet line by line
cat $inputfile | while read LINE; do
    
    echo "Working on $LINE"
    #get all files and sorts them for the R1 file
    cat_files_R1=$(find . -type f -name "${LINE}*R1_001.fastq.gz" | sort)
    echo "$cat_files_R1"
    #pulls the basename of the files
    r1_base=$(basename -a -s "_R1_001.fastq.gz" $cat_files_R1)
    
    #concatenates multiple files together
    echo "$cat_files_R1" | xargs cat > ${inputdir}/01_Samples/${LINE}_R1.fastq.gz
   
   #work on R2 grab files and sort them
   cat_files_R2=$(find . -type f -name "${LINE}*R2_001.fastq.gz" | sort)
 echo "$cat_files_R2"
 r2_base=$(basename -a -s "_R2_001.fastq.gz" $cat_files_R2)
    echo "$cat_files_R2" | xargs cat > ${inputdir}/01_Samples/${LINE}_R2.fastq.gz

#check if the file names in R1 and R2 are in correct order. If this is not equal proble.
if [ "$r1_base" = "$r2_base" ]; then
    echo "Strings are equal."
else
    echo "Strings are not equal."
fi



done


echo "Finished JOB on `date`";

