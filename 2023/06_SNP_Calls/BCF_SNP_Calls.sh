#!/bin/bash -l

#################################################
## Shell script to call SNPs				 ##
#################################################
####Modified from Ying Hu#####

#SBATCH --job-name=BCF_SNPs
#SBATCH --mail-type=ALL   # same as =BEGIN,FAIL,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=09-00:00:00   # Use the form DD-HH:MM:SS
#SBATCH --mem-per-cpu=16G   # Memory per core, use --mem= for memory per node
#SBATCH --output="bcftools_%A_%a.out"
#SBATCH --error="bcftools_%A_%a.err"
#SBATCH --array=1-21

INPUTDIR=$1 #get working directory for files
region_file=$2

ulimit -u 8192
mkdir -p ${INPUTDIR}/06_VCF ${INPUTDIR}/out_error/06_VCF

echo ${INPUTDIR}

ls ${INPUTDIR}/05_Indexed/*.bam > ${INPUTDIR}/BAMlist.txt


module load BCFtools/1.10.2-foss-2019b
module list

ref="/bulk/jpoland/genome/intermedium/C4-5353T1/assembly/v3_current/Thinopyrum_intermedium_V3_release/Thinopyrum_intermedium/sequences/Thinopyrum_intermedium.mainGenome.fasta"

region=$(head -n $SLURM_ARRAY_TASK_ID $region_file | tail -n 1)


echo "Doing JOB  on |" $(date) " |  $INPUTDIR "

bcftools mpileup --annotate AD,DP,INFO/AD \
                 --skip-indels \
                 -f $ref \
                 --regions $region \
                  -b ${INPUTDIR}/BAMlist.txt -B | bcftools call -m \
                 --variants-only \
                 --skip-variants indels \
                 --output-type v \
                 -o ${INPUTDIR}/06_VCF/${region}.vcf \
                 --group-samples -
                 
echo "Finished JOB on `date`";



#########To RUN###########
#######%sbatch --array=1-21 BCF_SNP_Calls.sh INPUTDIR region_file





