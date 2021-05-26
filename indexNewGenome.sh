#!/bin/bash -l        
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem=20g
#SBATCH --mail-type=END
#SBATCH --job-name=indexGenome


## use this script to create all of the genome indexs ect needed to run UMI MDS fgbio pipeline

module load samtools/1.9
module load bwa
module load picard/2.18.16

picard="java -Xmx10g -jar /panfs/roc/msisoft/picard/2.18.16/picard.jar"

cd /home/mansklm/shared/genomes/
fasta=HIV2_Edit_mask3urt.fa

bwa index $fasta
$picard CreateSequenceDictionary R=$fasta
samtools faidx $fasta


