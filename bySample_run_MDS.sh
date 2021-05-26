#!/bin/bash -l      

# this script generates a file that contains all of the commands needed to run the pipeline on each sample
# one line per sample
# the wrapSLURM.pl script will then wrap these commands in a sbatch command that can be submitted to SLURM


fgbio=<path to fgbio install>
picard="java -Xmx10g -jar <path to picard jar file>"
seqtk="<path to seqtk install>"
TRIMMOMATIC="<path to trimmomatic install> "
genome=<genome FASTA file>
fastq=<FASTQ dir>
wkdir=<output dir>

# if you shared compute cluster uses the module system 
# we will write the module commands to a file to be included with each 
# sample's job submission 

echo "module load samtools/1.9; \
module load picard/2.18.16; \
module load bwa/0.7.17; \
module load bcftools/1.9; \
module load trimmomatic/0.33; \
module load java/openjdk-8_202" > $wkdir/wrapMods

mkdir -p $wkdir
rm $wkdir/bySampleCommands.txt

cd $fastq
for f in `ls *R1*.fastq.gz`
do
    r1=`ls ${f}`
    r2=${r1/R1_001/R2_001}
    echo $r1
    echo $r2
    name=${r1/_R1_001.fastq.gz/} # specific to the formatting of the fastq file names will need to be edited 
    echo $name
    
    #these are the individual commands needed to run a single sample though the pipeline
    echo "java -version; \
    java -jar $TRIMMOMATIC/trimmomatic.jar PE -phred33 $fastq/$r1 $fastq/$r2 \
    $wkdir/$name.R1.PE.fastq $wkdir/$name.R1.SE.fastq \
    $wkdir/$name.R2.PE.fastq $wkdir/$name.R2.SE.fastq \
    ILLUMINACLIP:/home/umii/public/gopher-pipelines/1.8/resources/all_adapters.fa:2:30:10:2:true \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50;\
    \
    $seqtk sample -s100 $wkdir/$name.R1.PE.fastq 1500000 > $wkdir/$name.R1.PE.1.5M.fastq; \
    $seqtk sample -s100 $wkdir/$name.R2.PE.fastq 1500000 > $wkdir/$name.R2.PE.1.5M.fastq; \
    \
    $fgbio FastqToBam -i $wkdir/$name.R1.PE.1.5M.fastq $wkdir/$name.R2.PE.1.5M.fastq \
    -o $wkdir/$name.fqtobam.bam --sample=$name --library=$name  -r 14M+T +T; \
    \
    $picard SamToFastq I=$wkdir/$name.fqtobam.bam FASTQ=$wkdir/$name.fqtobam.1.fastq \
    SECOND_END_FASTQ=$wkdir/$name.fqtobam.2.fastq; \
    \
    bwa mem -t 2 $genome $wkdir/$name.fqtobam.1.fastq $wkdir/$name.fqtobam.2.fastq | \
    samtools view -b -o $wkdir/$name.mapped.bam; \
    \
    $picard AddOrReplaceReadGroups I=$wkdir/$name.mapped.bam \
    O=$wkdir/$name.mapped_RG.bam RGID=4 RGLB=$name RGPL=illumina RGPU=unit1 RGSM=20; \
    \
    $picard AddOrReplaceReadGroups I=$wkdir/$name.fqtobam.bam \
    O=$wkdir/$name.unmapped_RG.bam RGID=4 RGLB=$name RGPL=illumina RGPU=unit1 RGSM=20;\
    \
    $picard SortSam I=$wkdir/$name.mapped_RG.bam O=$wkdir/$name.sorted_mapped.bam SORT_ORDER=queryname; \
    $picard SortSam I=$wkdir/$name.unmapped_RG.bam O=$wkdir/$name.sorted_unmapped.bam SORT_ORDER=queryname; \
    \
    $picard MergeBamAlignment ALIGNED=$wkdir/$name.sorted_mapped.bam \
    UNMAPPED=$wkdir/$name.sorted_unmapped.bam \
    O=$wkdir/$name.merged_mapped.bam \
    R=$genome; \
    \
    samtools index $wkdir/$name.merged_mapped.bam; \
    \
    $fgbio GroupReadsByUmi -s adjacency -n true -i $wkdir/$name.merged_mapped.bam \
    -o $wkdir/$name.groupedoutput.bam -l 14 -f $wkdir/$name.hist; \
    \
    $fgbio SortBam -i $wkdir/$name.groupedoutput.bam -o $wkdir/$name.sortedoutput.bam -s TemplateCoordinate; \
    \
    $fgbio  CallMolecularConsensusReads -i $wkdir/$name.sortedoutput.bam -o $wkdir/$name.consensus_reads.bam -M 4; \
    \
    $picard SamToFastq I=$wkdir/$name.consensus_reads.bam FASTQ=$wkdir/$name.consensus.1.fastq \
    SECOND_END_FASTQ=$wkdir/$name.consensus.2.fastq; \
    \
    bwa mem -t 2 $genome $wkdir/$name.consensus.1.fastq $wkdir/$name.consensus.2.fastq | \
    samtools view -b -o $wkdir/$name.consensus.mapped.bam; \
    \
    samtools sort  $wkdir/$name.consensus.mapped.bam > $wkdir/$name.consensus.mapped.sort.bam; \
    \
    bcftools mpileup \
    -f $genome \
    -d 100000 \
    -a AD,DP \
    -o $name.vcf \
    $wkdir/$name.consensus.mapped.sort.bam" >> $wkdir/bySampleCommands.txt 
    
done

perl wrapSLURM.pl \
--commands $wkdir/bySampleCommands.txt  \
--ppn 1 \
--mem 150GB \
--walltime 96:00:00 \
--queue amdsmall \
--jobname MDS_HIV1 \
--modules $wkdir/wrapMods > $wkdir/runMDS
