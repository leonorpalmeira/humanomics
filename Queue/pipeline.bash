#!/usr/bin/env bash
# -*- coding: utf-8 -*-                                                                     
#
#SBATCH --job-name=Humanomics
#SBATCH --time=48:00:00
#SBATCH --mail-user=leonor.palmeira@gmail.com 
#SBATCH --mail-type=END
#SBATCH --output=/scratch/ulg/genan/palmeira/tmp/slurm-Humanomics-%j.out

echo "============== Setting ENV =============="
mem=4
ncores=4
BIN=$GLOBALSCRATCH/bin
JAVAcustom=$BIN"/java-1.7.0_25 -Xmx"$mem"g -XX:ParallelGCThreads="$ncores" -jar"
GATK=$GLOBALSCRATCH"/src/Queue-3.2-2"
Rscript=$HOME"/bin/Rscript-3.1.1"
#input=$GLOBALSCRATCH"/work/rodrigo/R120807_MG_L004_sorted_RG_rmd_indelreal_BQSR.bam"
#scalascript=$HOME"/scripts/work/rodrigo/gatk.qpp.scala"
echo "============== Started =================="
echo `date`
cd $GLOBALSCRATCH"/work/rodrigo/map/"
# NOTE to self:
#  Unsupported option: --cpus-per-task

# CAREFUL : in the bash script, I should check that the following dirs exist:
# fastqc
# samstat
# bam
# .qlog # pour le output de queue
scalascript=$HOME"/scripts/humanomics/scala/HumanomicsPipeline.scala"
jobpars=""
inputdir=$GLOBALSCRATCH"/work/rodrigo/map/"
#    --input_file_1 /scratch/ulg/genan/palmeira/work/rodrigo/map/ISDBM376558_1_small.fastq.gz \
#    --input_file_2 /scratch/ulg/genan/palmeira/work/rodrigo/map/ISDBM376558_2_small.fastq.gz \
outputdir=$GLOBALSCRATCH"/work/rodrigo/map/"
ref=$GLOBALSCRATCH"/genomes/homo_sapiens/hg19/genome/HG19.fasta"
refbwa=$GLOBALSCRATCH"/genomes/homo_sapiens/hg19/bwa_hash/HG19.fasta"
noET=$HOME"/.ssh/leonor.palmeira_ulg.ac.be.key"
$JAVAcustom ${GATK}"/Queue.jar" \
    -S $scalascript \
    -I1 $GLOBALSCRATCH"/work/rodrigo/map/ISDBM376558_1_small.fastq.gz" \
    -I2 $GLOBALSCRATCH"/work/rodrigo/map/ISDBM376558_2_small.fastq.gz" \
    -RG "@RG\\tID:Sample\\tSM:Sample\\tLB:Sample\\tPL:illumina\\tPU:hiseq\\tCN:null\\tDS:Instrument_HISEQ" \
    -snps "/scratch/ulg/genan/palmeira/genomes/homo_sapiens/hg19/variation/dbsnp_137.hg19_sorted.vcf" \
    -indels "/scratch/ulg/genan/palmeira/genomes/homo_sapiens/hg19/variation/Mills_and_1000G_gold_standard.indels.hg19_sorted.vcf" \
    -indels "/scratch/ulg/genan/palmeira/genomes/homo_sapiens/hg19/variation/1000G_phase1.indels.hg19_sorted.vcf" \
    --output_dir $outputdir \
    -l DEBUG \
    -R $ref \
    -Rbwa $refbwa \
    -jobRunner Drmaa \
    -jobNative "--time=12:00:00 --nodes=1 --ntasks-per-node="$ncores" --mem-per-cpu=4000" \
    -keepIntermediates \
    -run \
    -et NO_ET \
    -K $noET \
    -P 4 \
    -N $ncores \
    -C $ncores \
    -javacustom "/scratch/ulg/genan/palmeira/bin/java-1.7.0_25 -Xmx"$mem"g -XX:ParallelGCThreads="$ncores" -jar" \
    -bwa "/scratch/ulg/genan/palmeira/bin/bwa-0.7.7" \
    -samtools "/scratch/ulg/genan/palmeira/bin/samtools-0.1.18" \
    -samstat "/scratch/ulg/genan/palmeira/bin/samstat-1.09" \
    -fastqc "/scratch/ulg/genan/palmeira/bin/fastqc-0.10.1" \
    -picard "/scratch/ulg/genan/palmeira/src/picard-tools-1.73" \
    -startFromScratch
#    --input_dir $inputdir \

echo `date`
echo "============== Finished ================="


# # other ENV
# ref=$GLOBALSCRATCH"/genomes/homo_sapiens/hg19/genome/HG19.fasta"
# vcfdir=$GLOBALSCRATCH"/genomes/homo_sapiens/hg19/variation"
# INDELVCF=$vcfdir"/Mills_and_1000G_gold_standard.indels.hg19_sorted_chr1-9XYnoMT.vcf"
# DBSNP=$vcfdir"/dbsnp_137.hg19_sorted_chr1-9XYnoMT.vcf"

# # Input file at the comand line
# #input=$1

# java -Xmx500M -jar ${GATK}/Queue.jar \
#     -S $scalascript \
#     -P 5 \
#     -R $ref \
#     -I $input \
#     -known $INDELVCF \
#     -N 6 -C 6 \
#     -jobRunner Drmaa \
#     -jobNative "--time=12:00:00 --nodes=1 --ntasks-per-node=6 --mem-per-cpu=4000" \
#     -run \
#     -startFromScratch

# echo `date`
# echo "============== Finished =================="


