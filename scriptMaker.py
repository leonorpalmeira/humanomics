
import os

def header(param,step):
    """Outputs the beginning of each script w/ appropriate parameters"""

    for directory in ["SLURMlog","ScriptsDir"]:
        param[directory]=param[directory].rstrip(os.sep)+os.sep

    for directory in ["RawDataDir","ResultsDir"]:
        param[directory]=param[directory].rstrip(os.sep)

    str="""#!/bin/bash\n"""
    str+="""#\n"""
    str+="""#SBATCH --job-name="""+step+"""                                           #@@@ fill with appropriate value\n"""
    str+="""#SBATCH --mail-user="""+param["SLURMemailaddress"]+"""                    #@@@ fill with appropriate value: valid email address\n"""
    str+="""#SBATCH --mail-type="""+param["SLURMemailtype"]+"""                                          #@@@ fill with appropriate value: here, email sent at the END of job\n"""
    str+="""#SBATCH --output="""+param["SLURMlog"]+"""slurm-"""+step+"""-%j.out         #@@@ fill with appropriate value: output log file\n"""
    str+="""#\n"""
    str+="""#SBATCH --ntasks=1                                               #@@@ fill with appropriate value: here 1 task\n"""
    str+="""#SBATCH --cpus-per-task=1                                        #@@@ fill with appropriate value: here 1 core per task\n"""
    str+="""#SBATCH --mem-per-cpu=10000                                      #@@@ fill with appropriate value: here 10Gb of RAM\n"""
    str+="""#SBATCH --time=4:30:00                                           #@@@ fill with appropriate value: here 4h30\n"""
    str+="""#SBATCH --array="""+param["SLURMarray"]+"""                                              #@@@ fill with appropriate value: here samples 1 to 4\n"""
    str+="""\n"""
    str+="""echo "************** SLURM ENV ******************"\n"""
    str+="""echo "TASK_ID:" $TASK_ID\n"""
    str+="""echo "SLURM_ARRAY_TASK_ID:" $SLURM_ARRAY_TASK_ID\n"""
    str+="""echo "SLURM_JOB_NAME:" $SLURM_JOB_NAME\n"""
    str+="""echo "SLURM_NTASKS:" $SLURM_NTASKS\n"""
    str+="""echo "SLURM_CPUS_ON_NODE:" $SLURM_CPUS_ON_NODE\n"""
    str+="""echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST\n"""
    str+="""i=$SLURM_ARRAY_TASK_ID                                             \n"""
    str+="""ncores=$SLURM_CPUS_ON_NODE #@@@ this value should match --cpus-per-task\n"""
    str+="""nthreads=1\n"""
    str+="""mem=1 #@@@ should match in Gb --mem-per-cpu\n"""
    str+="""echo "*******************************************"\n"""
    str+="""echo ""\n"""
    str+="""echo "************** GLOBAL ENV ******************"\n"""
    str+="""export PATH=$GLOBALSCRATCH/bin:$PATH # get R from $BIN (not from the system)\n"""
    str+="""BIN=$GLOBALSCRATCH/bin\n"""
    str+="""SRC=$GLOBALSCRATCH/src\n"""
    str+="""JAVAcustom=$BIN"/java-1.7.0_25 -Xmx"$mem"g -XX:ParallelGCThreads="$ncores" -jar"\n"""
    str+="""echo "BIN:" $BIN\n"""
    str+="""echo "SRC:" $SRC\n"""
    str+="""echo "HOME:" $HOME\n"""
    str+="""echo "PATH:" $PATH\n"""
    str+="""echo "********************************************"\n"""
    str+="""echo ""\n"""
    str+="""echo "************** JOB ENV *********************"\n"""
    str+="""datadir="""+param["RawDataDir"]+""" #@@@ fill this in \n"""
    str+="""resdir="""+param["ResultsDir"]+"""       #@@@ fill this in\n"""
    str+="""fastqcdir=$resdir/fastqc\n"""
    str+="""bamdir=$resdir/bam\n"""
    str+="""samstatdir=$resdir/samstat\n"""
    str+="""picarddir=$resdir/picard\n"""
    str+="""fastqsuffix='"""+param["FastqGzSuffixPE"]+"""' #@@@ fill this in\n"""
    str+="""refname="""+param["ReferenceAssembly"]+"""\n"""
    str+="""ref=$GLOBALSCRATCH/genomes/homo_sapiens/hg19/genome/HG19.fasta\n"""
    str+="""db=$GLOBALSCRATCH/genomes/homo_sapiens/hg19/bwa_hash/HG19.fasta\n"""
    str+="""vcfdir=$GLOBALSCRATCH/genomes/homo_sapiens/hg19/variation\n"""
    str+="""targets="""+param["TargetFile"]+""" #@@@ fill this in\n"""
    str+="""baitNames="""+param["baitNames"]+""" #@@@ fill this in\n"""
    str+="""baitsPicard="""+param["BaitsFilePicard"]+""" #@@@ fill this in\n"""
    str+="""targetsPicard="""+param["TargetFilePicard"]+""" #@@@ fill this in\n"""
    str+="""echo "datadir:" $datadir\n"""
    str+="""echo "fastqcdir:" $fastqcdir\n"""
    str+="""echo "resdir:" $resdir\n"""
    str+="""echo "bamdir:" $bamdir\n"""
    str+="""echo "samstatdir:" $samstatdir\n"""
    str+="""echo "picarddir:" $picarddir\n"""
    str+="""echo "ref": $ref\n"""
    str+="""echo "********************************************"\n"""
    str+="""echo ""\n"""
    str+="""echo "************** This job is the $i th ******************"\n"""
    str+="""let "i -= 1"\n"""
    str+="""echo ""\n"""
    str+="""echo "************** Looking for the sample prefix *******************"\n"""
    str+="""cd $datadir\n"""
    str+="""all=( *$fastqsuffix ) # WARNING: interpreted on the fly\n"""
    str+="""echo ${all[$i]}\n"""
    str+="""prefix=${all[$i]//$fastqsuffix}\n"""
    str+="""\n"""
    str+="""query1=$datadir/$prefix$fastqsuffix\n"""
    str+="""query2=$datadir/$prefix${fastqsuffix/R1/R2}\n"""
    str+="""bam=$bamdir/$prefix.bam\n"""
    str+="""bam_sorted=$bamdir/$prefix"_sorted.bam"\n"""
    str+="""bam_sorted_RG=$bamdir/$prefix"_sorted_RG.bam"\n"""
    str+="""picard_bam=$bamdir/$prefix"_sorted_RG_rmd.bam"\n"""
    str+="""indel_bam=$bamdir/$prefix"_sorted_RG_rmd_indelreal.bam"\n"""
    str+="""bam_ready=$bamdir/$prefix"_sorted_RG_rmd_indelreal_BQSR.bam"\n"""
    str+="""\n"""
    str+="""echo "query1:" $query1\n"""
    str+="""echo "query2:" $query2\n"""
    str+="""\n"""
    str+="""echo "************** Checking for directories ************"\n"""
    str+="""for dir in $resdir $bamdir $samstatdir $picarddir ; do\n"""
    str+="""    if  [ ! -e $dir ] ; then\n"""
    str+="""    mkdir $dir\n"""
    str+="""    fi\n"""
    str+="""done\n"""
    str+="""echo "************** JOB ENV *********************"\n"""
    str+="""\n"""
    str+="""echo "************** Preparing the @RG ReadGroups ************"\n"""
    str+="""date=$(echo $prefix | awk 'BEGIN {FS="_"} {print $1}')\n"""
    str+="""rgdt=${date//R} # This is Iso8601date (2000, but removed in 2004...)\n"""
    str+="""info=$(echo $prefix | awk 'BEGIN {FS="_"} {OFS="_"} {print $1,$2,$3}') # in this experiment, one sample = one library, but several samples per run/lane\n"""
    str+="""rgsm=$(echo $prefix | awk 'BEGIN {FS="_"} {print $2}') # in this experiment, one sample = one library, but several samples per run/lane\n"""
    str+="""rglb=$(echo $prefix | awk 'BEGIN {FS="_"} {print $2}')\n"""
    str+="""lane=$(echo $prefix | awk 'BEGIN {FS="_"} {print $3}')\n"""
    str+="""rgid=$info #"_"$lane"_"$rglb"_"$rgdt\n"""
    str+="""rgpl=illumina\n"""
    str+="""rgpu=hiseq\n"""
    str+="""rgcn=null\n"""
    str+="""rgds="Instrument_HISEQ"\n"""
    str+="""RG="@RG\\\\tID:"$rgid"\\\\tSM:"$rgsm"\\\\tLB:"$rglb"\\\\tPL:"$rgpl"\\\\tPU:"$rgpu"\\\\tCN:"$rgcn"\\\\tDS:"$rgds\n"""
    str+="""echo "******************* Here is my @RG:" $RG "*******************"\n"""
    str+="""\n"""
    return str    

def FastQC(param):
    """Outputs the FastQC.bash script"""

    for directory in ["SLURMlog","ScriptsDir"]:
        param[directory]=param[directory].rstrip(os.sep)+os.sep

    for directory in ["RawDataDir","ResultsDir"]:
        param[directory]=param[directory].rstrip(os.sep)

    out=param["ScriptsDir"]+"0--FastQC.bash"

    str=header(param,"map")
    str+="""echo "************** Launching FastQC report PE1 ******************"\n"""
    str+="""echo "@INPUT" $query1\n"""
    str+="""echo "@OUTPUT"\n"""
    str+="""$BIN/fastqc-0.10.1 --noextract --outdir=$fastqcdir $query1\n"""
    str+="""\n"""
    str+="""echo "************** Launching FastQC report PE2 ******************"\n"""
    str+="""echo "@INPUT" $query2\n"""
    str+="""echo "@OUTPUT"\n"""
    str+="""$BIN/fastqc-0.10.1 --noextract --outdir=$fastqcdir $query2\n"""
    str+="""\n"""
    str+="""echo "************** Finished ******************"\n"""

    try:
        f=open(out,"w")
        f.write(str)
        f.close()
    except IOError, e:
        print "File not found: [", out, "]"
        sys.exit(2)
    return

def MappingAndPreProcessing(param):
    """Outputs the MappingAndPreProcessing.bash script"""

    for directory in ["SLURMlog","ScriptsDir"]:
        param[directory]=param[directory].rstrip(os.sep)+os.sep

    for directory in ["RawDataDir","ResultsDir"]:
        param[directory]=param[directory].rstrip(os.sep)

    out=param["ScriptsDir"]+"1--MappingAndPreProcessing.bash"

    str=header(param,"map")
    str+="""echo "************** Launching BWA-MEM alignments ******************"\n"""
    str+="""echo "@INPUT" $query1\n"""
    str+="""echo "@INPUT" $query2\n"""
    str+="""echo "@OUTPUT" $bam\n"""
    str+="""$BIN/bwa-0.7.7 mem \\\n"""
    str+="""    -t $ncores \\\n"""
    str+="""    -M \\\n"""
    str+="""    -R $RG \\\n"""
    str+="""    $db \\\n"""
    str+="""    $query1 \\\n"""
    str+="""    $query2 | $BIN/samtools-0.1.18 view -bSho $bam -\n"""
    str+="""#  -M: mark shorter split hits as secondary (for Picard/GATK compatibility)\n"""
    str+="""\n"""
    str+="""echo "************** Launching samtools > sorted BAM ******************"\n"""
    str+="""echo "@INPUT" $bam\n"""
    str+="""echo "@OUTPUT" $bam_sorted\n"""
    str+="""$BIN/samtools-0.1.18 sort $bam ${bam_sorted//.bam} # .bam suffix is appended by samtools sort...\n"""
    str+="""rm $bam # the sorted version is enough for (and demanded by) most tools       \n"""
    str+="""\n"""
    str+="""echo "************** Launching samstat > BAM.html ******************"         \n"""
    str+="""echo "@INPUT" $bam_sorted\n"""
    str+="""echo "@OUTPUT" $bam_sorted".html"\n"""
    str+="""$BIN/samstat-1.09 $bam_sorted # can be run before or after sorting\n"""
    str+="""mv $bam_sorted.html $samstatdir/.\n"""
    str+="""\n"""
    str+="""echo "************** Launching Picard MarkDuplicates ******************"\n"""
    str+="""rm $bam_sorted_RG\n"""
    str+="""ln -s $bam_sorted $bam_sorted_RG\n"""
    str+="""echo "@INPUT" $bam_sorted_RG\n"""
    str+="""echo "@OUTPUT" $picard_bam ${picard_bam//.bam}".rmd"\n"""
    str+="""$JAVAcustom $SRC/picard-tools-1.73/MarkDuplicates.jar \\\n"""
    str+="""    INPUT=$bam_sorted_RG  \\\n"""
    str+="""    OUTPUT=$picard_bam \\\n"""
    str+="""    METRICS_FILE=${picard_bam//.bam}".rmd" \\\n"""
    str+="""    REMOVE_DUPLICATES=false \\\n"""
    str+="""    ASSUME_SORTED=true \\\n"""
    str+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    str+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    str+="""\n"""
    str+="""echo "************** Launching samtools > index sorted BAM ******************"\n"""
    str+="""echo "@INPUT" $picard_bam\n"""
    str+="""echo "@OUTPUT" $picard_bam".bai"\n"""
    str+="""$BIN/samtools-0.1.18 index $picard_bam\n"""
    str+="""\n"""
    str+="""echo "************** Launching GATK RealignerTargetCreator ******************"\n"""
    str+="""echo "@INPUT" $picard_bam\n"""
    str+="""echo "@OUTPUT" $picard_bam".intervals"\n"""
    str+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    str+="""    -T RealignerTargetCreator \\\n"""
    str+="""    -nt $nthreads \\\n"""
    str+="""    -R $ref \\\n"""
    str+="""    -I $picard_bam \\\n"""
    str+="""    -o ${picard_bam//.bam}.intervals \\\n"""
    str+="""    --known $vcfdir"/1000G_phase1.indels.hg19_sorted.vcf" \\\n"""
    str+="""    --known $vcfdir"/Mills_and_1000G_gold_standard.indels.hg19_sorted.vcf"\n"""
    str+="""\n"""
    str+="""echo "************** Launching GATK IndelRealigner ******************"\n"""
    str+="""echo "@INPUT" $picard_bam\n"""
    str+="""echo "@INPUT" $picard_bam".intervals"\n"""
    str+="""echo "@OUTPUT" $indel_bam\n"""
    str+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    str+="""    -T IndelRealigner \\\n"""
    str+="""    -R $ref \\\n"""
    str+="""    -I $picard_bam \\\n"""
    str+="""    -targetIntervals ${picard_bam//.bam}.intervals \\\n"""
    str+="""    -o $indel_bam \\\n"""
    str+="""    --knownAlleles $vcfdir"/1000G_phase1.indels.hg19_sorted.vcf" \\\n"""
    str+="""    --knownAlleles $vcfdir"/Mills_and_1000G_gold_standard.indels.hg19_sorted.vcf"\n"""
    str+="""\n"""
    str+="""echo "************** Launching samtools > index sorted BAM ******************"\n"""
    str+="""echo "@INPUT" $indel_bam".bai"\n"""
    str+="""echo "@OUTPUT" $indel_bam\n"""
    str+="""$BIN/samtools-0.1.18 index $indel_bam\n"""
    str+="""\n"""
    str+="""echo "************** Launching GATK BaseRecalibrator 1st pass ******************"\n"""
    str+="""echo "@INPUT" $indel_bam\n"""
    str+="""echo "@OUTPUT" ${indel_bam//.bam}_recal1.grp\n"""
    str+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    str+="""    -T BaseRecalibrator \\\n"""
    str+="""    -nct $ncores \\\n"""
    str+="""    -I $indel_bam \\\n"""
    str+="""    -R $ref \\\n"""
    str+="""    -L $targets \\\n"""
    str+="""    -knownSites $vcfdir"/dbsnp_137.hg19_sorted.vcf" \\\n"""
    str+="""    -knownSites $vcfdir"/Mills_and_1000G_gold_standard.indels.hg19_sorted.vcf" \\\n"""
    str+="""    -knownSites $vcfdir"/1000G_phase1.indels.hg19_sorted.vcf" \\\n"""
    str+="""    -o ${indel_bam//.bam}_recal1.grp\n"""
    str+="""\n"""
    str+="""echo "************** Launching GATK BaseRecalibrator 2nd pass ******************"\n"""
    str+="""echo "@INPUT" $indel_bam\n"""
    str+="""echo "@INPUT" ${indel_bam//.bam}_recal1.grp\n"""
    str+="""echo "@OUTPUT" ${indel_bam//.bam}_recal2.grp\n"""
    str+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    str+="""    -T BaseRecalibrator \\\n"""
    str+="""    -BQSR ${indel_bam//.bam}_recal1.grp \\\n"""
    str+="""    -nct $ncores \\\n"""
    str+="""    -I $indel_bam \\\n"""
    str+="""    -R $ref \\\n"""
    str+="""    -L $targets \\\n"""
    str+="""    -knownSites $vcfdir"/dbsnp_137.hg19_sorted.vcf" \\\n"""
    str+="""    -knownSites $vcfdir"/Mills_and_1000G_gold_standard.indels.hg19_sorted.vcf" \\\n"""
    str+="""    -knownSites $vcfdir"/1000G_phase1.indels.hg19_sorted.vcf" \\\n"""
    str+="""    -o ${indel_bam//.bam}_recal2.grp\n"""
    str+="""\n"""
    str+="""echo "************** Launching GATK AnalyzeCovariates > Plots ******************"\n"""
    str+="""echo "@INPUT" ${indel_bam//.bam}_recal1.grp\n"""
    str+="""echo "@INPUT" ${indel_bam//.bam}_recal2.grp\n"""
    str+="""echo "@OUTPUT" ${indel_bam//.bam}_BQSR.csv\n"""
    str+="""echo "@OUTPUT" ${indel_bam//.bam}_BQSR.pdf\n"""
    str+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    str+="""    -T AnalyzeCovariates \\\n"""
    str+="""    -R $ref \\\n"""
    str+="""    -L $targets \\\n"""
    str+="""    -before ${indel_bam//.bam}_recal1.grp \\\n"""
    str+="""    -after ${indel_bam//.bam}_recal2.grp \\\n"""
    str+="""    -csv ${indel_bam//.bam}_BQSR.csv \\\n"""
    str+="""    -plots ${indel_bam//.bam}_BQSR.pdf\n"""
    str+="""\n"""
    str+="""echo "************** Launching GATK PrintReads ******************"\n"""
    str+="""echo "@INPUT" $indel_bam\n"""
    str+="""echo "@INPUT" ${indel_bam//.bam}_recal2.grp\n"""
    str+="""echo "@OUTPUT" $bam_ready\n"""
    str+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    str+="""    -T PrintReads \\\n"""
    str+="""    -nct $ncores \\\n"""
    str+="""    -R $ref \\\n"""
    str+="""    -L $targets \\\n"""
    str+="""    -I $indel_bam \\\n"""
    str+="""    -BQSR ${indel_bam//.bam}_recal2.grp \\\n"""
    str+="""    -o $bam_ready\n"""
    str+="""\n"""
    str+="""echo "************** Launching samtools > index sorted BAM ******************"\n"""
    str+="""echo "@INPUT" $bam_ready\n"""
    str+="""echo "@OUTPUT" $bam_ready.bai\n"""
    str+="""$BIN/samtools-0.1.18 index $bam_ready # .bai indexes positions on the reference genome and some tools can use this directly\n"""
    str+="""\n"""
    str+="""echo "************** Launching samstat > BAM.html ******************"\n"""
    str+="""echo "@INPUT" $bam_ready\n"""
    str+="""echo "@OUTPUT" $bam_ready.html\n"""
    str+="""$BIN/samstat-1.09 $bam_ready\n"""
    str+="""mv $bam_ready.html $samstatdir/.\n"""
    str+="""\n"""
    str+="""echo "************** Launching samtools idxstats ******************"\n"""
    str+="""echo "@INPUT" $bam_ready\n"""
    str+="""echo "@OUTPUT" ${bam_ready//.bam}.idxstats\n"""
    str+="""$BIN/samtools-0.1.18 idxstats $bam_ready > ${bam_ready//.bam}.idxstats\n"""
    str+="""\n"""
    str+="""echo "************** Finished ******************"\n"""

    try:
        f=open(out,"w")
        f.write(str)
        f.close()
    except IOError, e:
        print "File not found: [", out, "]"
        sys.exit(2)
    return

def QualityControl(param):
    """Outputs the QualityControl.bash script"""

    for directory in ["SLURMlog","ScriptsDir"]:
        param[directory]=param[directory].rstrip(os.sep)+os.sep

    for directory in ["RawDataDir","ResultsDir"]:
        param[directory]=param[directory].rstrip(os.sep)

    out=param["ScriptsDir"]+"2a--QualityControl.bash"

    str=header(param,"qc")
    str+="""echo 'BAM for QC:' $bam_ready\n"""
    str+="""\n"""
    str+="""echo '************** Checking for directories ************'\n"""
    str+="""for dir in $resdir $bamdir $samstatdir $picarddir ; do\n"""
    str+="""    if  [ ! -e $dir ] ; then\n"""
    str+="""	mkdir $dir\n"""
    str+="""    fi\n"""
    str+="""done\n"""
    str+="""echo '************** JOB ENV *********************'\n"""
    str+="""\n"""
    str+="""echo '************** Launching Samtools Flagstat ******************'\n"""
    str+="""$BIN/samtools-0.1.18 flagstat $bam_ready >  ${bam_ready//.bam}.flagstat\n"""
    str+="""\n"""
    str+="""echo '************** Launching Picard EstimateLibraryComplexity *****************'\n"""
    str+="""$JAVAcustom $SRC/picard-tools-1.73/EstimateLibraryComplexity.jar \\\n"""
    str+="""    INPUT=$bam_ready \\\n"""
    str+="""    OUTPUT=${bam_ready//.bam}'.lc' \\\n"""
    str+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    str+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    str+="""#ASSUME_SORTED=true         \n"""
    str+="""\n"""
    str+="""echo '************** Launching Picard MeanQualityByCycle *****************'\n"""
    str+="""$JAVAcustom $SRC/picard-tools-1.73/MeanQualityByCycle.jar \\\n"""
    str+="""    INPUT=$bam_ready \\\n"""
    str+="""    OUTPUT=$bam_ready.qbc \\\n"""
    str+="""    CHART_OUTPUT=${bam_ready//.bam}'_qbc.pdf' \\\n"""
    str+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    str+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    str+="""\n"""
    str+="""echo '************** Launching Picard QualityScoreDistribution ****************'\n"""
    str+="""$JAVAcustom $SRC/picard-tools-1.73/QualityScoreDistribution.jar \\\n"""
    str+="""    INPUT=$bam_ready \\\n"""
    str+="""    OUTPUT=${bam_ready//.bam}'.qsd' \\\n"""
    str+="""    CHART_OUTPUT=${bam_ready//.bam}'_qsd.pdf' \\\n"""
    str+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    str+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    str+="""\n"""
    str+="""echo '************** Launching Picard CollectGcBiasMetrics *****************'\n"""
    str+="""$JAVAcustom $SRC/picard-tools-1.73/CollectGcBiasMetrics.jar \\\n"""
    str+="""    INPUT=$bam_ready \\\n"""
    str+="""    REFERENCE_SEQUENCE=$ref \\\n"""
    str+="""    OUTPUT=${bam_ready//.bam}'.gcb' \\\n"""
    str+="""    CHART_OUTPUT=${bam_ready//.bam}'_gcb.pdf' \\\n"""
    str+="""    SUMMARY_OUTPUT=${bam_ready//.bam}'_gcb.sum' \\\n"""
    str+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    str+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    str+="""\n"""
    str+="""echo '************** Launching Picard CollectInsertSizeMetrics *****************'\n"""
    str+="""$JAVAcustom $SRC/picard-tools-1.73/CollectInsertSizeMetrics.jar \\\n"""
    str+="""    INPUT=$bam_ready \\\n"""
    str+="""    OUTPUT=${bam_ready//.bam}'.ism' \\\n"""
    str+="""    HISTOGRAM_FILE=${bam_ready//.bam}'_ism.pdf' \\\n"""
    str+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    str+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    str+="""\n"""
    str+="""echo '************** Launching Picard CollectAlignmentSummaryMetrics ****************'\n"""
    str+="""$JAVAcustom $SRC/picard-tools-1.73/CollectAlignmentSummaryMetrics.jar \\\n"""
    str+="""    INPUT=$bam_ready \\\n"""
    str+="""    OUTPUT=${bam_ready//.bam}'.asm' \\\n"""
    str+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    str+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    str+="""\n"""
    str+="""echo '************** Launching IGVtools count **********'\n"""
    str+="""$JAVAcustom $SRC/IGVTools-1.5.10/igvtools.jar count \\\n"""
    str+="""    $bam_ready  \\\n"""
    str+="""    ${bam_ready//.bam}.tdf \\\n"""
    str+="""    $refname\n"""
    str+="""\n"""
    str+="""# WARNING: careful with the manifest and target files here!\n"""
    str+="""echo '************** Launching Picard CalculateHSMetrics ******************'\n"""
    str+="""$JAVAcustom $SRC/picard-tools-1.73/CalculateHsMetrics.jar \\\n"""
    str+="""    BAIT_SET_NAME=$baitNames \\\n"""
    str+="""    BAIT_INTERVALS=$baitsPicard \\\n"""
    str+="""    TARGET_INTERVALS=$targetsPicard \\\n"""
    str+="""    INPUT=$bam_ready \\\n"""
    str+="""    OUTPUT=${bam_ready//.bam}'.hsm' \\\n"""
    str+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    str+="""    REFERENCE_SEQUENCE=$ref \\\n"""
    str+="""    METRIC_ACCUMULATION_LEVEL=SAMPLE \\\n"""
    str+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    str+="""\n"""
    str+="""echo '************** Launching GATK DiagnoseTargets ******************'\n"""
    str+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    str+="""    -T DiagnoseTargets \\\n"""
    str+="""    -R $ref \\\n"""
    str+="""    -I $bam_ready \\\n"""
    str+="""    -o ${bam_ready//.bam}.diagt \\\n"""
    str+="""    -L $targets\n"""
    str+="""\n"""
    str+="""echo '************** Launching GATK DepthOfCoverage ******************'\n"""
    str+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    str+="""    -T DepthOfCoverage \\\n"""
    str+="""    -R $ref \\\n"""
    str+="""    -I $bam_ready \\\n"""
    str+="""    -o ${bam_ready//.bam}'.cov' \\\n"""
    str+="""    -pt sample \\\n"""
    str+="""    --omitDepthOutputAtEachBase \\\n"""
    str+="""    -ct 20 \\\n"""
    str+="""    -L $targets\n"""
    str+="""\n"""
    str+="""echo '************** Finished ******************'\n"""

    try:
        f=open(out,"w")
        f.write(str)
        f.close()
    except IOError, e:
        print "File not found: [", out, "]"
        sys.exit(2)
    return

def HaplotypeCaller(param):
    """Outputs the HaplotypeCaller.bash script"""

    for directory in ["SLURMlog","ScriptsDir"]:
        param[directory]=param[directory].rstrip(os.sep)+os.sep

    for directory in ["RawDataDir","ResultsDir"]:
        param[directory]=param[directory].rstrip(os.sep)

    out=param["ScriptsDir"]+"2b--HaplotypeCaller.bash"

    str=header(param,"hc")
    str+="""middfix='_HC3.1' #@@@ fill this in\n"""
    str+="""\n"""
    str+="""echo 'BAM for HC:' $bam_ready\n"""
    str+="""\n"""
    str+="""echo '************** Checking for directories ************'\n"""
    str+="""for dir in $resdir $bamdir $samstatdir $picarddir ; do\n"""
    str+="""    if  [ ! -e $dir ] ; then\n"""
    str+="""	mkdir $dir\n"""
    str+="""    fi\n"""
    str+="""done\n"""
    str+="""echo '************** JOB ENV *********************'\n"""
    str+="""\n"""
    str+="""echo '************** Launching GATK HaplotypeCaller ******************'\n"""
    str+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    str+="""    -T HaplotypeCaller \\\n"""
    str+="""    -R $ref \\\n"""
    str+="""    -L $targets \\\n"""
    str+="""    -I $bamdir/$prefix'_sorted_RG_rmd_indelreal_BQSR.bam' \\\n"""
    str+="""    -o $bamdir/$prefix'_sorted_RG_rmd_indelreal_BQSR'$middfix'.gvcf' \\\n"""
    str+="""    --dbsnp $vcfdir'/dbsnp_137.hg19_sorted.vcf' \\\n"""
    str+="""    --emitRefConfidence GVCF \\\n"""
    str+="""    --variant_index_type LINEAR \\\n"""
    str+="""    --variant_index_parameter 128000 \\\n"""
    str+="""    --pair_hmm_implementation VECTOR_LOGLESS_CACHING \n"""
    str+="""\n"""
    str+="""echo '************** Finished ******************'\n"""

    try:
        f=open(out,"w")
        f.write(str)
        f.close()
    except IOError, e:
        print "File not found: [", out, "]"
        sys.exit(2)
    return
