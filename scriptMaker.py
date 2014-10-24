
import os

def noPhoneHome(param,sout):
    """Adds the noPhoneHome option for all GATK commands when a noET
    key is available in the configuration file"""

    stmp=sout.rstrip()+""" \\\n""" # replace w/ appropriate end of line
    stmp+="""    -et NO_ET \\\n"""
    stmp+="""    -K $noET\n"""
    stmp+="""\n"""
    return stmp

def header(param,step):
    """Outputs the beginning of each script w/ appropriate parameters"""

    for directory in ["SLURMlog","ScriptsDir"]:
        param[directory]=param[directory].rstrip(os.sep)+os.sep

    for directory in ["RawDataDir","ResultsDir"]:
        param[directory]=param[directory].rstrip(os.sep)

    sout="""#!/usr/bin/env bash\n"""
    sout+="""# -*- coding: utf-8 -*-\n"""
    sout+="""#\n"""
    sout+="""#SBATCH --job-name="""+step+"""\n"""
    sout+="""#SBATCH --mail-user="""+param["SLURMemailaddress"]+"""\n"""
    sout+="""#SBATCH --mail-type="""+param["SLURMemailtype"]+"""\n"""
    sout+="""#SBATCH --output="""+param["SLURMlog"]+"""slurm-"""+step+"""-%j.out\n"""
    sout+="""#\n"""
    sout+="""#SBATCH --ntasks=1                                               #@@@ fill with appropriate value: here 1 task\n"""
    sout+="""#SBATCH --cpus-per-task=1                                        #@@@ fill with appropriate value: here 1 core per task\n"""
    sout+="""#SBATCH --mem-per-cpu=10000                                      #@@@ fill with appropriate value: here 10Gb of RAM\n"""
    sout+="""#SBATCH --time=4:30:00                                           #@@@ fill with appropriate value: here 4h30\n"""
    sout+="""#SBATCH --array="""+param["SLURMarray"]+"""\n"""
    sout+="""\n"""
    sout+="""echo "************** SLURM ENV ******************"\n"""
    sout+="""echo "TASK_ID:" $TASK_ID\n"""
    sout+="""echo "SLURM_ARRAY_TASK_ID:" $SLURM_ARRAY_TASK_ID\n"""
    sout+="""echo "SLURM_JOB_NAME:" $SLURM_JOB_NAME\n"""
    sout+="""echo "SLURM_NTASKS:" $SLURM_NTASKS\n"""
    sout+="""echo "SLURM_CPUS_ON_NODE:" $SLURM_CPUS_ON_NODE\n"""
    sout+="""echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST\n"""
    sout+="""i=$SLURM_ARRAY_TASK_ID                                             \n"""
    sout+="""ncores=$SLURM_CPUS_ON_NODE #@@@ this value should match --cpus-per-task\n"""
    sout+="""nthreads=1\n"""
    sout+="""mem=10 #@@@ should match in Gb --mem-per-cpu\n"""
    sout+="""echo "*******************************************"\n"""
    sout+="""echo ""\n"""
    sout+="""echo "************** GLOBAL ENV ******************"\n"""
    sout+="""export PATH=$GLOBALSCRATCH/bin:$PATH # get R from $BIN (not from the system)\n"""
    sout+="""BIN=$GLOBALSCRATCH/bin\n"""
    sout+="""SRC=$GLOBALSCRATCH/src\n"""
    sout+="""JAVAcustom=$BIN"/java-1.7.0_25 -Xmx"$mem"g -XX:ParallelGCThreads="$ncores" -jar"\n"""
    if param.has_key("noPhoneHome"):
        sout+="""noET="""+param["noPhoneHome"]+"""\n"""
    sout+="""echo "BIN:" $BIN\n"""
    sout+="""echo "SRC:" $SRC\n"""
    sout+="""echo "HOME:" $HOME\n"""
    sout+="""echo "PATH:" $PATH\n"""
    if param.has_key("noPhoneHome"):
        sout+="""echo "No Phone Home (GATK):" $noET\n"""
    sout+="""echo "********************************************"\n"""
    sout+="""echo ""\n"""
    sout+="""echo "************** JOB ENV *********************"\n"""
    sout+="""datadir="""+param["RawDataDir"]+""" #@@@ fill this in \n"""
    sout+="""resdir="""+param["ResultsDir"]+"""       #@@@ fill this in\n"""
    sout+="""fastqcdir=$resdir/fastqc\n"""
    sout+="""bamdir=$resdir/bam\n"""
    sout+="""resvcfdir=$resdir/vcf\n"""
    sout+="""samstatdir=$resdir/samstat\n"""
    sout+="""picarddir=$resdir/picard\n"""
    sout+="""fastqsuffix='"""+param["FastqGzSuffixPE"]+"""' #@@@ fill this in\n"""
    sout+="""refname="""+param["ReferenceAssembly"]+"""\n"""
    sout+="""ref=$GLOBALSCRATCH/genomes/homo_sapiens/hg19/genome/HG19.fasta\n"""
    sout+="""dict=$GLOBALSCRATCH/genomes/homo_sapiens/hg19/genome/HG19.dict\n"""
    sout+="""db=$GLOBALSCRATCH/genomes/homo_sapiens/hg19/bwa_hash/HG19.fasta\n"""
    sout+="""vcfdir=$GLOBALSCRATCH/genomes/homo_sapiens/hg19/variation\n"""
    if param["AnalysisMode"]=="EXOME":
        sout+="""echo "*************** EXOME Analysis Mode *****************"\n"""
        sout+="""targets="""+param["TargetFile"]+""" #@@@ fill this in\n"""
        sout+="""baitNames="""+param["baitNames"]+""" #@@@ fill this in\n"""
        sout+="""baitsPicard="""+param["BaitsFilePicard"]+""" #@@@ fill this in\n"""
        sout+="""targetsPicard="""+param["TargetFilePicard"]+""" #@@@ fill this in\n"""
    else:
        sout+="""echo "*************** GENOME Analysis Mode *****************"\n"""
    sout+="""echo "datadir:" $datadir\n"""
    sout+="""echo "fastqcdir:" $fastqcdir\n"""
    sout+="""echo "resdir:" $resdir\n"""
    sout+="""echo "bamdir:" $bamdir\n"""
    sout+="""echo "resvcfdir:" $resvcfdir\n"""
    sout+="""echo "samstatdir:" $samstatdir\n"""
    sout+="""echo "picarddir:" $picarddir\n"""
    sout+="""echo "ref": $ref\n"""
    sout+="""echo "********************************************"\n"""
    sout+="""echo ""\n"""
    sout+="""echo "************** This job is the $i th ******************"\n"""
    sout+="""let "i -= 1"\n"""
    sout+="""echo ""\n"""
    sout+="""echo "************** Looking for the sample prefix *******************"\n"""
    sout+="""cd $datadir\n"""
    sout+="""all=( *$fastqsuffix ) # WARNING: interpreted on the fly\n"""
    sout+="""echo ${all[$i]}\n"""
    sout+="""prefix=${all[$i]//$fastqsuffix}\n"""
    sout+="""\n"""
    sout+="""query1=$datadir/$prefix$fastqsuffix\n"""
    sout+="""query2=$datadir/$prefix${fastqsuffix/R1/R2}\n"""
    sout+="""bam=$bamdir/$prefix.bam\n"""
    sout+="""bam_sorted=$bamdir/$prefix"_sorted.bam"\n"""
    sout+="""bam_sorted_RG=$bamdir/$prefix"_sorted_RG.bam"\n"""
    sout+="""picard_bam=$bamdir/$prefix"_sorted_RG_rmd.bam"\n"""
    sout+="""indel_bam=$bamdir/$prefix"_sorted_RG_rmd_indelreal.bam"\n"""
    sout+="""bam_ready=$bamdir/$prefix"_sorted_RG_rmd_indelreal_BQSR.bam"\n"""
    sout+="""\n"""
    sout+="""echo "query1:" $query1\n"""
    sout+="""echo "query2:" $query2\n"""
    sout+="""\n"""
    sout+="""echo "************** Checking for directories ************"\n"""
    sout+="""for dir in $fastqcdir $resdir $bamdir $resvcfdir $samstatdir $picarddir ; do\n"""
    sout+="""    if  [ ! -e $dir ] ; then\n"""
    sout+="""    mkdir $dir\n"""
    sout+="""    fi\n"""
    sout+="""done\n"""
    sout+="""echo "************** JOB ENV *********************"\n"""
    sout+="""\n"""
    return sout    

def FastQC(param):
    """Outputs the FastQC.bash script"""

    for directory in ["SLURMlog","ScriptsDir"]:
        param[directory]=param[directory].rstrip(os.sep)+os.sep

    for directory in ["RawDataDir","ResultsDir"]:
        param[directory]=param[directory].rstrip(os.sep)

    out=param["ScriptsDir"]+"0--FastQC.bash"

    sout=header(param,"fastqc")
    sout+="""echo "************** Launching FastQC report PE1 ******************"\n"""
    sout+="""echo "@INPUT" $query1\n"""
    sout+="""echo "@OUTPUT"\n"""
    sout+="""$BIN/fastqc-0.10.1 --noextract --outdir=$fastqcdir $query1\n"""
    sout+="""\n"""
    sout+="""echo "************** Launching FastQC report PE2 ******************"\n"""
    sout+="""echo "@INPUT" $query2\n"""
    sout+="""echo "@OUTPUT"\n"""
    sout+="""$BIN/fastqc-0.10.1 --noextract --outdir=$fastqcdir $query2\n"""
    sout+="""\n"""
    sout+="""echo "************** Finished ******************"\n"""

    try:
        f=open(out,"w")
        f.write(sout)
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

    sout=header(param,"map")
    sout+="""echo "************** Preparing the @RG ReadGroups ************"\n"""
    sout+="""date=$(echo $prefix | awk 'BEGIN {FS="_"} {print $1}')\n"""
    sout+="""rgdt=${date//R} # This is Iso8601date (2000, but removed in 2004...)\n"""
    sout+="""info=$(echo $prefix | awk 'BEGIN {FS="_"} {OFS="_"} {print $1,$2,$3}') # in this experiment, one sample = one library, but several samples per run/lane\n"""
    sout+="""rgsm=$(echo $prefix | awk 'BEGIN {FS="_"} {print $2}') # in this experiment, one sample = one library, but several samples per run/lane\n"""
    sout+="""rglb=$(echo $prefix | awk 'BEGIN {FS="_"} {print $2}')\n"""
    sout+="""lane=$(echo $prefix | awk 'BEGIN {FS="_"} {print $3}')\n"""
    sout+="""rgid=$info #"_"$lane"_"$rglb"_"$rgdt\n"""
    sout+="""rgpl=illumina\n"""
    sout+="""rgpu=hiseq\n"""
    sout+="""rgcn=null\n"""
    sout+="""rgds="Instrument_HISEQ"\n"""
    sout+="""RG="@RG\\\\tID:"$rgid"\\\\tSM:"$rgsm"\\\\tLB:"$rglb"\\\\tPL:"$rgpl"\\\\tPU:"$rgpu"\\\\tCN:"$rgcn"\\\\tDS:"$rgds\n"""
    sout+="""echo "******************* Here is my @RG:" $RG "*******************"\n"""
    sout+="""\n"""
    sout+="""echo "************** Launching BWA-MEM alignments ******************"\n"""
    sout+="""echo "@INPUT" $query1\n"""
    sout+="""echo "@INPUT" $query2\n"""
    sout+="""echo "@OUTPUT" $bam\n"""
    sout+="""$BIN/bwa-0.7.7 mem \\\n"""
    sout+="""    -t $ncores \\\n"""
    sout+="""    -M \\\n"""
    sout+="""    -R $RG \\\n"""
    sout+="""    $db \\\n"""
    sout+="""    $query1 \\\n"""
    sout+="""    $query2 | $BIN/samtools-0.1.18 view -bSho $bam -\n"""
    sout+="""#  -M: mark shorter split hits as secondary (for Picard/GATK compatibility)\n"""
    sout+="""\n"""
    sout+="""echo "************** Launching samtools > sorted BAM ******************"\n"""
    sout+="""echo "@INPUT" $bam\n"""
    sout+="""echo "@OUTPUT" $bam_sorted\n"""
    sout+="""$BIN/samtools-0.1.18 sort $bam ${bam_sorted%.bam} # .bam suffix is appended by samtools sort...\n"""
    sout+="""rm $bam # the sorted version is enough for (and demanded by) most tools       \n"""
    sout+="""\n"""
    sout+="""echo "************** Launching samstat > BAM.html ******************"         \n"""
    sout+="""echo "@INPUT" $bam_sorted\n"""
    sout+="""echo "@OUTPUT" $bam_sorted".html"\n"""
    sout+="""$BIN/samstat-1.09 $bam_sorted # can be run before or after sorting\n"""
    sout+="""mv $bam_sorted.html $samstatdir/.\n"""
    sout+="""\n"""
    sout+="""echo "************** Launching Picard MarkDuplicates ******************"\n"""
    sout+="""rm $bam_sorted_RG\n"""
    sout+="""ln -s $bam_sorted $bam_sorted_RG\n"""
    sout+="""echo "@INPUT" $bam_sorted_RG\n"""
    sout+="""echo "@OUTPUT" $picard_bam ${picard_bam/.bam/.rmd} \n"""
    sout+="""$JAVAcustom $SRC/picard-tools-1.73/MarkDuplicates.jar \\\n"""
    sout+="""    INPUT=$bam_sorted_RG  \\\n"""
    sout+="""    OUTPUT=$picard_bam \\\n"""
    sout+="""    METRICS_FILE=${picard_bam/.bam/.rmd} \\\n"""
    sout+="""    REMOVE_DUPLICATES=false \\\n"""
    sout+="""    ASSUME_SORTED=true \\\n"""
    sout+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    sout+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    sout+="""\n"""
    sout+="""echo "************** Launching samtools > index sorted BAM ******************"\n"""
    sout+="""echo "@INPUT" $picard_bam\n"""
    sout+="""echo "@OUTPUT" $picard_bam".bai"\n"""
    sout+="""$BIN/samtools-0.1.18 index $picard_bam\n"""
    sout+="""\n"""
    sout+="""echo "************** Launching GATK RealignerTargetCreator ******************"\n"""
    sout+="""echo "@INPUT" $picard_bam\n"""
    sout+="""echo "@OUTPUT" $picard_bam".intervals"\n"""
    sout+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    sout+="""    -T RealignerTargetCreator \\\n"""
    sout+="""    -nt $nthreads \\\n"""
    sout+="""    -R $ref \\\n"""
    sout+="""    -I $picard_bam \\\n"""
    sout+="""    -o ${picard_bam/.bam/.intervals} \\\n"""
    sout+="""    --known $vcfdir"/1000G_phase1.indels.hg19_sorted.vcf" \\\n"""
    sout+="""    --known $vcfdir"/Mills_and_1000G_gold_standard.indels.hg19_sorted.vcf"\n"""
    sout+="""\n"""
    if param.has_key("noPhoneHome"):
        sout=noPhoneHome(param,sout)
    sout+="""echo "************** Launching GATK IndelRealigner ******************"\n"""
    sout+="""echo "@INPUT" $picard_bam\n"""
    sout+="""echo "@INPUT" $picard_bam".intervals"\n"""
    sout+="""echo "@OUTPUT" $indel_bam\n"""
    sout+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    sout+="""    -T IndelRealigner \\\n"""
    sout+="""    -R $ref \\\n"""
    sout+="""    -I $picard_bam \\\n"""
    sout+="""    -targetIntervals ${picard_bam/.bam/.intervals} \\\n"""
    sout+="""    -o $indel_bam \\\n"""
    sout+="""    --knownAlleles $vcfdir"/1000G_phase1.indels.hg19_sorted.vcf" \\\n"""
    sout+="""    --knownAlleles $vcfdir"/Mills_and_1000G_gold_standard.indels.hg19_sorted.vcf"\n"""
    sout+="""\n"""
    if param.has_key("noPhoneHome"):
        sout=noPhoneHome(param,sout)
    sout+="""echo "************** Launching samtools > index sorted BAM ******************"\n"""
    sout+="""echo "@INPUT" $indel_bam".bai"\n"""
    sout+="""echo "@OUTPUT" $indel_bam\n"""
    sout+="""$BIN/samtools-0.1.18 index $indel_bam\n"""
    sout+="""\n"""
    sout+="""echo "************** Launching GATK BaseRecalibrator 1st pass ******************"\n"""
    sout+="""echo "@INPUT" $indel_bam\n"""
    sout+="""echo "@OUTPUT" ${indel_bam/.bam/_recal1.grp} \n"""
    sout+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    sout+="""    -T BaseRecalibrator \\\n"""
    sout+="""    -nct $ncores \\\n"""
    sout+="""    -I $indel_bam \\\n"""
    sout+="""    -R $ref \\\n"""
    if param["AnalysisMode"]=="EXOME":
        sout+="""    -L $targets \\\n"""
    sout+="""    -knownSites $vcfdir"/dbsnp_137.hg19_sorted.vcf" \\\n"""
    sout+="""    -knownSites $vcfdir"/Mills_and_1000G_gold_standard.indels.hg19_sorted.vcf" \\\n"""
    sout+="""    -knownSites $vcfdir"/1000G_phase1.indels.hg19_sorted.vcf" \\\n"""
    sout+="""    -o ${indel_bam/.bam/_recal1.grp} \n"""
    sout+="""\n"""
    if param.has_key("noPhoneHome"):
        sout=noPhoneHome(param,sout)
    sout+="""echo "************** Launching GATK PrintReads ******************"\n"""
    sout+="""echo "@INPUT" $indel_bam\n"""
    sout+="""echo "@INPUT" ${indel_bam/.bam/_recal1.grp} \n"""
    sout+="""echo "@OUTPUT" $bam_ready\n"""
    sout+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    sout+="""    -T PrintReads \\\n"""
    sout+="""    -nct $ncores \\\n"""
    sout+="""    -R $ref \\\n"""
    if param["AnalysisMode"]=="EXOME":
        sout+="""    -L $targets \\\n"""
    sout+="""    -I $indel_bam \\\n"""
    sout+="""    -BQSR ${indel_bam/.bam/_recal1.grp} \\\n"""
    sout+="""    -o $bam_ready\n"""
    sout+="""\n"""
    if param.has_key("noPhoneHome"):
        sout=noPhoneHome(param,sout)
    sout+="""echo "************** Launching GATK BaseRecalibrator 2nd pass ******************"\n"""
    sout+="""echo "@INPUT" $indel_bam\n"""
    sout+="""echo "@INPUT" ${indel_bam/.bam/_recal1.grp} \n"""
    sout+="""echo "@OUTPUT" ${indel_bam/.bam/_recal2.grp} \n"""
    sout+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    sout+="""    -T BaseRecalibrator \\\n"""
    sout+="""    -BQSR ${indel_bam/.bam/_recal1.grp} \\\n"""
    sout+="""    -nct $ncores \\\n"""
    sout+="""    -I $indel_bam \\\n"""
    sout+="""    -R $ref \\\n"""
    if param["AnalysisMode"]=="EXOME":
        sout+="""    -L $targets \\\n"""
    sout+="""    -knownSites $vcfdir"/dbsnp_137.hg19_sorted.vcf" \\\n"""
    sout+="""    -knownSites $vcfdir"/Mills_and_1000G_gold_standard.indels.hg19_sorted.vcf" \\\n"""
    sout+="""    -knownSites $vcfdir"/1000G_phase1.indels.hg19_sorted.vcf" \\\n"""
    sout+="""    -o ${indel_bam/.bam/_recal2.grp} \n"""
    sout+="""\n"""
    if param.has_key("noPhoneHome"):
        sout=noPhoneHome(param,sout)
    sout+="""echo "************** Launching GATK AnalyzeCovariates > Plots ******************"\n"""
    sout+="""echo "@INPUT" ${indel_bam/.bam/_recal1.grp} \n"""
    sout+="""echo "@INPUT" ${indel_bam/.bam/_recal2.grp} \n"""
    sout+="""echo "@OUTPUT" ${indel_bam/.bam/_BQSR.csv} \n"""
    sout+="""echo "@OUTPUT" ${indel_bam/.bam/_BQSR.pdf} \n"""
    sout+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    sout+="""    -T AnalyzeCovariates \\\n"""
    sout+="""    -R $ref \\\n"""
    if param["AnalysisMode"]=="EXOME":
        sout+="""    -L $targets \\\n"""
    sout+="""    -before ${indel_bam/.bam/_recal1.grp} \\\n"""
    sout+="""    -after ${indel_bam/.bam/_recal2.grp} \\\n"""
    sout+="""    -csv ${indel_bam/.bam/_BQSR.csv} \\\n"""
    sout+="""    -plots ${indel_bam/.bam/_BQSR.pdf} \n"""
    sout+="""\n"""
    if param.has_key("noPhoneHome"):
        sout=noPhoneHome(param,sout)
    sout+="""echo "************** Launching samtools > index sorted BAM ******************"\n"""
    sout+="""echo "@INPUT" $bam_ready\n"""
    sout+="""echo "@OUTPUT" $bam_ready.bai\n"""
    sout+="""$BIN/samtools-0.1.18 index $bam_ready # .bai indexes positions on the reference genome and some tools can use this directly\n"""
    sout+="""\n"""
    sout+="""echo "************** Launching samstat > BAM.html ******************"\n"""
    sout+="""echo "@INPUT" $bam_ready\n"""
    sout+="""echo "@OUTPUT" $bam_ready.html\n"""
    sout+="""$BIN/samstat-1.09 $bam_ready\n"""
    sout+="""mv $bam_ready.html $samstatdir/.\n"""
    sout+="""\n"""
    sout+="""echo "************** Launching samtools idxstats ******************"\n"""
    sout+="""echo "@INPUT" $bam_ready\n"""
    sout+="""echo "@OUTPUT" ${bam_ready/.bam/.idxstats} \n"""
    sout+="""$BIN/samtools-0.1.18 idxstats $bam_ready > ${bam_ready/.bam/.idxstats} \n"""
    sout+="""\n"""
    sout+="""echo "************** Finished ******************"\n"""

    try:
        f=open(out,"w")
        f.write(sout)
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

    sout=header(param,"qc")
    sout+="""echo 'BAM for QC:' $bam_ready\n"""
    sout+="""\n"""
    sout+="""echo '************** Launching Samtools Flagstat ******************'\n"""
    sout+="""$BIN/samtools-0.1.18 flagstat $bam_ready >  ${bam_ready/.bam/.flagstat} \n"""
    sout+="""\n"""
    sout+="""echo '************** Launching Picard EstimateLibraryComplexity *****************'\n"""
    sout+="""$JAVAcustom $SRC/picard-tools-1.73/EstimateLibraryComplexity.jar \\\n"""
    sout+="""    INPUT=$bam_ready \\\n"""
    sout+="""    OUTPUT=${bam_ready/.bam/.lc} \\\n"""
    sout+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    sout+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    sout+="""#ASSUME_SORTED=true         \n"""
    sout+="""\n"""
    sout+="""echo '************** Launching Picard MeanQualityByCycle *****************'\n"""
    sout+="""$JAVAcustom $SRC/picard-tools-1.73/MeanQualityByCycle.jar \\\n"""
    sout+="""    INPUT=$bam_ready \\\n"""
    sout+="""    OUTPUT=$bam_ready.qbc \\\n"""
    sout+="""    CHART_OUTPUT=${bam_ready/.bam/_qbc.pdf} \\\n"""
    sout+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    sout+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    sout+="""\n"""
    sout+="""echo '************** Launching Picard QualityScoreDistribution ****************'\n"""
    sout+="""$JAVAcustom $SRC/picard-tools-1.73/QualityScoreDistribution.jar \\\n"""
    sout+="""    INPUT=$bam_ready \\\n"""
    sout+="""    OUTPUT=${bam_ready/.bam/.qsd} \\\n"""
    sout+="""    CHART_OUTPUT=${bam_ready/.bam/_qsd.pdf} \\\n"""
    sout+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    sout+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    sout+="""\n"""
    sout+="""echo '************** Launching Picard CollectGcBiasMetrics *****************'\n"""
    sout+="""$JAVAcustom $SRC/picard-tools-1.73/CollectGcBiasMetrics.jar \\\n"""
    sout+="""    INPUT=$bam_ready \\\n"""
    sout+="""    REFERENCE_SEQUENCE=$ref \\\n"""
    sout+="""    OUTPUT=${bam_ready/.bam/.gcb} \\\n"""
    sout+="""    CHART_OUTPUT=${bam_ready/.bam/_gcb.pdf} \\\n"""
    sout+="""    SUMMARY_OUTPUT=${bam_ready/.bam/_gcb.sum} \\\n"""
    sout+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    sout+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    sout+="""\n"""
    sout+="""echo '************** Launching Picard CollectInsertSizeMetrics *****************'\n"""
    sout+="""$JAVAcustom $SRC/picard-tools-1.73/CollectInsertSizeMetrics.jar \\\n"""
    sout+="""    INPUT=$bam_ready \\\n"""
    sout+="""    OUTPUT=${bam_ready/.bam/.ism} \\\n"""
    sout+="""    HISTOGRAM_FILE=${bam_ready/.bam/_ism.pdf} \\\n"""
    sout+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    sout+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    sout+="""\n"""
    sout+="""echo '************** Launching Picard CollectAlignmentSummaryMetrics ****************'\n"""
    sout+="""$JAVAcustom $SRC/picard-tools-1.73/CollectAlignmentSummaryMetrics.jar \\\n"""
    sout+="""    INPUT=$bam_ready \\\n"""
    sout+="""    OUTPUT=${bam_ready/.bam/.asm} \\\n"""
    sout+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    sout+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    sout+="""\n"""
    sout+="""echo '************** Launching IGVtools count **********'\n"""
    sout+="""$JAVAcustom $SRC/IGVTools-1.5.10/igvtools.jar count \\\n"""
    sout+="""    $bam_ready  \\\n"""
    sout+="""    ${bam_ready/.bam/.tdf} \\\n"""
    sout+="""    $refname\n"""
    sout+="""\n"""
    if param["AnalysisMode"]=="EXOME":
        sout+="""# WARNING: careful with the manifest and target files here!\n"""
        sout+="""echo '************** Launching Picard CalculateHSMetrics ******************'\n"""
        sout+="""$JAVAcustom $SRC/picard-tools-1.73/CalculateHsMetrics.jar \\\n"""
        sout+="""    BAIT_SET_NAME=$baitNames \\\n"""
        sout+="""    BAIT_INTERVALS=$baitsPicard \\\n"""
        sout+="""    TARGET_INTERVALS=$targetsPicard \\\n"""
        sout+="""    INPUT=$bam_ready \\\n"""
        sout+="""    OUTPUT=${bam_ready/.bam/.hsm} \\\n"""
        sout+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
        sout+="""    REFERENCE_SEQUENCE=$ref \\\n"""
        sout+="""    METRIC_ACCUMULATION_LEVEL=SAMPLE \\\n"""
        sout+="""    TMP_DIR=$GLOBALSCRATCH\n"""
        sout+="""\n"""
        sout+="""echo '************** Launching GATK DiagnoseTargets ******************'\n"""
        sout+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
        sout+="""    -T DiagnoseTargets \\\n"""
        sout+="""    -R $ref \\\n"""
        sout+="""    -I $bam_ready \\\n"""
        sout+="""    -o ${bam_ready/.bam/.diagt} \\\n"""
        sout+="""    -L $targets\n"""
        sout+="""\n"""
        if param.has_key("noPhoneHome"):
            sout=noPhoneHome(param,sout)
    sout+="""echo '************** Launching GATK DepthOfCoverage ******************'\n"""
    sout+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    sout+="""    -T DepthOfCoverage \\\n"""
    sout+="""    -R $ref \\\n"""
    sout+="""    -I $bam_ready \\\n"""
    if param["AnalysisMode"]=="EXOME":
        sout+="""    -L $targets\\\n"""
    sout+="""    -o ${bam_ready/.bam/.cov} \\\n"""
    sout+="""    -pt sample \\\n"""
    sout+="""    --omitDepthOutputAtEachBase \\\n"""
    sout+="""    -ct 20 \n"""
    sout+="""\n"""
    if param.has_key("noPhoneHome"):
        sout=noPhoneHome(param,sout)
    sout+="""echo '************** Finished ******************'\n"""

    try:
        f=open(out,"w")
        f.write(sout)
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

    sout=header(param,"gvcf")
    sout+="""middfix='_HC3.1' #@@@ fill this in\n"""
    sout+="""\n"""
    sout+="""echo 'BAM for HC:' $bam_ready\n"""
    sout+="""\n"""
    sout+="""echo '************** Launching GATK HaplotypeCaller ******************'\n"""
    sout+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    sout+="""    -T HaplotypeCaller \\\n"""
    sout+="""    -R $ref \\\n"""
    if param["AnalysisMode"]=="EXOME":
        sout+="""    -L $targets \\\n"""
    sout+="""    -I $bam_ready \\\n"""
    sout+="""    -o ${resvcfdir}$( echo $bam_ready | sed -e 's#'$bamdir'##; s#.bam##' )${middfix}".gvcf" \\\n"""
    sout+="""    --dbsnp $vcfdir'/dbsnp_137.hg19_sorted.vcf' \\\n"""
    sout+="""    --emitRefConfidence GVCF \\\n"""
    sout+="""    --variant_index_type LINEAR \\\n"""
    sout+="""    --variant_index_parameter 128000 \\\n"""
    sout+="""    --pair_hmm_implementation VECTOR_LOGLESS_CACHING \n"""
    sout+="""\n"""
    if param.has_key("noPhoneHome"):
        sout=noPhoneHome(param,sout)
    sout+="""echo '************** Finished ******************'\n"""

    try:
        f=open(out,"w")
        f.write(sout)
        f.close()
    except IOError, e:
        print "File not found: [", out, "]"
        sys.exit(2)
    return

def GenotypingAndRecalibrating(param):
    """Outputs the GenotypingAndRecalibrating.bash script"""

    for directory in ["SLURMlog","ScriptsDir"]:
        param[directory]=param[directory].rstrip(os.sep)+os.sep

    for directory in ["RawDataDir","ResultsDir"]:
        param[directory]=param[directory].rstrip(os.sep)

    out=param["ScriptsDir"]+"3--GenotypingAndRecalibrating.bash"

    sout=header(param,"vcf")
    sout+="""middfix='_HC3.1'\n"""
    sout+="""\n"""
    sout+="""cd $resvcfdir\n"""
    sout+="""prefix='all'\n"""
    sout+="""\n"""
    sout+="""gvcflist=$resdir"/gvcfs"$middfix".list"\n"""
    sout+="""\ls $resvcfdir/*.gvcf > $gvcflist\n"""
    sout+="""echo "************* Launching GATK GenotypeGVCFs ******************"\n"""
    sout+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    sout+="""    -T GenotypeGVCFs \\\n"""
    sout+="""    -V $gvcflist \\\n"""
    if param["AnalysisMode"]=="EXOME":
        sout+="""    -L $targets\\\n"""
    sout+="""    -o $resvcfdir/$prefix$middfix.vcf \\\n"""
    sout+="""    -R $ref \\\n"""
    sout+="""    -nt $nthreads \\\n"""
    sout+="""    --dbsnp $vcfdir"/dbsnp_137.hg19_sorted.vcf" \n"""
    sout+="""\n"""
    if param.has_key("noPhoneHome"):
        sout=noPhoneHome(param,sout)
    sout+="""echo "************** Launching GATK VariantRecalibrator -- SNP pass ******************"\n"""
    sout+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    sout+="""    -T VariantRecalibrator \\\n"""
    sout+="""    -nt $nthreads \\\n"""
    sout+="""    -R $ref \\\n"""
    sout+="""    -input $resvcfdir/$prefix$middfix.vcf \\\n"""
    if param["AnalysisMode"]=="EXOME":
        sout+="""    -L $targets \\\n"""
    sout+="""    -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 $vcfdir"/hapmap_3.3.hg19_sorted.vcf" \\\n"""
    sout+="""    -resource:omni,VCF,known=false,training=true,truth=true,prior=12.0 $vcfdir"/1000G_omni2.5.hg19_sorted.vcf" \\\n"""
    sout+="""    -resource:1000G,known=false,training=true,truth=false,prior=10.0 $vcfdir"/1000G_phase1.snps.high_confidence.hg19_sorted.vcf" \\\n"""
    sout+="""    -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=2.0 $vcfdir"/dbsnp_137.hg19_sorted.vcf" \\\n"""
    sout+="""    -an QD -an MQRankSum -an ReadPosRankSum -an FS \\\n"""
    sout+="""    -mode SNP \\\n"""
    sout+="""    -recalFile $resvcfdir/$prefix$middfix.recal_snps \\\n"""
    sout+="""    -tranchesFile $resvcfdir/$prefix$middfix.tranches_snps \\\n"""
    sout+="""    -rscriptFile $resvcfdir/$prefix$middfix_VQSR_SNP.R \\\n"""
    sout+="""    --maxGaussians 4\n"""
    sout+="""\n"""
    if param.has_key("noPhoneHome"):
        sout=noPhoneHome(param,sout)
    sout+="""echo "************** Launching GATK ApplyRecalibration -- SNP pass ******************"\n"""
    sout+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    sout+="""    -T ApplyRecalibration \\\n"""
    sout+="""    -nt $nthreads \\\n"""
    sout+="""    -R $ref \\\n"""
    sout+="""    -input $resvcfdir/$prefix$middfix.vcf \\\n"""
    if param["AnalysisMode"]=="EXOME":
        sout+="""    -L $targets \\\n"""
    sout+="""    --ts_filter_level 99.0 \\\n"""
    sout+="""    -recalFile $resvcfdir/$prefix$middfix.recal_snps \\\n"""
    sout+="""    -tranchesFile $resvcfdir/$prefix$middfix.tranches_snps \\\n"""
    sout+="""    -mode SNP \\\n"""
    sout+="""    -o $resvcfdir/$prefix$middfix_recal_step1.vcf \n"""
    sout+="""\n"""
    if param.has_key("noPhoneHome"):
        sout=noPhoneHome(param,sout)
    sout+="""echo "************** Launching GATK VariantRecalibrator -- INDEL pass ******************"\n"""
    sout+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    sout+="""    -T VariantRecalibrator \\\n"""
    sout+="""    -nt $nthreads \\\n"""
    sout+="""    -R $ref \\\n"""
    sout+="""    -input $resvcfdir/$prefix$middfix_recal_step1.vcf \\\n"""
    if param["AnalysisMode"]=="EXOME":
        sout+="""    -L $targets \\\n"""
    sout+="""    -resource:mills,VCF,known=false,training=true,truth=true,prior=12.0 $vcfdir"/Mills_and_1000G_gold_standard.indels.hg19_sorted.vcf" \\\n"""
    sout+="""    -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=2.0 $vcfdir"/dbsnp_137.hg19_sorted.vcf" \\\n"""
    sout+="""    -an FS -an ReadPosRankSum -an MQRankSum \\\n"""
    sout+="""    -mode INDEL \\\n"""
    sout+="""    -recalFile $resvcfdir/$prefix$middfix.recal_indels \\\n"""
    sout+="""    -tranchesFile $resvcfdir/$prefix$middfix.tranches_indels \\\n"""
    sout+="""    -rscriptFile $resvcfdir/$prefix$middfix_VQSR_INDEL.R \\\n"""
    sout+="""    --maxGaussians 4\n"""
    sout+="""\n"""
    if param.has_key("noPhoneHome"):
        sout=noPhoneHome(param,sout)
    sout+="""echo "************** Launching GATK ApplyRecalibration -- INDEL pass ******************"\n"""
    sout+="""$JAVAcustom $SRC/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \\\n"""
    sout+="""    -T ApplyRecalibration \\\n"""
    sout+="""    -nt $nthreads \\\n"""
    sout+="""    -R $ref \\\n"""
    sout+="""    -input $resvcfdir/$prefix$middfix_recal_step1.vcf \\\n"""
    if param["AnalysisMode"]=="EXOME":
        sout+="""    -L $targets \\\n"""
    sout+="""    --ts_filter_level 99.0 \\\n"""
    sout+="""    -tranchesFile $resvcfdir/$prefix$middfix.tranches_indels \\\n"""
    sout+="""    -recalFile $resvcfdir/$prefix$middfix.recal_indels \\\n"""
    sout+="""    -mode INDEL \\\n"""
    sout+="""    -o $resvcfdir/$prefix$middfix_recal_final.vcf \n"""
    sout+="""\n"""
    if param.has_key("noPhoneHome"):
        sout=noPhoneHome(param,sout)
    sout+="""echo "************** Launching Variant Effect Predictor (VEP - Ensembl) annotation ******************"\n"""
    sout+="""perl $SRC/ensembl-tools-release-77/scripts/variant_effect_predictor/variant_effect_predictor.pl \\\n"""
    sout+="""    -i $resvcfdir/$prefix$middfix_recal_final.vcf \\\n"""
    sout+="""    --cache --everything --vcf --force_overwrite \\\n"""
    sout+="""    --output_file $resvcfdir/$prefix$middfix_recal_final_VEP.vcf \\\n"""
    sout+="""    --stats_file $resvcfdir/$prefix$middfix_recal_final_VEP_summary.html \n"""
    sout+="""\n"""
    sout+="""randomNumber=$RANDOM \n"""
    sout+="""$SRC/vcfsorter.pl $dict $resvcfdir/$prefix$middfix_recal_final_VEP.vcf > swap_$randomNumber \n"""
    sout+="""mv swap_$randomNumber $resvcfdir/$prefix$middfix_recal_final_VEP.vcf \n"""
    sout+="""\n"""
    sout+="""echo "************** Finished ******************"\n"""

    try:
        f=open(out,"w")
        f.write(sout)
        f.close()
    except IOError, e:
        print "File not found: [", out, "]"
        sys.exit(2)
    return
