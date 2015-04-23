
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
    sout+="""#SBATCH --output="""+param["SLURMlog"]+"""slurm-%j-"""+step+""".out\n"""
    sout+="""#\n"""
    sout+="""#SBATCH --ntasks=1                                               #@@@ fill with appropriate value: here 1 task\n"""
    sout+="""#SBATCH --cpus-per-task=1                                        #@@@ fill with appropriate value (max on HMEM=25): here 1 core per task\n"""
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
    sout+="""JAVAcustom=$BIN"/java-1.7.0_71 -Djava.io.tmpdir="$GLOBALSCRATCH"/tmp -Xmx"$mem"g -XX:ParallelGCThreads="$ncores" -jar"\n"""
    sout+="""FASTQC=$BIN/fastqc_v0.11.2\n"""
    sout+="""BWA=$BIN/bwa-0.7.10\n"""
    sout+="""SAMTOOLS=$BIN/samtools-1.1\n"""
    sout+="""SAMSTAT=$BIN/samstat-1.5\n"""
    sout+="""PICARD=$JAVAcustom" "$SRC/picard-tools-1.78\n"""
    sout+="""GATK=$JAVAcustom" "$SRC/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar\n"""
    sout+="""PERL=/usr/bin/perl\n"""
    sout+="""VEP=$PERL" "$SRC/ensembl-tools-release-78/scripts/variant_effect_predictor/variant_effect_predictor.pl\n"""
    sout+="""VCFSORTER=$BIN/vcfsorter.pl\n"""
    sout+="""BGZIP=$BIN/bgzip-1.1\n"""
    sout+="""TABIX=$BIN/tabix-1.1\n"""
    sout+="""SNPSIFT=$JAVAcustom" "$BIN/snpSift_4.0_E_2014-09-13.jar\n"""
    sout+="""module av R\n"""
    sout+="""module load R/3.1.2-goolf-1.6.10\n"""
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
    sout+="""datadir="""+param["RawDataDir"]+"""\n"""
    sout+="""resdir="""+param["ResultsDir"]+"""\n"""
    sout+="""fastqcdir=$resdir/QC/fastqc\n"""
    sout+="""bamdir=$resdir/bam\n"""
    sout+="""resvcfdir=$resdir/vcf\n"""
    sout+="""samstatdir=$resdir/QC/samstat\n"""
    sout+="""fastqsuffix='"""+param["FastqGzSuffixPE"]+"""'\n"""
    sout+="""refname="""+param["ReferenceAssembly"]+"""\n"""
    sout+="""ref=$GLOBALSCRATCH/genomes/homo_sapiens/hg19/genome/HG19.fasta\n"""
    sout+="""dict=$GLOBALSCRATCH/genomes/homo_sapiens/hg19/genome/HG19.dict\n"""
    sout+="""db=$GLOBALSCRATCH/genomes/homo_sapiens/hg19/bwa_hash/HG19.fasta\n"""
    sout+="""vcfdir=$GLOBALSCRATCH/genomes/homo_sapiens/hg19/variation\n"""
    sout+="""clinvarVCF=$GLOBALSCRATCH/genomes/homo_sapiens/hg19/allelic_frequencies/clinvar_20141202.vcf.gz\n"""
    sout+="""exacVCF=$GLOBALSCRATCH/genomes/homo_sapiens/hg19/allelic_frequencies/ExAC.r0.3.sites.vep.ExAC.vcf.gz\n"""
    sout+="""kgVCF=$GLOBALSCRATCH/genomes/homo_sapiens/hg19/allelic_frequencies/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz\n"""
    if param["AnalysisMode"]=="EXOME":
        sout+="""echo "*************** EXOME Analysis Mode *****************"\n"""
    elif param["AnalysisMode"]=="PANEL":
        sout+="""echo "*************** PANEL Analysis Mode *****************"\n"""
    if param["AnalysisMode"]=="EXOME" or param["AnalysisMode"]=="PANEL":
        sout+="""targets="""+param["TargetFile"]+"""\n"""
        sout+="""baitsPicard="""+param["BaitsFilePicard"]+"""\n"""
        sout+="""targetsPicard="""+param["TargetFilePicard"]+"""\n"""
    else:
        sout+="""echo "*************** GENOME Analysis Mode *****************"\n"""
    sout+="""echo "datadir:" $datadir\n"""
    sout+="""echo "fastqcdir:" $fastqcdir\n"""
    sout+="""echo "resdir:" $resdir\n"""
    sout+="""echo "bamdir:" $bamdir\n"""
    sout+="""echo "resvcfdir:" $resvcfdir\n"""
    sout+="""echo "samstatdir:" $samstatdir\n"""
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
    sout+="""for dir in $fastqcdir $resdir $bamdir $resvcfdir $samstatdir ; do\n"""
    sout+="""    if  [ ! -e $dir ] ; then\n"""
    sout+="""    mkdir -p $dir\n"""
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
    sout+="""$FASTQC --noextract --outdir=$fastqcdir $query1\n"""
    sout+="""\n"""
    sout+="""echo "************** Launching FastQC report PE2 ******************"\n"""
    sout+="""echo "@INPUT" $query2\n"""
    sout+="""echo "@OUTPUT"\n"""
    sout+="""$FASTQC --noextract --outdir=$fastqcdir $query2\n"""
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
    sout+="""$BWA mem \\\n"""
    sout+="""    -t $ncores \\\n"""
    sout+="""    -M \\\n"""
    sout+="""    -R $RG \\\n"""
    sout+="""    $db \\\n"""
    sout+="""    $query1 \\\n"""
    sout+="""    $query2 | $SAMTOOLS view -bSho $bam -\n"""
    sout+="""#  -M: mark shorter split hits as secondary (for Picard/GATK compatibility)\n"""
    sout+="""\n"""
    sout+="""echo "************** Launching samtools > sorted BAM ******************"\n"""
    sout+="""echo "@INPUT" $bam\n"""
    sout+="""echo "@OUTPUT" $bam_sorted\n"""
    sout+="""$SAMTOOLS sort $bam ${bam_sorted%.bam} # .bam suffix is appended by samtools sort...\n"""
    sout+="""rm $bam # the sorted version is enough for (and demanded by) most tools       \n"""
    sout+="""\n"""
    sout+="""echo "************** Launching samstat > BAM.html ******************"         \n"""
    sout+="""echo "@INPUT" $bam_sorted\n"""
    sout+="""echo "@OUTPUT" $bam_sorted".samstat.html"\n"""
    sout+="""$SAMSTAT $bam_sorted # can be run before or after sorting\n"""
    sout+="""mv $bam_sorted.samstat.html $samstatdir/.\n"""
    sout+="""\n"""
    sout+="""echo "************** Launching Picard MarkDuplicates ******************"\n"""
    sout+="""rm $bam_sorted_RG\n"""
    sout+="""ln -s $bam_sorted $bam_sorted_RG\n"""
    sout+="""echo "@INPUT" $bam_sorted_RG\n"""
    sout+="""echo "@OUTPUT" $picard_bam ${picard_bam/.bam/.rmd} \n"""
    sout+="""$PICARD/MarkDuplicates.jar \\\n"""
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
    sout+="""$SAMTOOLS index $picard_bam\n"""
    sout+="""\n"""
    sout+="""echo "************** Launching GATK RealignerTargetCreator ******************"\n"""
    sout+="""echo "@INPUT" $picard_bam\n"""
    sout+="""echo "@OUTPUT" $picard_bam".intervals"\n"""
    sout+="""$GATK \\\n"""
    sout+="""    -T RealignerTargetCreator \\\n"""
    sout+="""    -nt $nthreads \\\n"""
    sout+="""    -R $ref \\\n"""
    sout+="""    -I $picard_bam \\\n"""
    sout+="""    -o ${picard_bam/.bam/.intervals} \\\n"""
    sout+="""    --known $vcfdir"/1000G_phase1.indels.hg19.sites.vcf.gz" \\\n"""
    sout+="""    --known $vcfdir"/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"\n"""
    sout+="""\n"""
    if param.has_key("noPhoneHome"):
        sout=noPhoneHome(param,sout)
    sout+="""echo "************** Launching GATK IndelRealigner ******************"\n"""
    sout+="""echo "@INPUT" $picard_bam\n"""
    sout+="""echo "@INPUT" $picard_bam".intervals"\n"""
    sout+="""echo "@OUTPUT" $indel_bam\n"""
    sout+="""$GATK \\\n"""
    sout+="""    -T IndelRealigner \\\n"""
    sout+="""    -R $ref \\\n"""
    sout+="""    -I $picard_bam \\\n"""
    sout+="""    -targetIntervals ${picard_bam/.bam/.intervals} \\\n"""
    sout+="""    -o $indel_bam \\\n"""
    sout+="""    --knownAlleles $vcfdir"/1000G_phase1.indels.hg19.sites.vcf.gz" \\\n"""
    sout+="""    --knownAlleles $vcfdir"/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"\n"""
    sout+="""\n"""
    if param.has_key("noPhoneHome"):
        sout=noPhoneHome(param,sout)
    sout+="""echo "************** Launching samtools > index sorted BAM ******************"\n"""
    sout+="""echo "@INPUT" $indel_bam".bai"\n"""
    sout+="""echo "@OUTPUT" $indel_bam\n"""
    sout+="""$SAMTOOLS index $indel_bam\n"""
    sout+="""\n"""
    sout+="""echo "************** Launching GATK BaseRecalibrator 1st pass ******************"\n"""
    sout+="""echo "@INPUT" $indel_bam\n"""
    sout+="""echo "@OUTPUT" ${indel_bam/.bam/_recal1.grp} \n"""
    sout+="""$GATK \\\n"""
    sout+="""    -T BaseRecalibrator \\\n"""
    sout+="""    -nct $ncores \\\n"""
    sout+="""    -I $indel_bam \\\n"""
    sout+="""    -R $ref \\\n"""
    if param["AnalysisMode"]=="EXOME" or param["AnalysisMode"]=="PANEL":
        sout+="""    -L $targets \\\n"""
    sout+="""    -knownSites $vcfdir"/dbsnp_138.hg19.vcf.gz" \\\n"""
    sout+="""    -knownSites $vcfdir"/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz" \\\n"""
    sout+="""    -knownSites $vcfdir"/1000G_phase1.indels.hg19.sites.vcf.gz" \\\n"""
    sout+="""    -o ${indel_bam/.bam/_recal1.grp} \n"""
    sout+="""\n"""
    if param.has_key("noPhoneHome"):
        sout=noPhoneHome(param,sout)
    sout+="""echo "************** Launching GATK PrintReads ******************"\n"""
    sout+="""echo "@INPUT" $indel_bam\n"""
    sout+="""echo "@INPUT" ${indel_bam/.bam/_recal1.grp} \n"""
    sout+="""echo "@OUTPUT" $bam_ready\n"""
    sout+="""$GATK \\\n"""
    sout+="""    -T PrintReads \\\n"""
    sout+="""    -nct $ncores \\\n"""
    sout+="""    -R $ref \\\n"""
    if param["AnalysisMode"]=="EXOME" or param["AnalysisMode"]=="PANEL":
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
    sout+="""$GATK \\\n"""
    sout+="""    -T BaseRecalibrator \\\n"""
    sout+="""    -BQSR ${indel_bam/.bam/_recal1.grp} \\\n"""
    sout+="""    -nct $ncores \\\n"""
    sout+="""    -I $indel_bam \\\n"""
    sout+="""    -R $ref \\\n"""
    if param["AnalysisMode"]=="EXOME" or param["AnalysisMode"]=="PANEL":
        sout+="""    -L $targets \\\n"""
    sout+="""    -knownSites $vcfdir"/dbsnp_138.hg19.vcf.gz" \\\n"""
    sout+="""    -knownSites $vcfdir"/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz" \\\n"""
    sout+="""    -knownSites $vcfdir"/1000G_phase1.indels.hg19.sites.vcf.gz" \\\n"""
    sout+="""    -o ${indel_bam/.bam/_recal2.grp} \n"""
    sout+="""\n"""
    if param.has_key("noPhoneHome"):
        sout=noPhoneHome(param,sout)
    sout+="""echo "************** Launching GATK AnalyzeCovariates > Plots ******************"\n"""
    sout+="""echo "@INPUT" ${indel_bam/.bam/_recal1.grp} \n"""
    sout+="""echo "@INPUT" ${indel_bam/.bam/_recal2.grp} \n"""
    sout+="""echo "@OUTPUT" ${indel_bam/.bam/_BQSR.csv} \n"""
    sout+="""echo "@OUTPUT" ${indel_bam/.bam/_BQSR.pdf} \n"""
    sout+="""$GATK \\\n"""
    sout+="""    -T AnalyzeCovariates \\\n"""
    sout+="""    -R $ref \\\n"""
    if param["AnalysisMode"]=="EXOME" or param["AnalysisMode"]=="PANEL":
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
    sout+="""$SAMTOOLS index $bam_ready # .bai indexes positions on the reference genome and some tools can use this directly\n"""
    sout+="""\n"""
    sout+="""echo "************** Launching samstat > BAM.html ******************"\n"""
    sout+="""echo "@INPUT" $bam_ready\n"""
    sout+="""echo "@OUTPUT" $bam_ready.samstat.html\n"""
    sout+="""$SAMSTAT $bam_ready\n"""
    sout+="""mv $bam_ready.samstat.html $samstatdir/.\n"""
    sout+="""\n"""
    sout+="""echo "************** Launching samtools idxstats ******************"\n"""
    sout+="""echo "@INPUT" $bam_ready\n"""
    sout+="""echo "@OUTPUT" ${bam_ready/.bam/.idxstats} \n"""
    sout+="""$SAMTOOLS idxstats $bam_ready > ${bam_ready/.bam/.idxstats} \n"""
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
    sout+="""$SAMTOOLS flagstat $bam_ready >  ${bam_ready/.bam/.flagstat} \n"""
    sout+="""\n"""
    sout+="""echo '************** Launching Picard EstimateLibraryComplexity *****************'\n"""
    sout+="""$PICARD/EstimateLibraryComplexity.jar \\\n"""
    sout+="""    INPUT=$bam_ready \\\n"""
    sout+="""    OUTPUT=${bam_ready/.bam/.lc} \\\n"""
    sout+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    sout+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    sout+="""#ASSUME_SORTED=true         \n"""
    sout+="""\n"""
    sout+="""echo '************** Launching Picard MeanQualityByCycle *****************'\n"""
    sout+="""$PICARD/MeanQualityByCycle.jar \\\n"""
    sout+="""    INPUT=$bam_ready \\\n"""
    sout+="""    OUTPUT=$bam_ready.qbc \\\n"""
    sout+="""    CHART_OUTPUT=${bam_ready/.bam/_qbc.pdf} \\\n"""
    sout+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    sout+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    sout+="""\n"""
    sout+="""echo '************** Launching Picard QualityScoreDistribution ****************'\n"""
    sout+="""$PICARD/QualityScoreDistribution.jar \\\n"""
    sout+="""    INPUT=$bam_ready \\\n"""
    sout+="""    OUTPUT=${bam_ready/.bam/.qsd} \\\n"""
    sout+="""    CHART_OUTPUT=${bam_ready/.bam/_qsd.pdf} \\\n"""
    sout+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    sout+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    sout+="""\n"""
    sout+="""echo '************** Launching Picard CollectGcBiasMetrics *****************'\n"""
    sout+="""$PICARD/CollectGcBiasMetrics.jar \\\n"""
    sout+="""    INPUT=$bam_ready \\\n"""
    sout+="""    REFERENCE_SEQUENCE=$ref \\\n"""
    sout+="""    OUTPUT=${bam_ready/.bam/.gcb} \\\n"""
    sout+="""    CHART_OUTPUT=${bam_ready/.bam/_gcb.pdf} \\\n"""
    sout+="""    SUMMARY_OUTPUT=${bam_ready/.bam/_gcb.sum} \\\n"""
    sout+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    sout+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    sout+="""\n"""
    sout+="""echo '************** Launching Picard CollectInsertSizeMetrics *****************'\n"""
    sout+="""$PICARD/CollectInsertSizeMetrics.jar \\\n"""
    sout+="""    INPUT=$bam_ready \\\n"""
    sout+="""    OUTPUT=${bam_ready/.bam/.ism} \\\n"""
    sout+="""    HISTOGRAM_FILE=${bam_ready/.bam/_ism.pdf} \\\n"""
    sout+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    sout+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    sout+="""\n"""
    sout+="""echo '************** Launching Picard CollectAlignmentSummaryMetrics ****************'\n"""
    sout+="""$PICARD/CollectAlignmentSummaryMetrics.jar \\\n"""
    sout+="""    INPUT=$bam_ready \\\n"""
    sout+="""    OUTPUT=${bam_ready/.bam/.asm} \\\n"""
    sout+="""    VALIDATION_STRINGENCY=LENIENT \\\n"""
    sout+="""    TMP_DIR=$GLOBALSCRATCH\n"""
    sout+="""\n"""
    sout+="""echo '************** Launching IGVtools count **********'\n"""
    sout+="""$JAVAcustom $SRC/IGVTools-1.5.10/igvtools.jar count \\\n"""
    sout+="""    $bam_ready  \\\n"""
    sout+="""    ${bam_ready/.bam/.tdf},${bam_ready/.bam/.wig} \\\n"""
    sout+="""    $refname\n"""
    sout+="""\n"""
    if param["AnalysisMode"]=="EXOME" or param["AnalysisMode"]=="PANEL":
        sout+="""# WARNING: careful with the manifest and target files here!\n"""
        sout+="""echo '************** Launching Picard CalculateHSMetrics ******************'\n"""
        sout+="""$PICARD/CalculateHsMetrics.jar \\\n"""
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
        sout+="""$GATK \\\n"""
        sout+="""    -T DiagnoseTargets \\\n"""
        sout+="""    -R $ref \\\n"""
        sout+="""    -I $bam_ready \\\n"""
        sout+="""    -o ${bam_ready/.bam/.diagt} \\\n"""
        sout+="""    -L $targets\n"""
        sout+="""\n"""
        if param.has_key("noPhoneHome"):
            sout=noPhoneHome(param,sout)
    sout+="""echo '************** Launching GATK DepthOfCoverage ******************'\n"""
    sout+="""$GATK \\\n"""
    sout+="""    -T DepthOfCoverage \\\n"""
    sout+="""    -R $ref \\\n"""
    sout+="""    -I $bam_ready \\\n"""
    if param["AnalysisMode"]=="EXOME" or param["AnalysisMode"]=="PANEL":
        sout+="""    -L $targets\\\n"""
    sout+="""    -o ${bam_ready/.bam/.cov} \\\n"""
    sout+="""    -pt sample \\\n"""
    sout+="""    --omitDepthOutputAtEachBase \\\n"""
    sout+="""    -ct 10 -ct 20 -ct 30 \n"""
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
    sout+="""middfix='_HC3.3'\n"""
    sout+="""\n"""
    sout+="""echo 'BAM for HC:' $bam_ready\n"""
    sout+="""\n"""
    sout+="""echo '************** Launching GATK HaplotypeCaller ******************'\n"""
    sout+="""$GATK \\\n"""
    sout+="""    -T HaplotypeCaller \\\n"""
    sout+="""    -R $ref \\\n"""
    if param["AnalysisMode"]=="EXOME" or param["AnalysisMode"]=="PANEL":
        sout+="""    -L $targets \\\n"""
    sout+="""    -I $bam_ready \\\n"""
    sout+="""    -o ${resvcfdir}$( echo $bam_ready | sed -e 's#'$bamdir'##; s#.bam##' )${middfix}".gvcf" \\\n"""
    sout+="""    --dbsnp $vcfdir'/dbsnp_138.hg19.vcf.gz' \\\n"""
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
    sout+="""\n"""
    sout+="""middfix='_HC3.3'\n"""
    sout+="""\n"""
    sout+="""cd $resvcfdir\n"""
    sout+="""prefix='all'\n"""
    sout+="""\n"""
    sout+="""gvcflist=$resdir"/gvcfs"$middfix".list"\n"""
    sout+="""\ls $resvcfdir/*.gvcf > $gvcflist\n"""
    sout+="""\n"""
    sout+="""echo "************* Launching GATK GenotypeGVCFs ******************"\n"""
    sout+="""$GATK \\\n"""
    sout+="""    -T GenotypeGVCFs \\\n"""
    sout+="""    -V $gvcflist \\\n"""
    if param["AnalysisMode"]=="EXOME" or param["AnalysisMode"]=="PANEL":
        sout+="""    -L $targets\\\n"""
    sout+="""    -o $resvcfdir/$prefix$middfix.vcf \\\n"""
    sout+="""    -R $ref \\\n"""
    sout+="""    -nt $nthreads \\\n"""
    sout+="""    --dbsnp $vcfdir"/dbsnp_138.hg19.vcf.gz" \n"""
    if param["AnalysisMode"]=="PANEL":
        sout+="""    -A MappingQualityRankSumTest \\\n"""
        sout+="""    -A ReadPosRankSumTest \\\n"""
    sout+="""\n"""
    if param.has_key("noPhoneHome"):
        sout=noPhoneHome(param,sout)
    if param["AnalysisMode"]=="EXOME" or param["AnalysisMode"]=="GENOME":
        sout+="""echo "************** Launching GATK VariantRecalibrator -- SNP pass ******************"\n"""
        sout+="""$GATK \\\n"""
        sout+="""    -T VariantRecalibrator \\\n"""
        sout+="""    -nt $nthreads \\\n"""
        sout+="""    -R $ref \\\n"""
        sout+="""    -input $resvcfdir/$prefix$middfix.vcf \\\n"""
        if param["AnalysisMode"]=="EXOME":
            sout+="""    -L $targets \\\n"""
        sout+="""    -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 $vcfdir"/hapmap_3.3.hg19.sites.vcf.gz" \\\n"""
        sout+="""    -resource:omni,VCF,known=false,training=true,truth=true,prior=12.0 $vcfdir"/1000G_omni2.5.hg19.sites.vcf.gz" \\\n"""
        sout+="""    -resource:1000G,known=false,training=true,truth=false,prior=10.0 $vcfdir"/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz" \\\n"""
        sout+="""    -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=2.0 $vcfdir"/dbsnp_138.hg19.vcf.gz" \\\n"""
        sout+="""    -an QD -an MQRankSum -an ReadPosRankSum -an FS \\\n"""
        sout+="""    -mode SNP \\\n"""
        sout+="""    -recalFile $resvcfdir/$prefix$middfix.recal_snps \\\n"""
        sout+="""    -tranchesFile $resvcfdir/$prefix$middfix.tranches_snps \\\n"""
        sout+="""    -rscriptFile $resvcfdir/$prefix$middfix"_VQSR_SNP.R" \\\n"""
        sout+="""    --maxGaussians 4\n"""
        sout+="""\n"""
        if param.has_key("noPhoneHome"):
            sout=noPhoneHome(param,sout)
        sout+="""echo "************** Launching GATK ApplyRecalibration -- SNP pass ******************"\n"""
        sout+="""$GATK \\\n"""
        sout+="""    -T ApplyRecalibration \\\n"""
        sout+="""    -nt $nthreads \\\n"""
        sout+="""    -R $ref \\\n"""
        sout+="""    -input $resvcfdir/$prefix$middfix.vcf \\\n"""
        if param["AnalysisMode"]=="EXOME" or param["AnalysisMode"]=="PANEL":
            sout+="""    -L $targets \\\n"""
        sout+="""    --ts_filter_level 99.0 \\\n"""
        sout+="""    -recalFile $resvcfdir/$prefix$middfix.recal_snps \\\n"""
        sout+="""    -tranchesFile $resvcfdir/$prefix$middfix.tranches_snps \\\n"""
        sout+="""    -mode SNP \\\n"""
        sout+="""    -o $resvcfdir/$prefix$middfix"_recal_step1.vcf" \n"""
        sout+="""\n"""
        if param.has_key("noPhoneHome"):
            sout=noPhoneHome(param,sout)
        sout+="""echo "************** Launching GATK VariantRecalibrator -- INDEL pass ******************"\n"""
        sout+="""$GATK \\\n"""
        sout+="""    -T VariantRecalibrator \\\n"""
        sout+="""    -nt $nthreads \\\n"""
        sout+="""    -R $ref \\\n"""
        sout+="""    -input $resvcfdir/$prefix$middfix"_recal_step1.vcf" \\\n"""
        if param["AnalysisMode"]=="EXOME":
            sout+="""    -L $targets \\\n"""
        sout+="""    -resource:mills,VCF,known=false,training=true,truth=true,prior=12.0 $vcfdir"/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz" \\\n"""
        sout+="""    -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=2.0 $vcfdir"/dbsnp_138.hg19.vcf.gz" \\\n"""
        sout+="""    -an FS -an ReadPosRankSum -an MQRankSum \\\n"""
        sout+="""    -mode INDEL \\\n"""
        sout+="""    -recalFile $resvcfdir/$prefix$middfix.recal_indels \\\n"""
        sout+="""    -tranchesFile $resvcfdir/$prefix$middfix.tranches_indels \\\n"""
        sout+="""    -rscriptFile $resvcfdir/$prefix$middfix"_VQSR_INDEL.R" \\\n"""
        sout+="""    --maxGaussians 4\n"""
        sout+="""\n"""
        if param.has_key("noPhoneHome"):
            sout=noPhoneHome(param,sout)
        sout+="""echo "************** Launching GATK ApplyRecalibration -- INDEL pass ******************"\n"""
        sout+="""$GATK \\\n"""
        sout+="""    -T ApplyRecalibration \\\n"""
        sout+="""    -nt $nthreads \\\n"""
        sout+="""    -R $ref \\\n"""
        sout+="""    -input $resvcfdir/$prefix$middfix"_recal_step1.vcf" \\\n"""
        if param["AnalysisMode"]=="EXOME":
            sout+="""    -L $targets \\\n"""
        sout+="""    --ts_filter_level 99.0 \\\n"""
        sout+="""    -tranchesFile $resvcfdir/$prefix$middfix.tranches_indels \\\n"""
        sout+="""    -recalFile $resvcfdir/$prefix$middfix.recal_indels \\\n"""
        sout+="""    -mode INDEL \\\n"""
        sout+="""    -o $resvcfdir/$prefix$middfix"_recal_final.vcf" \n"""
        sout+="""\n"""
        if param.has_key("noPhoneHome"):
            sout=noPhoneHome(param,sout)
    elif param["AnalysisMode"]=="PANEL":
        sout+="""echo "************* Extract SNPs ******************"\n"""
        sout+="""$GATK \\\n"""
        sout+="""    -T SelectVariants \\\n"""
        sout+="""    -R $ref \\\n"""
        sout+="""    -V $resvcfdir/$prefix$middfix.vcf \\\n"""
        sout+="""    -L $targets \\\n"""
        sout+="""    -selectType SNP \\\n"""
        sout+="""    -nt $nthreads \\\n"""
        sout+="""    -o $resvcfdir/$prefix$middfix"_snps.vcf" \n"""
        sout+="""\n"""
        sout+="""echo "************* Hard-filter SNPs ******************"\n"""
        sout+="""$GATK \\\n"""
        sout+="""    -T VariantFiltration \\\n"""
        sout+="""    -R $ref \\\n"""
        sout+="""    -V $resvcfdir/$prefix$middfix"_snps.vcf" \\\n"""
        sout+="""    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0" \\\n"""
        sout+="""    --filterName "my_snp_filter" \\\n"""
        sout+="""    -nt $nthreads \\\n"""
        sout+="""    -o $resvcfdir/$prefix$middfix"_snps_HardFiltered.vcf" \n"""
        sout+="""\n"""
        sout+="""echo "************* Extract INDELs ******************"\n"""
        sout+="""$GATK \\\n"""
        sout+="""    -T SelectVariants \\\n"""
        sout+="""    -R $ref \\\n"""
        sout+="""    -V $resvcfdir/$prefix$middfix.vcf \\\n"""
        sout+="""    -L $targets \\\n"""
        sout+="""    -selectType INDEL \\\n"""
        sout+="""    -nt $nthreads \\\n"""
        sout+="""    -o $resvcfdir/$prefix$middfix"_indels.vcf" \n"""
        sout+="""\n"""
        sout+="""echo "************* Hard-filter INDELs ******************"\n"""
        sout+="""$GATK \\\n"""
        sout+="""    -T VariantFiltration \\\n"""
        sout+="""    -R $ref \\\n"""
        sout+="""    -V $resvcfdir/$prefix$middfix"_indels.vcf" \\\n"""
        sout+="""    --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \\\n"""
        sout+="""    --filterName "my_indel_filter" \\\n"""
        sout+="""    -nt $nthreads \\\n"""
        sout+="""    -o $resvcfdir/$prefix$middfix"_indels_HardFiltered.vcf" \n"""
        sout+="""\n"""
    if param["AnalysisMode"]!="PANEL":
        sout+="""echo "************** Launching Variant Effect Predictor (VEP - Ensembl) annotation ******************"\n"""
        sout+="""$VEP \\\n"""
        sout+="""    -i $resvcfdir/$prefix$middfix"_recal_final.vcf" \\\n"""
        sout+="""    --cache \\\n"""
        sout+="""    --everything \\\n"""
        sout+="""    --vcf \\\n"""
        sout+="""    --force_overwrite \\\n"""
        sout+="""    --dir_cache $GLOBALSCRATCH/genomes/homo_sapiens/hg19/annotations/VEP \\\n"""
        sout+="""    --dir_plugins $GLOBALSCRATCH/genomes/homo_sapiens/hg19/annotations/VEP/Plugins \\\n"""
        sout+="""    --plugin LoF \\\n"""
        sout+="""    --offline \\\n"""
        sout+="""    --output_file $resvcfdir/$prefix$middfix"_recal_final_VEP.vcf" \\\n"""
        sout+="""    --stats_file $resvcfdir/$prefix$middfix"_recal_final_VEP_summary.html" \n"""
        sout+="""\n"""
        sout+="""randomNumber=$RANDOM \n"""
        sout+="""$VCFSORTER $dict $resvcfdir/$prefix$middfix"_recal_final_VEP.vcf" > swap_$randomNumber \n"""
        sout+="""mv swap_$randomNumber $resvcfdir/$prefix$middfix"_recal_final_VEP.vcf" \n"""
        sout+="""\n"""
        # sout+="""echo "************** Prepare VCF (remove 'chr') ******************"\n"""
        # sout+="""vcf=$resvcfdir/$prefix$middfix"_recal_final_VEP.vcf"\n"""
        # sout+="""randomNumber=$RANDOM\n"""
        # sout+="""sed 's/^chr//' $vcf | sed 's/ID=chr/ID=/' > ${vcf/.vcf/_noCHR.vcf}\n"""
        # sout+="""$BGZIP -f ${vcf/.vcf/_noCHR.vcf}\n"""
        # sout+="""$TABIX -f -p vcf ${vcf/.vcf/_noCHR.vcf.gz}\n"""
        # sout+="""\n"""
        # sout+="""echo "************** Launching Allelic Frequencies Annotation -- 1KG ******************"\n"""
        # sout+="""$SNPSIFT \\\n"""
        # sout+="""    annotateMem \\\n"""
        # sout+="""    -info EUR_AF \\\n"""
        # sout+="""    -v $kgVCF \\\n"""
        # sout+="""    ${vcf/.vcf/_noCHR.vcf.gz} | $BGZIP -f > ${vcf/.vcf/_noCHR_1KG.vcf.gz}\n"""
        # sout+="""\n"""
        # sout+="""echo "************** Launching Allelic Frequencies Annotation -- ExAC ******************"\n"""
        # sout+="""$SNPSIFT \\\n"""
        # sout+="""    annotateMem \\\n"""
        # sout+="""    -info ExAC_AF \\\n"""
        # sout+="""    -v $exacVCF \\\n"""
        # sout+="""    ${vcf/.vcf/_noCHR.vcf_1KG.gz} | $BGZIP -f > ${vcf/.vcf/_noCHR_1KG_ExAC.vcf.gz}\n"""
        # sout+="""\n"""
        # sout+="""echo "************** Launching ClinVar Annotation ******************"\n"""
        # sout+="""$SNPSIFT \\\n"""
        # sout+="""    annotateMem \\\n"""
        # sout+="""    -v $clinvarVCF \\\n"""
        # sout+="""    ${vcf/.vcf/_noCHR_1KG_ExAC.vcf.gz} | $BGZIP -f > ${vcf/.vcf/_noCHR_1KG_ExAC_ClinVar.vcf.gz}\n"""
        # sout+="""\n"""
        sout+="""echo "**************** Preparing for SnpSift ***************************"\n"""
        sout+="""$BGZIP $resvcfdir/$prefix$middfix"_recal_final_VEP.vcf"\n"""
        sout+="""$TABIX -p vcf $resvcfdir/$prefix$middfix"_recal_final_VEP.vcf.gz"\n"""
        sout+="""\n"""
        sout+="""VCFGZ=$resvcfdir/$prefix$middfix"_recal_final_VEP.vcf.gz"\n"""
        sout+="""\n"""
        sout+="""echo "**************** Annotate w/ 1KG frequencies *********************"\n"""
	sout+="""randomNumber=$RANDOM\n"""
	sout+="""$SNPSIFT \\\n"""
	sout+="""    annotateMem \\\n"""
	sout+="""    -info EUR_AF \\\n"""
	sout+="""    -v $GLOBALSCRATCH/genomes/homo_sapiens/hg19/allelic_frequencies/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz \\\n"""
	sout+="""    $VCFGZ | $BGZIP > swap_$randomNumber\n"""
	sout+="""if (( `wc -l  swap_$randomNumber | awk '{print $1}'` > 0)) ; then\n"""
	sout+="""    mv swap_$randomNumber $VCFGZ\n"""
	sout+="""fi\n"""
	sout+="""\n"""
	sout+="""echo "*************** Annotate w/ ExAC frequencies *********************"\n"""
	sout+="""randomNumber=$RANDOM\n"""
	sout+="""$SNPSIFT \\\n"""
	sout+="""    annotateMem \\\n"""
	sout+="""    -info ExAC_AF \\\n"""
	sout+="""    -v $GLOBALSCRATCH/genomes/homo_sapiens/hg19/allelic_frequencies/ExAC.r0.3.sites.vep.ExAC_AF.vcf.gz \\\n"""
	sout+="""    $VCFGZ | $BGZIP > swap_$randomNumber\n"""
	sout+="""if (( `wc -l  swap_$randomNumber | awk '{print $1}'` > 0)) ; then\n"""
	sout+="""    mv swap_$randomNumber $VCFGZ\n"""
	sout+="""fi\n"""
	sout+="""\n"""
	sout+="""echo "************** Annotate w/ ClinVar *******************************"\n"""
	sout+="""randomNumber=$RANDOM\n"""
	sout+="""$SNPSIFT \\\n"""
	sout+="""    annotateMem \\\n"""
	sout+="""    -v $GLOBALSCRATCH/genomes/homo_sapiens/hg19/allelic_frequencies/clinvar_20141202.vcf.gz \\\n"""
	sout+="""\\\n"""
	sout+="""    $VCFGZ | $BGZIP > swap_$randomNumber\n"""
	sout+="""if (( `wc -l  swap_$randomNumber | awk '{print $1}'` > 0)) ; then\n"""
	sout+="""    mv swap_$randomNumber $VCFGZ\n"""
	sout+="""fi\n"""
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
