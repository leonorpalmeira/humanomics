import scala.sys.process._ // for system operations
import java.io.File

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.extensions.picard._
import org.broadinstitute.gatk.queue.function._
import org.broadinstitute.gatk.utils.commandline._

class HumanomicsPipeline extends QScript {
  qscript =>
    
  /****************************************************************************
   * Required Parameters
   *****************************************************************************/

  // NOTE TO SELF: _ is shorthand for null

  // Input/Output/Resources
  @Input(doc="Fasta.gz file to use as query 1 (paired-end)", shortName="I1", fullName="input_file_1")
  var query1: File = _
  @Input(doc="Fasta.gz file to use as query 2 (paired-end)", shortName="I2", fullName="input_file_2")
  var query2: File = _
  @Argument(doc="Global output directory", shortName="O", fullName="output_dir")
  var outdir: File = _
  @Argument(doc="Read Group string (@RG) for this sample.", shortName="RG", fullName="read_group")
  var RG: String = _
  @Argument(doc="The reference BWA-indexed files.", shortName="Rbwa", fullName="reference_bwa")
  var referenceBWA: String = _
  @Input(doc="The reference FASTA file for the GATK tools.", shortName="R", fullName="reference_sequence")
  var referenceFasta: File = _
  @Input(doc="Known indels", shortName="indels", fullName="knownIndels")
  var knownIndels: Seq[File] = Nil
  @Input(doc="Known SNPs", shortName="snps", fullName="knownSNPs")
  var knownSNPs: File = _ //Seq[File] = Nil

  // Binaries
  @Argument(doc="Java invocation to use", shortName="javacustom")
  var javacustom: String = _
  @Argument(doc="bwa binary to use", shortName="bwa")
  var bwa: String = _
  @Argument(doc="samtools binary to use", shortName="samtools")
  var samtools: String = _
  @Argument(doc="samstat binary to use", shortName="samstat")
  var samstat: String = _
  @Argument(doc="fastqc binary to use", shortName="fastqc")
  var fastqc: String = _
  @Argument(doc="picard java invocation to use", shortName="picard")
  var picard: String = _
  
  /****************************************************************************
   * Optional Parameters
   *****************************************************************************/

  @Argument(doc="How many chunks of scatter/gather", shortName="P", fullName="scatter_gather", required=false)
  var scatter:  Int = -1
  @Argument(doc="Number of data threads to allocate to this analysis (GATK-specific: data multithreading)", shortName="N", fullName="num_threads", required=false)
  var nt:  Int = 1
  @Argument(doc="Number of CPU threads to allocate per data thread (GATK and non-GATK: CPU multithreading)", shortName="C", fullName="num_cpu_threads_per_data_thread", required=false)
  var nct:  Int = 1
  @Input(doc="Intervals to realign", shortName="L", fullName="intervals", required = false)
  var intervals: File = _
  @Input(doc="dbSNP file", shortName="dbsnp", fullName="knownDbsnp", required=false)
  var dbsnp: File = _
//  @Input(doc="Known in/del VCF file", shortName="indel", fullName="knownIndels", required=false)
//  var known: File = _

  // noET option
  @Argument(doc="Run reporting mode", shortName="et", fullName="phone_home", required=false)
  var ET: org.broadinstitute.gatk.engine.phonehome.GATKRunReport.PhoneHomeOption = _
  @Argument(doc="GATK key file required to run with -et NO_ET (NoPhoneHome option)", shortName="K", fullName="gatk_key", required=false)
  var ETkey: File = _
  
  /****************************************************************************
  * Global Variables
  ****************************************************************************/

  val queueLogDir: String = outdir + "/.qlog/"  // Gracefully hide Queue's output

  /****************************************************************************
  * Helper classes and methods
  ****************************************************************************/

  case class FastQC (inFastqz: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="Fastq.gz file (will not be extracted)") var query = inFastqz
    def commandLine = fastqc + " --noextract --outdir=" + outdir +"/fastqc " + query
    this.analysisName = queueLogDir + query + ".fastqc"
    this.jobName = queueLogDir + query + ".fastqc"
  }

  case class BwaMem (inFastq1: File, inFastq2: File, inReference: File, outBam: File, RG: String) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="Fastq (1) of a paired-end") var query1 = inFastq1
    @Input(doc="Fastq (2) of a paired-end") var query2 = inFastq2
    @Input(doc="Reference database") var reference = inReference
    @Output(doc="Bam output file") var bam = outBam
    def commandLine = bwa + " mem -M " + optional(" -t ", nct) + optional(" -R ", RG) + required(" ", reference) + required(" ", query1) + required(" ", query2) + " | " + samtools + " view -bSho " + required(bam) + " - "
    this.analysisName = queueLogDir + outBam + ".bwamem"
    this.jobName = queueLogDir + outBam + ".bwamem"
  }

  case class SamtoolsSort (inBam: File, outBam: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="Bam to be sorted") var inbam = inBam
    @Output(doc="Sorted bam") var outbam = outBam
    def commandLine = samtools + " sort " + inbam + " " + swapExt(outbam, ".bam", " " ) // Removes .bam suffix, which is appended by samtools sort...
    this.analysisName = queueLogDir + outBam + ".samsort"
    this.jobName = queueLogDir + outBam + ".samsort"
  }

  case class Samstat (inBam: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="Bam to be analysed") var inbam = inBam
//    @Output(doc="Html file") var outhtml=outHtml //$bam_sorted".html"
    def commandLine = samstat + " " + inbam + " && " + " mv " + inbam + ".html " + outdir + "/samstat/." // mv bamsorted to samstatdir
    this.analysisName = queueLogDir + inBam + ".samstat"
    this.jobName = queueLogDir + inBam + ".samstat"
  }

  case class DeDuplicate (inBam: File, outBam: File, metricsFile: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="Bam to be de-duplicated") var inbam = inBam
    @Output(doc="De-duplicated bam") var outbam = outBam
    @Output(doc="Metrics file") var metrics = metricsFile
    def commandLine = javacustom + " " + picard + "/MarkDuplicates.jar INPUT=" + inBam + " OUTPUT=" + outBam + " METRICS_FILE=" + metrics + " REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=/scratch/ulg/genan/palmeira/tmp"
    this.analysisName = queueLogDir + outBam + ".rmd"
    this.jobName = queueLogDir + outBam + ".rmd"
  }

  case class SamIndex (inBam: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="Bam to be indexed") var inbam = inBam
    def commandLine = samtools + " index " + inbam
    this.analysisName = queueLogDir + inBam + ".index"
    this.jobName = queueLogDir + inBam + ".index"
  }

  // /****************************************************************************
  //  * CommonArguments
  //  ******************************************************************************/

  // General arguments to non-GATK tools
  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = 4
    // NOTE TO SELF: this will erase all files made during the pipeline once the pipeline has finished successfully. Circumvent by -keepIntermediates when calling the Queue script
    this.isIntermediate = true
  }
  
  // General arguments to GATK walkers
  trait GATKCommonArgs extends CommandLineGATK {
    this.reference_sequence = qscript.referenceFasta
    this.intervals = if (qscript.intervals == null) Nil else List(qscript.intervals)
    // NOTE TO SELF: Nil is equivalent to List()
    // Set the memory limit to 2 gigabytes on each command.
    this.memoryLimit = 20
    // added for NO_ET:
    this.phone_home = qscript.ET
    this.gatk_key = qscript.ETkey
  }

  /****************************************************************************  
   * Main script
   *****************************************************************************/

  def script() {
    
    // the following should also catch errors when empty or no fastq.gz:
//    val queries = new java.io.File(querydir).listFiles.filter(_.getName.endsWith(".fastq.gz"))
//    val query1 = queries(0) // get first value -> TO BE MODIFIED for more files
//    val query2 = queries(1) // get second value -> TO BE MODIFIED for more files

    val tmp = query1.split("/")
    val bam = swapExt(outdir + "/bam/" + tmp(tmp.length-1), ".fastq.gz", "_RG.bam")

    val sorted_bam = swapExt(bam, "_RG.bam", "_sorted_RG.bam") 
    val dedup_bam = swapExt(bam, "_RG.bam", "_sorted_RG_rmd.bam")
    val dedup_metrics = swapExt(bam, "_RG.bam", "_sorted_RG_rmd.rmd")

    // defining all pipelining steps:
    //-------------------------------
    val fastqc1 = new FastQC(qscript.query1)
    val fastqc2 = new FastQC(qscript.query2)
    val map = new BwaMem(qscript.query1, qscript.query2, qscript.referenceBWA, bam, qscript.RG)
    val samsort = new SamtoolsSort(bam, sorted_bam)
    val samstat = new Samstat(sorted_bam)
    val dedup = new DeDuplicate(sorted_bam, dedup_bam, dedup_metrics)
    val samindex = new SamIndex(dedup_bam)
    val targetCreator = new RealignerTargetCreator with GATKCommonArgs
    val indelRealigner = new IndelRealigner with GATKCommonArgs
    val bqsr = new BaseRecalibrator with GATKCommonArgs
    val bqsr2 = new BaseRecalibrator with GATKCommonArgs
    val analyzeCovariates = new AnalyzeCovariates with GATKCommonArgs
    val applyRecalibration = new PrintReads with GATKCommonArgs
    val hc = new HaplotypeCaller with GATKCommonArgs
 
    targetCreator.input_file +:= dedup_bam
    targetCreator.out = swapExt(bam, "_RG.bam", "_sorted_RG_rmd.intervals")
    targetCreator.nt = qscript.nt
    targetCreator.known = qscript.knownIndels
    targetCreator.scatterCount = qscript.scatter
    targetCreator.isIntermediate = true

    indelRealigner.input_file +:= dedup_bam
    indelRealigner.targetIntervals = targetCreator.out
    indelRealigner.out = swapExt(bam, "_RG.bam", "_sorted_RG_rmd_indelreal.bam")
    indelRealigner.known = qscript.knownIndels
    indelRealigner.scatterCount = qscript.scatter
    indelRealigner.isIntermediate = true

    bqsr.input_file +:= indelRealigner.out
    bqsr.out = swapExt(bam, "_RG.bam", "_sorted_RG_rmd_indelreal_recal1.grp")
    bqsr.knownSites = qscript.knownIndels
    bqsr.knownSites +:= qscript.knownSNPs
    bqsr.nct = qscript.nct
    bqsr.scatterCount = qscript.scatter
    bqsr.isIntermediate = true

    bqsr2.input_file +:= indelRealigner.out
    bqsr2.out = swapExt(bqsr.out, "_recal1.grp", "_recal2.grp")
    bqsr2.knownSites = qscript.knownIndels
    bqsr2.knownSites +:= qscript.knownSNPs
    bqsr2.BQSR = bqsr.out
    bqsr2.nct = qscript.nct
    bqsr2.scatterCount = qscript.scatter
    bqsr2.isIntermediate = true

    analyzeCovariates.before = bqsr.out
    analyzeCovariates.after = bqsr2.out
    analyzeCovariates.csv = swapExt(bqsr.out, "_recal1.grp", "_bqsr.csv")
    analyzeCovariates.plots = swapExt(bqsr.out, "_recal1.grp", "_bqsr.pdf")
    analyzeCovariates.scatterCount = qscript.scatter

    applyRecalibration.input_file +:= indelRealigner.out
    applyRecalibration.BQSR = bqsr2.out
    applyRecalibration.out = swapExt(bam, "_RG.bam", "_sorted_RG_rmd_indelreal_BQSR.bam")
    applyRecalibration.nct = qscript.nct
    applyRecalibration.scatterCount = qscript.scatter

    hc.input_file +:= applyRecalibration.out
    hc.nct = qscript.nct
    hc.pairHMM = org.broadinstitute.gatk.utils.pairhmm.PairHMM.HMM_IMPLEMENTATION.VECTOR_LOGLESS_CACHING
    hc.pcr_indel_model = org.broadinstitute.gatk.tools.walkers.haplotypecaller.PairHMMLikelihoodCalculationEngine.PCR_ERROR_MODEL.CONSERVATIVE
    hc.out = swapExt(bam, "_RG.bam", "_sorted_RG_rmd_indelreal_BQSR.g.vcf") // NOTE TO SELF: has to have the .vcf suffix for the Queue aggregation to work
    hc.emitRefConfidence = org.broadinstitute.gatk.tools.walkers.haplotypecaller.ReferenceConfidenceMode.GVCF
    hc.variant_index_type = org.broadinstitute.gatk.utils.variant.GATKVCFIndexType.LINEAR
    hc.variant_index_parameter = 128000
    hc.scatterCount = qscript.scatter

    add(fastqc1, fastqc2, map, samsort, samstat, dedup, samindex, targetCreator, indelRealigner, bqsr, bqsr2, analyzeCovariates, applyRecalibration, hc)

  }
}

