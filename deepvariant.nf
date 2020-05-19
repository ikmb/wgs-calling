#!/usr/bin/env nextflow

folder = params.folder
run_mode = "WGS"

params.outdir = "deepvariant"

OUTDIR = file(params.outdir)

// Specifies the underlying genome assembly
params.assembly = "hg38"

// *****************************************
// Assembly-specific variables and resources
// *****************************************

if (params.genomes.containsKey(params.assembly) == false) {
   exit 1, "Specified unknown genome assembly, please consult the documentation for valid assemblies."
}

REF = file(params.genomes[ params.assembly ].fasta)
INTERVAL_LIST = params.interval_list ? file(params.interval_list) : file(params.genomes[params.assembly].interval_list)

// ******************
// Misc
// ******************

params.email = false

// Whether to use a local scratch disc
use_scratch = params.scratch

VERSION = "0.1"

// Header log info
log.info "========================================="
log.info "GATK Best Practice for Genome-Seq calling v${VERSION}"
log.info "Section:             		DeepVariant"
log.info "Nextflow Version:		$workflow.nextflow.version"
log.info "Assembly version:		${params.assembly}"
log.info "Command Line:			$workflow.commandLine"
log.info "========================================="

(fastaToIndexCh, fastaToGzCh, fastaToGzFaiCh, fastaToGziCh) = Channel.fromPath(REF).into(4)

if(params.fai){
faiToExamples = Channel
    .fromPath(params.fai)
    .ifEmpty{exit 1, "Fai file not found: ${params.fai}"}
}

if(params.fastagz){
fastaGz = Channel
    .fromPath(params.fastagz)
    .ifEmpty{exit 1, "Fastagz file not found: ${params.fastagz}"}
    .into {fastaGzToExamples; fastaGzToVariants }
}

if(params.gzfai){
gzFai = Channel
    .fromPath(params.gzfai)
    .ifEmpty{exit 1, "gzfai file not found: ${params.gzfai}"}
    .into{gzFaiToExamples; gzFaiToVariants }
}

if(params.gzi){
gzi = Channel
    .fromPath(params.gzi)
    .ifEmpty{exit 1, "gzi file not found: ${params.gzi}"}
    .into {gziToExamples; gziToVariants}
}

if(params.bed) {
	bedToExamples = Channel
	.fromPath(params.bed)
	.ifEmpty {exit 1, "bed file not found: ${params.bed}"}
}

Channel.fromPath("${params.folder}/*.*am")
	.set { inputDv }

// ****************************
// START WORKFLOW - MAKE INDICES
// ****************************

if(!params.fai) {
  process preprocess_fai {
      tag "${fasta}.fai"
      //publishDir "$baseDir/sampleDerivatives"

      input:
      file(fasta) from fastaToIndexCh

      output:
      file("${fasta}.fai") into faiToExamples

      script:
      """
      samtools faidx $fasta
      """
  }
}

if(!params.fastagz) {
  process preprocess_fastagz {
      tag "${fasta}.gz"
      //publishDir "$baseDir/sampleDerivatives"

      input:
      file(fasta) from fastaToGzCh

      output:
      file("*.gz") into (tmpFastaGzCh, fastaGzToExamples, fastaGzToVariants)

      script:
      """
      bgzip -c ${fasta} > ${fasta}.gz
      """
  }
}

if(!params.gzfai) {
  process preprocess_gzfai {
    tag "${fasta}.gz.fai"
    //publishDir "$baseDir/sampleDerivatives"

    input:
    file(fasta) from fastaToGzFaiCh
    file(fastagz) from tmpFastaGzCh

    output:
    file("*.gz.fai") into (gzFaiToExamples, gzFaiToVariants)

    script:
    """
    samtools faidx $fastagz
    """
  }
}

if(!params.gzi){
  process preprocess_gzi {
    tag "${fasta}.gz.gzi"
    //publishDir "$baseDir/sampleDerivatives"

    input:
    file(fasta) from fastaToGziCh

    output:
    file("*.gz.gzi") into (gziToExamples, gziToVariants)

    script:
    """
    bgzip -c -i ${fasta} > ${fasta}.gz
    """
  }
}

if(!params.bed) {
  process runMakeBed {
	input:
	file(intervals) from Intervals
	
	output:
	file(bed_file) into bedToExamples

	script:

	bed_file = intervals.getBaseName() + ".bed"

	"""
		picard IntervalListToBed I=$INTERVAL_LIST O=$bed_file
	"""
  }

}
// **************************
// RUN DEEPVARIANT
// **************************

process runDV {

  publishDir "${OUTDIR}/variants", mode: 'copy'

  scratch true

  input:
  file(bam) from inputDv
  file fai from faiToExamples.collect()
  file fastagz from fastaGzToExamples.collect()
  file gzfai from gzFaiToExamples.collect()
  file gzi from gziToExamples.collect()
  file bed from bedToExamples.collect()

  output:
  set file(vcf),file(tbi) into DvVCF

  script:

  vcf = bam.getBaseName() + ".vcf.gz"
  gvcf = bam.getBaseName() + ".g.vcf.gz" 
  vcf_tbi = vcf + ".tbi"
  gvcf_tbi = gvcf + ".tbi"

  """
	run_deepvariant \
	--model_type=$run_mode \
	--ref=$fastagz \
  	--reads=$bam \
  	--regions $bed \
  	--output_vcf=$vcf \
  	--output_gvcf=$gvcf \
  	--num_shards=${task.cpus} 	
  """

}

if (params.email) {
	workflow.onComplete {
	    def subject = 'WGS DeepVariant step finished.'
	    def recipient = params.email

	    ['mail', '-s', subject, recipient].execute() << """

	    Pipeline execution summary
	    ---------------------------
	    Completed at: ${workflow.complete}
	    Duration    : ${workflow.duration}
	    Success     : ${workflow.success}
	    workDir     : ${workflow.workDir}
	    exit status : ${workflow.exitStatus}
	    Error report: ${workflow.errorReport ?: '-'}
	    """
	}
}



