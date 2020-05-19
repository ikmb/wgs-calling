#!/usr/bin/env nextflow

folder = params.folder
model = "WGS"

params.outdir = "deepvariant"

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

bams = Channel.fromPath(params.folder + "/*am") 
/**
**/
bais = Channel.fromPath(params.folder + "/*ai")
/**
**/

inputDv = bams.merge(bais)

// ****************************
 START WORKFLOW - MAKE INDICES
// ****************************/

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
		gatk IntervalListToBed -I $INTERVAL_LIST -O $bed_file
	"""
  }

}
// **************************
// RUN DEEPVARIANT
// **************************

/********************************************************************
  process make_examples
  Getting bam files and converting them to images ( named examples )
********************************************************************/

process make_examples{

  tag "${bam}"
  publishDir "${params.outdir}/make_examples", mode: 'copy',
  saveAs: {filename -> "logs/log"}

  input:
  file fai from faiToExamples.collect()
  file fastagz from fastaGzToExamples.collect()
  file gzfai from gzFaiToExamples.collect()
  file gzi from gziToExamples.collect()
  file bed from bedToExamples.collect()
  set file(bam), file(bai) from inputDv

  output:
  set file("${bam}"),file(bai),file('*_shardedExamples') into examples

  script:
  """
  mkdir logs
  mkdir ${bam.baseName}_shardedExamples
  dv_make_examples.py \
  --cores ${task.cpus} \
  --sample ${bam} \
  --ref ${fastagz} \
  --reads ${bam} \
  --regions ${bed} \
  --logdir logs \
  --examples ${bam.baseName}_shardedExamples
  """
}
/********************************************************************
  process call_variants
  Doing the variant calling based on the ML trained model.
********************************************************************/

process call_variants{

  tag "${bam}"

  input:
  set file(bam),file(bai),file(shardedExamples) from examples

  output:
  set file(bam),file(bai),file('*_call_variants_output.tfrecord') into called_variants

  script:
  """
  dv_call_variants.py \
    --cores ${task.cpus} \
    --sample ${bam} \
    --outfile ${bam.baseName}_call_variants_output.tfrecord \
    --examples $shardedExamples \
    --model ${model}
  """
}

/********************************************************************
  process postprocess_variants
  Trasforming the variant calling output (tfrecord file) into a standard vcf file.
********************************************************************/

process postprocess_variants{

  tag "${bam}"

  publishDir params.outdir, mode: 'copy'

  input:
  file fastagz from fastaGzToVariants.collect()
  file gzfai from gzFaiToVariants.collect()
  file gzi from gziToVariants.collect()
  set file(bam),file(bai),file('call_variants_output.tfrecord') from called_variants

  output:
  set val("${vcf_gz}"),file("${vcf_gz_tbi}") into postout

  script:
  vcf = bam.getBaseName() + ".vcf"
  vcf_gz = vcf + ".gz"
  vcf_gz_tbi = vcf +".tbi"

  """
  dv_postprocess_variants.py \
  --ref ${fastagz} \
  --infile call_variants_output.tfrecord \
  --outfile "${bam}.vcf"

  bgzip $vcf && tabix $vcfgz
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



