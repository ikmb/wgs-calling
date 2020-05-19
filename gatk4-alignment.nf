#!/usr/bin/env nextflow

inputFile = file(params.samples)

params.outdir = "alignments"

OUTDIR = file(params.outdir)

// Specifies the underlying genome assembly
params.assembly = "hg38"

// *****************************************
// Assembly-specific variables and resources
// *****************************************

if (params.genomes.containsKey(params.assembly) == false) {
   exit 1, "Specified unknown genome assembly, please consult the documentation for valid assemblies."
}

REF = params.fasta ?: file(params.genomes[ params.assembly ].fasta)
DICT = params.dict ?: file(params.genomes[ params.assembly ].dict)
DBSNP = params.dbsnp ?: file(params.genomes[ params.assembly ].dbsnp )
MILLS = params.mills ?: file(params.genomes[ params.assembly ].mills )
if (params.intervals) {
	INTERVALS = file(params.intervals)
} else {
	INTERVALS = file(params.genomes[ params.assembly ].intervals )
}
INTERVAL_LIST = params.interval_list ?: file(params.genomes[ params.assembly ].interval_list )

// This is a bit simplistic and uses each interval , instead of pooling smaller intervals into one job. 
regions = []
// Make the region files
INTERVALS.eachLine { str ->
        if(! str.startsWith("@") ) {
                regions << str.trim()
        }
}

summary = [:]

summary['Assembly'] = params.assembly
summary['SampleList']=  params.samples
summary['Intervals'] = INTERVALS

// ******************
// Misc
// ******************

params.format = "cram"

// Whether to use a local scratch disc
use_scratch = params.scratch

// Make sure the Nextflow version is current enough

// Header log info
log.info "========================================="
log.info "GATK Best Practice for Genome-Seq Preprocessing v${workflow.manifest.version}"
log.info "Section:             		Read alignment"
log.info "Commit hash:			$workflow.commitId"
log.info "Nextflow Version:		$workflow.nextflow.version"
log.info "Assembly version:		${params.assembly}"
log.info "Command Line:			$workflow.commandLine"
log.info "========================================="

Channel.from(inputFile)
       .splitCsv(sep: ';', header: true)
       .set { inputFastp }

process runFastp {

  scratch true

  input:
  set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, center, run_date, fastqR1, fastqR2 from inputFastp

  output:
  set val(indivID), val(sampleID), val(libraryID), val(rgID), val(platform_unit), val(platform), val(platform_model), val(center), val(run_date),file(left),file(right) into inputBwa
  set val(indivID), val(sampleID), val(libraryID), file(json),file(html) into outputReportTrimming

  script:
  left = file(fastqR1).getBaseName() + ".trimmed.fastq.gz"
  right = file(fastqR2).getBaseName() + ".trimmed.fastq.gz"
  json = indivID + "_" + sampleID + "_" + libraryID + "-" + rgID + ".fastp.json"
  html = indivID + "_" + sampleID + "_" + libraryID + "-" + rgID + ".fastp.html"

  """
	fastp --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right --detect_adapter_for_pe -w ${task.cpus} -j $json -h $html
  """

}

// Run BWA on each trimmed chunk
process runBwa {

    scratch true

    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, center, run_date,file(fastqR1),file(fastqR2) from inputBwa

    output:
    set indivID, sampleID, file(outfile) into runBWAOutput

    script:
    outfile = sampleID + "_" + libraryID + "_" + rgID + ".aligned.bam"
    outfile_index = outfile + ".bai"

    """
	bwa mem -H $DICT -M -R "@RG\\tID:${rgID}\\tPL:ILLUMINA\\tPU:${platform_unit}\\tSM:${indivID}_${sampleID}\\tLB:${libraryID}\\tDS:${REF}\\tCN:${center}" -t ${task.cpus} ${REF} $fastqR1 $fastqR2 | samtools sort -m 4G -@ 4 -o $outfile - 
	samtools index $outfile
    """
}

process runFixTags {

	scratch true

	input:
    	set val(indivID), val(sampleID), file(aligned_bam) from runBWAOutput

	output:
	set indivID, sampleID, file(outfile_bam),file(outfile_bai) into runMarkDuplicatesOutput
	file(outfile_metrics) into runMarkDuplicatesOutput_QC
        file(outfile_md5) into MarkDuplicatesMD5

    	script:
    	outfile_bam = indivID + "_" + sampleID + ".dedup.bam"
	outfile_bai = indivID + "_" + sampleID + ".dedup.bai"
	outfile_md5 = indivID + "_" + sampleID + ".dedup.bam.md5"

	outfile_metrics = sampleID + "_duplicate_metrics.txt"

	"""
		gatk SetNmMdAndUqTags \
                -I $aligned_bam \
                -O alignment.clean.bam \
                -R $REF \
                --IS_BISULFITE_SEQUENCE false 
		
		gatk --java-options "-Xms4G -Xmx${task.memory.toGiga()-3}G" MarkDuplicates \
                -I alignment.clean.bam \
                -O ${outfile_bam} \
                -R ${REF} \
                -M ${outfile_metrics} \
                --CREATE_INDEX true \
                --TMP_DIR \$TMPDIR \
                --MAX_RECORDS_IN_RAM 1000000 \
                --ASSUME_SORT_ORDER coordinate \
                --CREATE_MD5_FILE true

		rm alignment.clean.bam

	"""
}

// Generate a model for base recalibration within target intervals
process runBaseRecalibrator {

    	// publishDir "${OUTDIR}/${params.assembly}/${indivID}/${sampleID}/Processing/BaseRecalibrator/batches", mode: 'copy'

    	scratch true
	    
    	input:
    	set indivID, sampleID, dedup_bam, dedup_bai from runMarkDuplicatesOutput
	each region from regions
    
    	output:
    	set indivID, sampleID, file(recal_table) into outputBaseRecalibrator
	set indivID, sampleID, dedup_bam into BamForBQSR
    
    	script:
        region_tag = region.trim().replace(/:/, '_')
    	recal_table = indivID + "." + sampleID + "." + region_tag + ".recal_table.txt" 
       
    	"""
		gatk --java-options "-Xmx${task.memory.toGiga()}G" BaseRecalibrator \
		--reference ${REF} \
		--input ${dedup_bam} \
		--known-sites ${MILLS} \
		--known-sites ${DBSNP} \
		--use-original-qualities \
		-L $region \
		-ip 150 \
		--output ${recal_table}
	"""
}

ReportsBySample = outputBaseRecalibrator.groupTuple(by: [0,1])

process runGatherBQSRReports {

        // publishDir "${OUTDIR}/${params.assembly}/${indivID}/${sampleID}/Processing/BaseRecalibrator/"

	input:
	set indivID,sampleID,file(reports) from ReportsBySample

	output:
	set indivID,sampleID,file(merged_report) into MergedReport

	script:
	sorted_reports = []
	regions.each { region ->
	        region_tag = region.trim().replace(/:/, '_')
                this_report = sampleID + "." + region_tag + ".recal_table.txt" 
                sorted_reports << reports.find { it =~ this_report }
        }	
	merged_report = indivID + "_" + sampleID + ".recal_table.txt"

	"""
		gatk --java-options "-Xmx${task.memory.toGiga()}G" GatherBQSRReports \
		-I ${sorted_reports.join(' -I ')} \
		-O $merged_report 
	"""
}

inputForApplyBQSR = MergedReport.join(BamForBQSR, by: [0,1])

process runApplyBQSR {

	publishDir "${OUTDIR}/${params.assembly}/cram/", mode: 'copy'

	scratch true
	    
	input:
	set indivID, sampleID, file(recal_table), file(realign_bam) from inputForApplyBQSR

	output:
	set indivID, sampleID, file(outfile_bam), file(outfile_bai) into (BamForWGSStats,BamForAlignStats)
	            
    	script:
    	outfile_bam = indivID + "." + sampleID + ".clean.cram"
    	outfile_bai = indivID + "." + sampleID + ".clean.cram.bai"
	outfile_md5 = indivID + "." + sampleID + ".clean.cram.md5"
          
    	"""
          gatk --java-options "-Xmx${task.memory.toGiga()}G" ApplyBQSR \
             --reference ${REF} \
             --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
             --use-original-qualities \
             --input ${realign_bam} \
             -bqsr ${recal_table} \
             --output ${outfile_bam} \
             -OBM true \
	     -OVM true
    	"""	
} 


// ------------------------------------------------------------------------------------------------------------
//
// Perform a several tasks to assess QC:
// 1) Depth of coverage over targets
// 2) Generate alignment stats, insert size stats, quality score distribution
// 3) Generate hybrid capture stats
// 4) Run FASTQC to assess read quality
//
// ------------------------------------------------------------------------------------------------------------

// ------------------------------------------------------------------------------------------------------------
//
// Plot results with multiqc
//
// ------------------------------------------------------------------------------------------------------------

process runAlignmentStats {

        publishDir "${OUTDIR}/${params.assembly}/${indivID}/${sampleID}/picard_stats", mode: 'copy'

	scratch true

	input:
	set val(indivID),val(sampleID),val(bam),val(bai) from BamForAlignStats

	output:
	file(align_stats) into AlignStats

	script:
	align_stats = indivID + "_" + sampleID + "_alignment_stats.txt"

	"""
		picard CollectAlignmentSummaryMetrics \
		R=$REF \
		I=$bam \
		O=$align_stats
	"""
	
}

process runWgsCoverage {

	publishDir "${OUTDIR}/${params.assembly}/${indivID}/${sampleID}/picard_stats", mode: 'copy'

	scratch true

	input:
	set val(indivID),val(sampleID),val(bam),val(bai) from BamForWGSStats

	output:
	file(coverage_stats) into CoverageStats

	script:
	coverage_stats = indivID + "_" + sampleID + "_wgs_stats.txt"

	"""
		picard -Xmx${task.memory.toGiga()}G CollectWgsMetrics \
		I=$bam \
		O=$coverage_stats \
		R=$REF \
		INTERVALS=$INTERVAL_LIST
	"""
}

process runMultiQCLibrary {

    tag "Generating library level summary and QC plots"
	publishDir "${OUTDIR}/Summary/Library", mode: 'copy'
	    
    input:
    file('*') from runMarkDuplicatesOutput_QC.flatten().toList()
    file('*') from AlignStats.flatten().toList()
    file('*') from outputReportTrimming.flatten().toList()
    file('*') from CoverageStats.flatten().toList()

    output:
    file("library_multiqc*") into multiqc_report

    script:

    """
    cp $baseDir/config/multiqc_config.yaml multiqc_config.yaml
    multiqc -n library_multiqc *
    """
}


workflow.onComplete {
	log.info "========================================="
	log.info "Duration:		$workflow.duration"
	log.info "========================================="

	def email_fields = [:]
	email_fields['version'] = workflow.manifest.version
	email_fields['session'] = workflow.sessionId
	email_fields['runName'] = run_name
	email_fields['Samples'] = params.samples
	email_fields['success'] = workflow.success
	email_fields['dateStarted'] = workflow.start
	email_fields['dateComplete'] = workflow.complete
	email_fields['duration'] = workflow.duration
	email_fields['exitStatus'] = workflow.exitStatus
	email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
	email_fields['errorReport'] = (workflow.errorReport ?: 'None')
	email_fields['commandLine'] = workflow.commandLine
	email_fields['projectDir'] = workflow.projectDir
	email_fields['script_file'] = workflow.scriptFile
	email_fields['launchDir'] = workflow.launchDir
	email_fields['user'] = workflow.userName
	email_fields['Pipeline script hash ID'] = workflow.scriptId
	email_fields['kit'] = TARGETS
	email_fields['assembly'] = REF
	email_fields['manifest'] = workflow.manifest
	email_fields['summary'] = summary

	email_info = ""
  for (s in email_fields) {
	email_info += "\n${s.key}: ${s.value}"
  }

  def output_d = new File( "${params.outdir}/pipeline_info/" )
  if( !output_d.exists() ) {
      output_d.mkdirs()
  }

  def output_tf = new File( output_d, "pipeline_report.txt" )
  output_tf.withWriter { w -> w << email_info }	

 // make txt template
  def engine = new groovy.text.GStringTemplateEngine()

  def tf = new File("$baseDir/assets/email_template.txt")
  def txt_template = engine.createTemplate(tf).make(email_fields)
  def email_txt = txt_template.toString()

  // make email template
  def hf = new File("$baseDir/assets/email_template.html")
  def html_template = engine.createTemplate(hf).make(email_fields)
  def email_html = html_template.toString()
  
  def subject = "WGS analysis finished ($run_name)."

  if (params.email) {

  	def mqc_report = null
  	try {
        	if (workflow.success && !params.skip_multiqc) {
            		mqc_report = multiqc_report.getVal()
            		if (mqc_report.getClass() == ArrayList){
                		log.warn "[IKMB ExoSeq] Found multiple reports from process 'multiqc', will use only one"
                		mqc_report = mqc_report[0]
                	}
        	}
    	} catch (all) {
        	log.warn "[IKMB ExoSeq] Could not attach MultiQC report to summary email"
  	}

	def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
	def sf = new File("$baseDir/assets/sendmail_template.txt")	
    	def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    	def sendmail_html = sendmail_template.toString()

	try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
        }

  }

}
