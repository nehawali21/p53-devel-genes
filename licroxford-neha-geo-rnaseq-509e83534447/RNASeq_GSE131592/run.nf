#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Nextflow workflow script to process RNASeq fastq.gz

// run from script (workflow) directory with e.g.

// NXF_VER=20.10.0 nextflow -log <path_to_logs>/nextflow.log run run.nf -with-report -with-trace -with-timeline

// e.g.

// cd ~/Code/nxf_workflows/ramess/OAC_Trial_RNASeq_CD

// with /mnt/lustre/archive/Lu_lab/OAC_Immuno_Trial/RNASeq/P230128
// NXF_VER=20.10.0 nextflow -log /mnt/lustre/users/ramess/OAC_Trial_RNASeq_CD/Log/test_1_nextflow.log run run.nf


// Notes

// Building of reference indexes is conditional and only occurs if not at location specified in params

// When MultiQC is provided with featureCounts output from two aligners e.g. STAR and HISAT2, it thinks it has duplicate data and ignores one set of counts
// to avoid this the name of the alignment tool is included in the names of the outputs from featureCounts and path_filters were added to module_order blocks in multiqc_config.yaml 

// order in module_order section of multiqc_config.yaml is not followed, for unknown reason


// include modules

include { trimReads } from './modules/trim_galore' params(params)
include { runFastQC } from './modules/fastqc' params(params)
include { starIndex; starQuant } from './modules/star' params(params)
include { featureCounts } from './modules/featurecounts' params(params)
include { collectMultipleMetrics; collectRnaSeqMetrics; markDuplicates; rseqc; samtools_rrna } from './modules/rnaseq' params(params)
include { runMultiQC } from './modules/multiqc' params(params)
include { runMultiQC as runMultiQC2 } from './modules/multiqc' params(params)

include { pickleScriptVersion } from './modules/utilities'

// define variables and initialize files

def star_index = new File( "${params.starIndexDir}/${params.starIndex}" )

//def nlanes_data_file = new File( "${params.resultDir}/STAR/nlanes_mqc.tsv" )

def result_dir = new File( "${params.resultDir}" )

if( !result_dir.exists() ) {result_dir.mkdir()}


// define subworkflowss

workflow STAR {
		
	take: 
		sample_fastqs
		
	main:
	
		if( "${params.paired_end}" == true ) {
			sample_fastqs \
				| multiMap { sample, id, reads, options ->
					read1: [sample, reads[0]]
					read2: [sample, reads[1]]
					options: [sample, options]
					} \
				| set { sample_reads }
		
			sample_reads.read1 \
				| groupTuple(sort: true) \
				| set { sample_read1 }
				
			sample_reads.read2 \
				| groupTuple(sort: true) \
				| set { sample_read2 }
				
			sample_read1.join(sample_read2) \
				| join(sample_reads.options) \
				| set { sample_group_reads }
				
		} else {
                        // sample_fastqs | view

			// Note that reads here is just a single item rather than
			// a list, as the path output in trimReads only returns
			// a list if there is more than one.
			sample_fastqs \
			| multiMap { sample, id, reads, options ->
				read1: [sample, reads]
				options: [sample, options]
				} \
			| set { sample_reads }
	
			// sample_reads.read1 \
			// 	| groupTuple(sort: true) \
			// 	| join(sample_reads.options) \
			// 	| set { sample_group_reads }
		
			sample_reads.read1 \
				| groupTuple(sort: true) \
				| map {sample, reads -> [sample, reads, []] } \
				| join(sample_reads.options) | view \
				| set { sample_group_reads }
		}
		
	 	if( !star_index.exists() ) {
			// create index to genome
			starIndex(["${params.starIndex}", [fa: "${params.genomeFastaFiles}", GTF: "${params.sjdbGTFfile}", overhang: "${params.sjdbOverhang}"]])
			// align to genome
			starQuant(starIndex.out, sample_group_reads)
	 	} else {
			// align to genome
			 starQuant("${params.starIndexDir}/${params.starIndex}", sample_group_reads)
	 	}
		 
		starQuant.out.tab | collect | ifEmpty { [] } | set { tab }
		starQuant.out.log | collect | ifEmpty { [] } | set { logs }
		
		// mark duplicates  - need to pass aligner name
		starQuant.out.alignment \
			| map { sample, bam, bai -> [sample, "STAR", bam, bai] } \
			| markDuplicates
			
		markDuplicates.out.alignment \
			| map { sample, bam, bai -> [sample, "STAR", bam, bai] } \
			| set { md_bam_bai }
			
	emit:
		starQuant.out.alignment; md_bam_bai; tab; logs; starQuant.out.lanes
}


workflow FEATURECOUNTS {
	
	take:
		md_bam_bai
	
	main:	

		md_bam_bai \
			| map { sample, aligner, bam, bai -> [sample, aligner, bam] } \
			| set { tup_aligner_md_bam }
		
		tup_aligner_md_bam \
			| map { sample, aligner, bam -> [sample, aligner, bam, [GTFfile: "${params.GTFfile}", direction: "${params.direction}", type: "${params.type}", group: "${params.group}"]] } \
			| set { all_bam_fc }
			
		featureCounts(all_bam_fc)
		
		featureCounts.out.counts | collect | ifEmpty { [] } | set { counts }
		featureCounts.out.summary | collect | ifEmpty { [] } | set { summary }
		
	emit:
		counts; summary	
}



workflow COLLECTRNASEQMETRICS {
	
	take:
		aligner_md_bam_bai
	
	main:
	
		aligner_md_bam_bai \
			| map { id, aligner, bam, bai -> [id, aligner, bam, [ribo_intervals: "${params.ribo_intervals}", ref_flat: "${params.ref_flat}", strand: "${params.strand}"]] } \
			| set { tup_aligner_md_bam_bai }
		
		collectRnaSeqMetrics(tup_aligner_md_bam_bai) \
			| collect | ifEmpty { [] } | set { metrics }
		
	emit:
		metrics
}



workflow COLLECTMULTIPLEMETRICS {
	
	take:
		aligner_md_bam_bai

	main:

		aligner_md_bam_bai \
			| map { id, aligner, bam, bai -> [id, aligner, bam, bai, [dbSNP: "${params.dbSNP}"]] } \
			| set { tup_aligner_md_bam_bai }
			
		collectMultipleMetrics(tup_aligner_md_bam_bai) \
			| map { it[1] } \
			| filter { it =~ /txt/ } \
			| collect | ifEmpty { [] } | set { metrics }
		
	emit:
		metrics
}


workflow RSEQC {
	
	take:
		aligner_md_bams
	
	main:
	
		aligner_md_bams \
			| map { id, aligner, bam, bai -> [id, aligner, bam, bai, "${params.all_rRNA_bedFile}"] } \
			| set { aligner_md_bams_bed }
	
		rseqc(aligner_md_bams_bed)
	
	emit:
		rseqc.out.report; rseqc.out.values
}


workflow FASTQC {
	
	take:
		fastqs
		
	main:
	
		fastqs \
			| map { tuple(it[0], it[1], [kmerLength: "${params.kmerLength}"]) } \
			| runFastQC
	
		runFastQC.out.zip \
			| map { it[1] } \
			| collect | ifEmpty { [] } | set { fastqc_zip }
		
	emit:
		fastqc_zip
}



workflow SAMTOOLS_RRNA {
	
	take:
		aligner_md_bams

	main:
			
		aligner_md_bams
			| map { id, aligner, bam, bai -> [id, aligner, bam, bai, "${params.all_rRNA_bedFile}", "${params.rRNA_bedFile}", "${params.mt_rRNA_bedFile}", "${params.mt_rRNA_1_bedFile}", "${params.mt_rRNA_2_bedFile}", "${params.pseudo_rRNA_bedFile}"] } \
			| set { aligner_md_bams_bed }
	
		samtools_rrna(aligner_md_bams_bed)
	
	emit:
		samtools_rrna.out.report; samtools_rrna.out.values
}


workflow MULTIQC {

	take:
		pubdir
		change_names
		fastqc_out
		star_out
		featurecounts_out
		rnaseqmetrics_star_out
		multiplemetrics_star_out
		rseqc_metrics
		samtools_rrna_metrics
				
	main:
		runMultiQC(pubdir, change_names, fastqc_out, star_out, featurecounts_out, rnaseqmetrics_star_out, multiplemetrics_star_out, rseqc_metrics, samtools_rrna_metrics)
}

workflow MULTIQC_FQC {
	
		take:
			pubdir
			change_names
			fastqc_out
			star_out
			featurecounts_out
			rnaseqmetrics_star_out
			multiplemetrics_star_out
			rseqc_metrics
			samtools_rrna_metrics
						
		main:
			runMultiQC2(pubdir, change_names, fastqc_out, star_out, featurecounts_out, rnaseqmetrics_star_out, multiplemetrics_star_out, rseqc_metrics, samtools_rrna_metrics)
}


// define main workflow

workflow {
	main:
	
		// Save the script and important files in the resultDir
		Channel.value([params.version, [projectDir, params.inputTable]]) \
			| pickleScriptVersion

		// read in metadata table and select fields/columns (for inputTable prepared for wf by CombineMetadata/combine_metadata_v3.R only)
		Channel.fromPath("${params.inputTable}") \
			| splitCsv(sep: "\t", header: true, strip: true) \
			| multiMap { it ->
				sample_names: [it["Read.Group"], it["Submission.Number"] + '_' + it.Sample, it.Sample, it["Source.Label"]]
				read_group_sample: [it["Read.Group"], it["Source.Label"], it.Centre]
				fastqs: [it["Read.Group"], ["${params.dataDir}/" + it["R1.File"]]]
				} \
			| set { metadata }
			
	  	metadata.sample_names \
			| map { read, experiment, sample, subject_id -> ["trimmed_" + read, "trimmed_" + experiment, "trimmed_" + sample, "trimmed_" + subject_id] } \
			| set { trimmed_sample_names }
			 
		metadata.sample_names.concat(trimmed_sample_names) \
			| map { read, experiment, sample, subject_id -> "$read\t$experiment\t$sample\t$subject_id"  } \
			| collectFile(name: "${params.sample_names}", seed: 'Run\tExperiment\tSample\tSubject_ID', newLine: true)
			
		// add options required by trimReads and run
		metadata.fastqs \
			| map { id, paths -> tuple(id, paths, [twoColour: params.minQuality, minOverlap: params.minOverlap, fwAdapter: params.fwAdapter, clipR1: params.trimLength.toInteger()])} \
			| trimReads	
	
		metadata.read_group_sample \
			| join(trimReads.out.fastqs) \
			| map { id, sample, centre, fastqs -> [sample, id, fastqs, [sample: "$sample", library: "$sample", platform: "ILLUMINA", center: "$centre"]] } \
			| set { sample_fastqs }
		
		// update id with prefix "trimmed_" or suffix "_trimmed"- if this is not included within the file name, MultiQC ignores thinking they are duplicates of the raw untrimmed files
		trimReads.out.fastqs \
			| map { id, paths  -> tuple("$id" + "_trimmed", paths) } \
			| set { trimmed_fastqs }
		
		trimmed_fastqs.concat(metadata.fastqs) \
			| set { all_fastqs  }
			
		if( "${params.fastqc}" == true ) {
			// run FastQC	
			FASTQC(all_fastqs)
			
			FASTQC.out \
				| set { fastqc }
				
		} else {
			Channel.empty() \
				| collect | ifEmpty { [] } | set { fastqc }
		}
		
		if( "${params.star}" == true ) {
			// align fastqs to genome with STAR
			STAR(sample_fastqs)

		}	
	
		if( "${params.featurecounts}" == true ) {
			// featureCounts from subRead	
			if( "${params.star}" == true ) {
				
				FEATURECOUNTS(STAR.out[1])
			}
		}
	
		if( "${params.rnaseqmetrics}" == true ) {
			COLLECTRNASEQMETRICS(STAR.out[1])
		}
			
		if( "${params.multiplemetrics}" == true ) {
			COLLECTMULTIPLEMETRICS(STAR.out[1])
		}
			
		if( "${params.rseqc}" == true ) {
			RSEQC(STAR.out[1])
		}
		
		if( "${params.samtools_rrna}" == true ) {
			SAMTOOLS_RRNA(STAR.out[1])
		}		
			
		if( "${params.star}" == true ) {
			
			if( "${params.rseqc}" == true ) {
			
				STAR.out[4] \
					| join(RSEQC.out[1]) \
					| map { id, nlanes, aligner, Total_Count, non_rRNA_Count, non_rRNA_pct, rRNA_Count, rRNA_pct -> "$id\t$nlanes\t$aligner\t$Total_Count\t$non_rRNA_Count\t$non_rRNA_pct\t$rRNA_Count\t$rRNA_pct" } \
					| collectFile(name: "${params.resultDir}/RSeQC/rRNA_QC_mqc.tsv", seed: "Sample\tLanes\tAligner\tTotal Count\tnon-rRNA Count\tnon-rRNA pct\trRNA Count\trRNA pct", sort: true, newLine: true) \
					| set { rseqc_metrics_data }
			
			} else {
				Channel.empty() \
					| collect | ifEmpty { [] } | set { rseqc_metrics_data }
			}
			
			
			if( "${params.samtools_rrna}" == true ) {
				
				STAR.out[4] \
					| join(SAMTOOLS_RRNA.out[1]) \
					| map { id, nlanes, aligner, Total_COUNT, non_rRNA_COUNT, all_rRNA_COUNT, rRNA_COUNT, mt_rRNA_COUNT, mt_rRNA_1_COUNT, mt_rRNA_2_COUNT, pseudo_rRNA_COUNT, non_rRNA_pct, all_rRNA_pct, rRNA_pct, mt_rRNA_pct, mt_rRNA_1_pct, mt_rRNA_2_pct, pseudo_rRNA_pct ->
						"$id\t$nlanes\t$aligner\t$Total_COUNT\t$non_rRNA_COUNT\t$all_rRNA_COUNT\t$rRNA_COUNT\t$mt_rRNA_COUNT\t$mt_rRNA_1_COUNT\t$mt_rRNA_2_COUNT\t$pseudo_rRNA_COUNT\t$non_rRNA_pct\t$all_rRNA_pct\t$rRNA_pct\t$mt_rRNA_pct\t$mt_rRNA_1_pct\t$mt_rRNA_2_pct\t$pseudo_rRNA_pct" } \
					| collectFile(name: "${params.resultDir}/samtools_qc/st_rRNA_QC_mqc.tsv", seed: 'Sample\tLanes\tAligner\tTotal Count\tnon-rRNA Count\tall rRNA Count\trRNA Count\tmt rRNA Count\tmt rRNA 1 Count\tmt rRNA 2 Count\tpseudo rRNA Count\tnon-rRNA pct\tall rRNA pct\trRNA pct\tmt rRNA pct\tmt rRNA 1 pct\tmt rRNA 2 pct\tpseudo rRNA pct', sort: true, newLine: true) \
					| set { samtools_rrna_metrics }
			
			} else {
				Channel.empty() \
					| collect | ifEmpty { [] } | set { samtools_rrna_metrics }
			}
			
			if( "${params.nlanes}" == true ) {
				STAR.out[4] \
					| map { sample, nlanes -> "$sample\t$nlanes" } \
					| collectFile(name: "${params.resultDir}/STAR/nlanes//nlanes_mqc.tsv", seed: "Sample\tLanes", sort: true, newLine: true) \
					| set { nlanes_data }
					
			} else {
				Channel.empty() \
					| collect | ifEmpty { [] } | set { nlanes_data }
			}

		}
		
		
		if( "${params.multiqc}" == true ) {
			
			// generate reports
			MULTIQC("multiqc", true, [], STAR.out[3], FEATURECOUNTS.out[1], COLLECTRNASEQMETRICS.out, COLLECTMULTIPLEMETRICS.out, rseqc_metrics_data, samtools_rrna_metrics)
			MULTIQC_FQC("multiqc_fastqc", true, fastqc, [], [], [], [], [], [])
			
		}

}
