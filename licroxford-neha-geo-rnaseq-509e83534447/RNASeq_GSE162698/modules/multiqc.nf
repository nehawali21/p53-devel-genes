process runMultiQC {
	
	// Input paths are not referenced directly in script block as MultiQC will look for appropriate files in current directory recursively
	// Subdirectories are included in input paths to avoid any possible filename clashes when Nextflow stages these to the work subdirectory
	// Channels used for input should be collections i.e. use collect() operator so that MultiQC only runs when all instances of a process linked to each input channel have completed
	// When preparing data in workflow, use an ifEmpty { [] } closure after collect() to avoid halting workflow in case any channel is empty
	// If a tool is not being used then provide empty list/array with [] in call - note, keep order as in input block
	
	// Include publishDir block here as ${pubdir} is passed to process via input
	
	publishDir pattern: "{multiqc_data,multiqc_report.html}", mode: "link", path: "${params.resultDir}/${pubdir}"
	
	input:
		val pubdir
		val change_names
		path 'fastqc/*'
		path 'star/*'
		path 'featureCounts/*'
		path 'rnaseqmetrics_star/*'
		path 'multiplemetrics_star/*'
		path 'rseqc_rrna/*'
		path 'samtools_rrna/*'

	output:
		path "multiqc_report.html", emit: html
		path "multiqc_data", emit: data

	script:
	// include --sample-names ${params.sample_names} to provide alternative sample names, selectable via buttons, in tab separated file in which headers are button names
	// and first column contains the sample names initally displayed
	"""
		if [[ $change_names = true ]]; then
		
			python3 -m multiqc --config ${params.multiqc_config} --sample-names ${params.sample_names} .
			
		else
	
			python3 -m multiqc --config ${params.multiqc_config} .
		fi
	"""
}
