process featureCounts {
		
	input:
		tuple val(id), val(aligner), path("${id}.bam"), val(options)
		
	output:
		path "${id}_${aligner}_gene.featureCounts.txt", emit: counts
		path "${id}_${aligner}_gene.featureCounts.txt.summary", emit: summary
	
	script:
	"""
		if [[ "${params.paired_end}" == "true" ]]
		then
			featureCounts -p \\
				-T ${task.cpus} \\
				${options.GTFfile?"-a ${options.GTFfile}":''} \\
				${options.type?"-t ${options.type}":'exon'} \\
				${options.group?"-g ${options.group}":'gene_id'} \\
				${options.direction?"-s ${options.direction}":'2'} \\
				-o ${id}_${aligner}_gene.featureCounts.txt \\
				${id}.bam
		else
				featureCounts -T ${task.cpus} \\
				${options.GTFfile?"-a ${options.GTFfile}":''} \\
				${options.type?"-t ${options.type}":'exon'} \\
				${options.group?"-g ${options.group}":'gene_id'} \\
				${options.direction?"-s ${options.direction}":'2'} \\
				-o ${id}_${aligner}_gene.featureCounts.txt \\
				${id}.bam
	
		fi
	"""
}
