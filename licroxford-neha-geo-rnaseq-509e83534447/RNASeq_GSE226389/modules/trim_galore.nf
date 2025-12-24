process trimReads {
	
	input:
		tuple val(id), path("${id}_?.fastq.gz"), val(options)
	output:
		tuple val(id), path("*.fq.gz"), emit: fastqs
		tuple val(id), path("*.txt"), emit: report
	script:
	"""
		trim_galore --gzip \\
					-j ${task.cpus} \\
					${options.minQuality==null?'':"-q ${options.minQuality}"} \\
					${options.twoColour?"--2colour ${options.twoColour}":""} \\
					${params.minOverlap==null?'':"--stringency ${params.minOverlap}"} \\
					${options.fwAdapter?"-a ${options.fwAdapter}":''} \\
					${options.rvAdapter?"-a2 ${options.rvAdapter}":''} \\
					${options.clipR1 && options.clipR1>0?"--clip_R1 ${options.clipR1}":""} \\
					${options.clipR2 && options.clipR2>0?"--clip_R2 ${options.clipR2}":""} \\
					${options.clip3primeR1 && options.clip3primeR1>0?"--three_prime_clip_R1 ${options.clip3primeR1}":""} \\
					${options.clip3primeR2 && options.clip3primeR2>0?"--three_prime_clip_R2 ${options.clip3primeR2}":""} \\
					`if [ -e "${id}_2.fastq.gz" ]; then echo --paired; fi` \\
					${id}_*.fastq.gz
	"""
}