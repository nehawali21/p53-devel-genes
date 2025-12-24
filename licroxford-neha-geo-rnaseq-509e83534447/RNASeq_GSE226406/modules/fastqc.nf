process runFastQC {
	input:
		tuple val(id), path("${id}_*.fastq.gz"), val(options) // using ${id} in input path and call to fastqc avoids having to rename output files in script block (BA)
	
	output:
		tuple val(id), path("*.html"), emit: html
		tuple val(id), path("*.zip"), emit: zip
	
	script:
	"""
		mkdir -p tempDir
		fastqc -o . \\
				 -d tempDir \\
				 ${task.cpus?"-t ${task.cpus}":""} \\
				 ${(options && options.contaminantsFile)?"-c ${options.contaminantsFile}":""} \\
				 ${(options && options.adapterFile)?"-a ${options.adapterFile}":""} \\
				 ${(options && options.kmerLength)?"-k ${options.kmerLength}":""} \\
				 ${(options && options.limitFile)?"-l ${options.limitFile}":""} \\
				 ${id}_*.fastq.gz
	"""
}
