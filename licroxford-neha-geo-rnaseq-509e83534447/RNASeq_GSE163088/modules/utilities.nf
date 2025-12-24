process gzip {
	input:
		tuple val(id), path(inputFile)
	output:
		tuple val(id), path("${inputFile}.gz")
	when:
		!inputFile.name.endsWith(".gz")
	script:
	"""
		if [ \$(file $inputFile | fgrep -c "gzip compressed") -gt "0" ]; then
			ln -s $inputFile ${inputFile}.gz
		else
			if [[ "${task.cpus}" > 1 ]]; then
				pigz -f -p ${task.cpus} $inputFile
			else
				gzip $inputFile
			fi
		fi
	"""
}

process pickleScriptVersion {
	input:
		tuple val(version), path("*")
	output:
		path("version_${version}.tar.gz")
	script:
	"""
		tar --exclude-vcs --exclude='.*' -hczf version_${version}.tar.gz *
	"""
}