process starIndex {
	
	input:
		tuple val(indexDir), val(options)
	
	output:
		path "${indexDir}"
		
	script:
	"""
		STAR --runThreadN ${task.cpus} \\
			--runMode genomeGenerate \\
			--genomeDir ${indexDir} \\
			--genomeFastaFiles ${options.fa?options.fa:''} \\
			--sjdbGTFfile ${options.GTF?options.GTF:''} \\
			--sjdbOverhang ${options.overhang?options.overhang:74}
	"""	
}

process starQuant {
		
	// merges multiple lanes per sample if present and reports number of lanes
	
	input:
		path starIndex
		tuple val(sample), path("run?_1.fq.gz"), path("run?_2.fq.gz"), val(options) // works for single- and paired-end providing empty vector [] is passed to starQuant for single-end read2
	
	output:
		tuple val(sample), path("${sample}_Aligned.sortedByCoord.out.bam"), path("${sample}_Aligned.sortedByCoord.out.bam.bai"), emit: alignment
		path "*.tab", emit: tab
		path "*.out", emit: log
		tuple val(sample), env(NL), emit: lanes
		tuple val(sample), env(READ1), env(READ2), env(ID), env(PU), env(RG), emit: header
			
	script:
	"""	
		# count number of flow cell lanes to merge
		NL=\$(ls run*_1.fq.gz | wc -l)
		
		LIBRARY=${options.library?options.library:"$sample"}
		SAMPLE=${options.sample?options.sample:"$sample"}
		PL=${options.platform?options.platform:'ILLUMINA'}
		CN=${options.center?options.center:'WTCHG'}
		
		RG=""
		
		for fq_gz in ./run*_1.fq.gz
		do
			READID=\$(zcat \$fq_gz | head -n 1 | cut -d ' ' -f 1 | sed 's/\\@//')
			INSTRUMENTID=\$(echo \$READID | cut -d : -f 1)
			RUNID=\$(echo \$READID | cut -d : -f 2)
			FCID=\$(echo \$READID | cut -d : -f 3)
			LANE=\$(echo \$READID | cut -d : -f 4)
		
			ID="ID:\$RUNID:\$FCID:\$LANE"
			PU="PU:\$FCID:\$LANE"
		
			RG="\${RG}\$ID\tSM:\$SAMPLE\tLB:\$LIBRARY\tPL:\$PL\tPU:\$PU\tCN:\$CN , "
		done
		
		# remove trailing " , "
		RG=\${RG% , }
		
		# remove spaces in expansion of run*_1.fq.gz and run*_2.fq.gz but retain comma separation
		# reads must be sorted, or at least in matching order for READ1 and READ2
		# note different meaning of * in bash than in Nextflow input path where ? must be used, otherwise if only one file found "1" is not inserted
		READ1=\$(ls -mv run*_1.fq.gz | tr -d " " | tr -d "\\n")
		
		if [ -e run1_2.fq.gz ]
		then
			READ2=\$(ls -mv run*_2.fq.gz | tr -d " " | tr -d "\\n")
		else
			READ2=""
		fi
		
		if [ -z "\${READ2}" ]
		then
			STAR --genomeDir $starIndex \\
			--genomeLoad NoSharedMemory \\
			--outSAMattrRGline \${RG} \\
			--limitBAMsortRAM 12000000000 \\
			--readFilesIn \${READ1} \\
			--readFilesCommand gunzip -c \\
			--quantMode GeneCounts \\
			--runThreadN ${task.cpus} \\
			--outSAMunmapped Within \\
			--outFileNamePrefix "${sample}_" \\
			--outSAMtype BAM SortedByCoordinate
			
			READ2="not used" # needed otherwise an empty string will break the output and the process will fail, although it still gives exit status == 0
		else
		
			STAR --genomeDir $starIndex \\
				--genomeLoad NoSharedMemory \\
				--outSAMattrRGline \${RG} \\
				--limitBAMsortRAM 12000000000 \\
				--readFilesIn \${READ1} \${READ2} \\
				--readFilesCommand gunzip -c \\
				--quantMode GeneCounts \\
				--runThreadN ${task.cpus} \\
				--outSAMunmapped Within KeepPairs \\
				--outFileNamePrefix "${sample}_" \\
				--outSAMtype BAM SortedByCoordinate
		fi
			
		samtools index -@ ${task.cpus} "${sample}_Aligned.sortedByCoord.out.bam"
	"""
}
