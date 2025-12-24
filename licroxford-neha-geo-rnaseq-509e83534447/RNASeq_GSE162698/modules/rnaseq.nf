// calculate percentage of rRNA reads - method 1 - alternative to CollectRnaSeqMetrics which calculates percentage of rRNA bases wrt to input bed file 
process rseqc {
	
	publishDir pattern: "${id}_rRNA_QC.txt", mode: "link", path: "${params.resultDir}/RSeQC"
	
	input:
		tuple val(id), val(aligner),  path("input.bam"), path("input.bam.bai"), path ("input.bed")
	
	output:
		path "${id}_rRNA_QC.txt", emit: report
		tuple val(id), val(aligner), env(TOTALCOUNT), env(NONRIBOCOUNT), env(non_rRNA_pct), env(RIBOCOUNT), env(rRNA_pct), emit: values

	shell:
		
	'''
		if [ ! -f !{params.resultDir}/RSeQC/!{aligner}_rRNA_QC.txt ]; then
		
			if [ ! -d !{params.resultDir}/RSeQC ]; then
				mkdir !{params.resultDir}/RSeQC
			fi
				
			echo -e "Sample\tnon-rRNA Count\trRNA Count\trRNA pct" > !{params.resultDir}/RSeQC/!{aligner}_rRNA_QC.txt
		fi
	
		split_bam.py -i input.bam -r input.bed -o !{id}_output
		
		NONRIBOCOUNT=`samtools view -@ !{task.cpus} -c -q 255 !{id}_output.ex.bam`
		RIBOCOUNT=`samtools view -@ !{task.cpus} -c -q 255 !{id}_output.in.bam`
		TOTALCOUNT=$((RIBOCOUNT + NONRIBOCOUNT))
	
		non_rRNA_pct=`python3 -c "print('{:.3f}'.format(100 * $NONRIBOCOUNT / ($TOTALCOUNT)))"`
		rRNA_pct=`python3 -c "print('{:.3f}'.format(100 * $RIBOCOUNT / ($TOTALCOUNT)))"`

		
		env printf "!{id}\t%u\t%u\t%s\t%u\t%s" $TOTALCOUNT $NONRIBOCOUNT $non_rRNA_pct $RIBOCOUNT $rRNA_pct >> !{params.resultDir}/RSeQC/!{aligner}_rRNA_QC.txt
		echo "" >> !{params.resultDir}/RSeQC/!{aligner}_rRNA_QC.txt
	
		echo -e "Aligner\tTotal Count\tnon-rRNA Count\tnon-rRNA pct\trRNA Count\trRNA pct" > !{id}_rRNA_QC.txt
		env printf "!{aligner}\t%u\t%u\t%s\t%u\t%s" $TOTALCOUNT $NONRIBOCOUNT $non_rRNA_pct $RIBOCOUNT $rRNA_pct >> !{id}_rRNA_QC.txt
		echo ""  >> !{id}_rRNA_QC.txt
	'''
}


// calculate percentage of rRNA reads - method 2 - alternative to CollectRnaSeqMetrics which calculates percentage of rRNA bases wrt to input bed file 
process samtools_rrna {

	publishDir pattern: "${id}_rRNA_QC.txt", mode: "link", path: "${params.resultDir}/samtools_qc"
	
	input:
		tuple val(id), val(aligner),  path("input.bam"), path("input.bam.bai"), path("all_rRNA_bedFile"), path("rRNA_bedFile"), path("mt_rRNA_bedFile"), path("mt_rRNA_1_bedFile"), path("mt_rRNA_2_bedFile"), path("pseudo_rRNA_bedFile")
	
	output:
		path "${id}_rRNA_QC.txt", emit: report
		tuple val(id), val(aligner), env(TOTALCOUNT), env(NONRIBOCOUNT), env(all_rRNA_COUNT), env(rRNA_COUNT), env(mt_rRNA_COUNT), env(mt_rRNA_1_COUNT), env(mt_rRNA_2_COUNT), env(pseudo_rRNA_COUNT), \
					env(non_RIBO_pct), env(all_rRNA_pct), env(rRNA_pct), env(mt_rRNA_pct), env(mt_rRNA_1_pct), env(mt_rRNA_2_pct), env(pseudo_rRNA_pct), emit: values
		
	shell:
	'''	
		if [ ! -f !{params.resultDir}/samtools_qc/!{aligner}_rRNA_QC.txt ]; then
		
			if [ ! -d !{params.resultDir}/samtools_qc ]; then
				mkdir !{params.resultDir}/samtools_qc
			fi
				
			echo -e "Sample\tTotal_Count\tnon-rRNA_Count\tall_rRNA_Count\trRNA_Count\tmt_rRNA_Count\tmt_rRNA_1_Count\tmt_rRNA_2_Count\tpseudo_rRNA_Count\tnon_rRNA_pct\tall_rRNA_pct\trRNA_pct\tmt_rRNA_pct\tmt_rRNA_1_pct\tmt_rRNA_2_pct\tpseudo_rRNA_pct" > !{params.resultDir}/samtools_qc/!{aligner}_rRNA_QC.txt
		fi
	
		TOTALCOUNT=`samtools view -q 255 -@ !{task.cpus} -h -c input.bam`
	
		all_rRNA_COUNT=`samtools view -q 255 -L all_rRNA_bedFile -@ !{task.cpus} -h -c input.bam`
	
		rRNA_COUNT=`samtools view -q 255 -L rRNA_bedFile -@ !{task.cpus} -h -c input.bam`
	
		mt_rRNA_COUNT=`samtools view -q 255 -L mt_rRNA_bedFile -@ !{task.cpus} -h -c input.bam`
	
		mt_rRNA_1_COUNT=`samtools view -q 255 -L mt_rRNA_1_bedFile -@ !{task.cpus} -h -c input.bam`
	
		mt_rRNA_2_COUNT=`samtools view -q 255 -L mt_rRNA_2_bedFile -@ !{task.cpus} -h -c input.bam`
	
		pseudo_rRNA_COUNT=`samtools view -q 255 -L pseudo_rRNA_bedFile -@ !{task.cpus} -h -c input.bam`
		
		NONRIBOCOUNT=$((TOTALCOUNT - all_rRNA_COUNT))
	
	
		non_RIBO_pct=`python3 -c "print('{:.4f}'.format(100 * $NONRIBOCOUNT / $TOTALCOUNT))"`
	
		all_rRNA_pct=`python3 -c "print('{:.4f}'.format(100 * $all_rRNA_COUNT / $TOTALCOUNT))"`
	
		rRNA_pct=`python3 -c "print('{:.4f}'.format(100 * $rRNA_COUNT / $TOTALCOUNT))"`
	
		mt_rRNA_pct=`python3 -c "print('{:.4f}'.format(100 * $mt_rRNA_COUNT / $TOTALCOUNT))"`
	
		mt_rRNA_1_pct=`python3 -c "print('{:.4f}'.format(100 * $mt_rRNA_1_COUNT / $TOTALCOUNT))"`
	
		mt_rRNA_2_pct=`python3 -c "print('{:.4f}'.format(100 * $mt_rRNA_2_COUNT / $TOTALCOUNT))"`
	
		pseudo_rRNA_pct=`python3 -c "print('{:.4f}'.format(100 * $pseudo_rRNA_COUNT / $TOTALCOUNT))"`
	
	
		env printf "!{id}\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%s\t%s\t%s\t%s\t%s\t%s\t%s" $TOTALCOUNT $NONRIBOCOUNT $all_rRNA_COUNT $rRNA_COUNT $mt_rRNA_COUNT $mt_rRNA_1_COUNT $mt_rRNA_2_COUNT $pseudo_rRNA_COUNT $non_RIBO_pct $all_rRNA_pct $rRNA_pct $mt_rRNA_pct $mt_rRNA_1_pct $mt_rRNA_2_pct $pseudo_rRNA_pct >> !{params.resultDir}/samtools_qc/!{aligner}_rRNA_QC.txt
		echo "" >> !{params.resultDir}/samtools_qc/!{aligner}_rRNA_QC.txt
	
	
		echo -e "Aligner\tTotal_Count\tnon_rRNA_Count\tall_rRNA_Count\trRNA_Count\tmt_rRNA_Count\tmt_rRNA_1_Count\tmt_rRNA_2_Count\tpseudo_rRNA_Count\tnon_Ribo_pct\tall_rRNA_pct\trRNA_pct\tmt_rRNA_pct\tmt_rRNA_1_pct\tmt_rRNA_2_pct\tpseudo_rRNA_pct" > !{id}_rRNA_QC.txt
		env printf "!{aligner}\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%s\t%s\t%s\t%s\t%s\t%s\t%s" $TOTALCOUNT $NONRIBOCOUNT $all_rRNA_COUNT $rRNA_COUNT $mt_rRNA_COUNT $mt_rRNA_1_COUNT $mt_rRNA_2_COUNT $pseudo_rRNA_COUNT $non_RIBO_pct $all_rRNA_pct $rRNA_pct $mt_rRNA_pct $mt_rRNA_1_pct $mt_rRNA_2_pct $pseudo_rRNA_pct >> !{id}_rRNA_QC.txt
		echo ""  >> !{id}_rRNA_QC.txt
	'''
}

// mark duplicates - publishes metrics, merged and deduplicated bam/bai
process markDuplicates {
	
	// publishDir blocks needed here because nextflow.config will not know value of $aligner
	publishDir pattern: "${id}_md_report.txt", mode: "link", path: "${params.resultDir}/$aligner/QCMetrics"
	publishDir pattern: "${id}_md.{bam, bam.bai}", mode: "link", path: "${params.resultDir}/$aligner/merged_md_bam"
	
	input:
		tuple val(id), val(aligner),  path("input.bam"), path("input.bam.bai")

	output:
		tuple val(id), path("${id}_md.bam"), path("${id}_md.bam.bai"), emit: alignment
		path "${id}_md_report.txt", emit: report
	
	script:
	"""
		set -euo pipefail
		picard --java-options "-Djava.io.tmpdir=\$PWD" MarkDuplicates -I input.bam -O ${id}_md.bam -M ${id}_md_report.txt --ASO coordinate --TMP_DIR . -R ${params.genomeFile}
		samtools index ${id}_md.bam
	"""
}

process collectRnaSeqMetrics {

	publishDir pattern: "${id}_${aligner}_output.RNA_Metrics", mode: 'link', path: "${params.resultDir}/RnaSeqMetrics"
	
	input:
		tuple val(id), val(aligner), path("${id}.bam"), val(options)
	
	output:
		path "${id}_${aligner}_output.RNA_Metrics"
	
	script:
	"""
		picard --java-options "-Xmx14g -Djava.io.tmpdir=\$PWD" \\
				CollectRnaSeqMetrics \\
				-I ${id}.bam \\
				-O ${id}_${aligner}_output.RNA_Metrics \\
				${options.ref_flat?"--REF_FLAT ${options.ref_flat}":''} \\
				${options.strand?"--STRAND ${options.strand}":''} \\
				${options.ribo_intervals?"--RIBOSOMAL_INTERVALS ${options.ribo_intervals}":''}
				
	
		if [ ! -f ${params.resultDir}/RnaSeqMetrics/${aligner}_RNA_Metrics.tsv ]; then
		
			if [ ! -d ${params.resultDir}/RnaSeqMetrics ]; then
				mkdir ${params.resultDir}/RnaSeqMetrics
			fi
				
			echo -e "PF_BASES\tPF_ALIGNED_BASES\tRIBOSOMAL_BASES\tCODING_BASES\tUTR_BASES\tINTRONIC_BASES\tINTERGENIC_BASES\tIGNORED_READS\tCORRECT_STRAND_READS\tINCORRECT_STRAND_READS\tNUM_R1_TRANSCRIPT_STRAND_READS\tNUM_R2_TRANSCRIPT_STRAND_READS\tNUM_UNEXPLAINED_READS\tPCT_R1_TRANSCRIPT_STRAND_READS\tPCT_R2_TRANSCRIPT_STRAND_READS\tPCT_RIBOSOMAL_BASES\tPCT_CODING_BASES\tPCT_UTR_BASES\tPCT_INTRONIC_BASES\tPCT_INTERGENIC_BASES\tPCT_MRNA_BASES\tPCT_USABLE_BASES\tPCT_CORRECT_STRAND_READS\tMEDIAN_CV_COVERAGE\tMEDIAN_5PRIME_BIAS\tMEDIAN_3PRIME_BIAS\tMEDIAN_5PRIME_TO_3PRIME_BIAS\tSAMPLE\tLIBRARY\tREAD_GROUP" \\
				> ${params.resultDir}/RnaSeqMetrics/${aligner}_RNA_Metrics.tsv
		fi
				
		NUM=`awk '/## METRICS CLASS/{ print NR; exit }' ${id}_${aligner}_output.RNA_Metrics`
		NUM=\$((NUM + 2))
		LINE=`sed "\${NUM}q;d" ${id}_${aligner}_output.RNA_Metrics`
		
		echo -e \$LINE >> ${params.resultDir}/RnaSeqMetrics/${aligner}_RNA_Metrics.tsv
	"""
}


process collectMultipleMetrics {

	input:
	  tuple val(id), val(aligner), path("${id}.bam"), path("${id}.bam.bai"), val(options)
	  
	output:
	  tuple val(id), path("${id}_${aligner}*")
	
	script:
	"""
		picard --java-options "-Xmx14g -Djava.io.tmpdir=\$PWD" \\
				CollectMultipleMetrics \\
				-I ${id}.bam \\
				-O ${id}_${aligner} \\
				--FILE_EXTENSION ".txt" \\
				-R ${params.genomeFile} \\
				${options.intervals?"--INTERVALS ${options.intervals}":''} \\
				${options.dbSNP?"--DB_SNP ${options.dbSNP}":''} \\
				--PROGRAM null \\
				--PROGRAM CollectSequencingArtifactMetrics \\
				--PROGRAM CollectAlignmentSummaryMetrics \\
				--PROGRAM CollectInsertSizeMetrics \\
				--PROGRAM QualityScoreDistribution
	"""
}
