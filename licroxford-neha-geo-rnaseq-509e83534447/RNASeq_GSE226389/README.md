RNAseq processing pipeline
==========================

## Summary

This version is a clone of the paper_OAC_Trial_RNASeq pipeline created for the first trial publication and reatisn its version number.  It has been edited to work with an input table
prepared by CombineMetadata/combine_metadata_v3 for the C and D cohorts and referenced in nextflow.config as params.inputTable.

This nextflow pipeline will process fastq.gz files from bulk RNA sequencing, producing counts tables using STAR and featureCounts and
- per-read and per-alignment QC files
- quality- and adapter-trimmed, and deduplicated bam files
- MultiQC report with general statistics, FastQC for raw and quality- and adapter-trimmed *.fastq.gz (sequence counts, sequence quality histograms, 
per sequence quality scores, per base sequence content, per sequence GC content, per base N content, sequence length distribution, 
sequence duplication levels, overrepresented sequences, adapter content, status checks)
- MultiQC report with general statistics, STAR alignment Scores, RSeQC %rRNA, Samtools %rRNA, Subread featureCounts, RnaSeqMetrics (Picard: primary alignments,
strand mapping, gene coverage), MultipleMetrics (Picard: alignment summary, mean read length, base quality distribution)


## Requirements

### Software
- [Nextflow](https://www.nextflow.io/)
- [TrimGalore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- [samtools](https://www.htslib.org/)
- [GATK tools](https://gatk.broadinstitute.org/hc/en-us)
- [Picard tools](https://broadinstitute.github.io/picard/) (now part of GATK tools anyway)
- [STAR](https://github.com/alexdobin/STAR)
- [featureCounts](http://subread.sourceforge.net )
- [MultiQC](https://multiqc.info/)
- [R](https://www.r-project.org/)

### Data
The fastq files should be named as downloaded from OGC. Paths to fastqs are constructed from `params.dataDir` defined in nextflow.config and the R1.File and R2.File 
columns in the input table, with `params.dataDir` containing either the actual fastq files, or symlinks to them.

### Running the pipeline

The nextflow.config reference genome files and associated bed files in 
/mnt/lustre/users/ramess/Code/nxf_workflows/ramess/paper_OAC_Trial_RNASeq/ref

The pipeline should be run from the script (workflow) directory with

`nextflow -log <your_path_to_logs>/nextflow.log run run.nf`

The pipeline would likely run faster after removing, or increasing, the maxForks lines in nextflow.config.

Please note that the building of STAR's reference index is conditional and only occurs if not already at the location specified by `"${params.starIndexDir}/${params.starIndex}"`.

