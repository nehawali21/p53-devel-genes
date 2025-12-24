process cramToBam {
  input:
    tuple val(id), path(cramFile), path("reference*")
  output:
    tuple val(id), path('${id}.bam'), path('${id}.bam.bai')
  script:
  """
    samtools view -T $reference -b -@ ${task.cpus} -o ${id}.bam $cramFile
    samtools index ${id}.bam
  """
}

process sortBam {
  input:
    tuple val(id), path('input.bam'), val(options)
  output:
    tuple val(id), path("${id}_sorted.bam")
  script:
  """
    samtools sort -@ ${task.cpus} -m ${task.memory?"${task.memory}".replaceAll(' ',''):"2G"} \\
             -T ${id}.tmp -O bam -o ${id}_sorted.bam \\
             ${options.byName?"-n":""} \\
             input.bam
  """
}

process indexBam {
  input:
    tuple val(id), path(input)
  output:
    tuple val(id), path(input), path("${input}.bai")
  script:
  """
    samtools index -@ ${task.cpus} $input
  """
}

process mergeBam {
  input:
    tuple val(id), path("input_?.bam"), path("input_?.bam.bai")
  output:
    tuple val(id), path("${id}_merged.bam"), path("${id}_merged.bam.bai")
  script:
  """
    set -euo pipefail
    if [ -f input_2.bam ]; then
      samtools merge -cp ${task.cpus>1?"-@ ${task.cpus}":""} ${id}_merged.bam input*.bam
      samtools index -@ ${task.cpus} ${id}_merged.bam
    else
      ln -s input_1.bam ${id}_merged.bam
      ln -s input_1.bam.bai ${id}_merged.bam.bai
    fi
    
  """
}

process viewBam {
  input:
    tuple val(id), path('input.bam'), path('input.bam.bai'), val(options)
  output:
    tuple val(id), path('output.bam')
  script:
  """
    samtools view -o output.bam -O bam -@ ${task.cpus} \
      ${options.filterBedFile?"-L ${options.filterBedFile}":""} \
      ${options.readGroup?"-r ${options.readGroup}":""} \
      ${options.readGroupFile?"-R ${options.readGroupFile}":""} \
      ${options.minMapQ?"-q ${options.minMapQ}":""} \
      ${options.library?"-l ${options.library}":""} \
      ${options.cigarOptCount?"-m ${options.cigarOptCount}":""} \
      ${options.includeFlags?"-f ${options.includeFlags}":""} \
      ${options.excludeFlags?"-F ${options.excludeFlags}":""} \
      ${options.excludeExactFlags?"-G ${options.excludeExactFlags}":""} \
      ${options.downsample?"-s ${options.downsample}":""} \
      ${options.stripTag?"-x ${options.stripTag}":""} \
      ${options.collapseBackwardCigar?"-B":""} \
      ${options.reference?"-T ${options.reference}":""} \
      input.bam ${options.region?options.region:''}
  """
}
