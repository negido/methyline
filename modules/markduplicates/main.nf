process markDuplicates {
    tag "${sampleId}"
    label 'high'
    container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'
    input:
    tuple val(sampleId), path(sorted_bam)
    tuple val(sampleId), path(bai)
    output:
    tuple val(sampleId), path("${sampleId}_dedup.bam"), path("${sampleId}_dedup.bai"), emit: MARKED_DUP
    path "${sampleId}_dup-metrics.txt", emit: MD_METRICS
    script:
    """
	gatk --java-options "-Xmx${task.memory.giga}g" MarkDuplicates -I ${sorted_bam} -O ${sampleId}_dedup.bam -M ${sampleId}_dup-metrics.txt --CREATE_INDEX true
    """
}
