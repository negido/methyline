process bwamethAlign {
    tag "$sampleId"
    label 'high'
    container 'quay.io/biocontainers/bwameth:0.2.9--pyh7e72e81_0'
    input:
    tuple val(sampleId), path(reads)
    path genome_fasta
    val referenceGenome
    output:
    tuple val(sampleId), path("${sampleId}.bam"), emit: bam
    tuple val(sampleId), path("${sampleId}.bam.bai"), emit: bai
    script:
    """
    /usr/local/bin/bwameth.py  -t ${task.cpus} --reference ${referenceGenome}.fa ${reads} | samtools sort -o ${sampleId}.bam -
    samtools index ${sampleId}.bam
    samtools flagstat ${sampleId}.bam > ${sampleId}_flagstat_report.txt
    samtools stats ${sampleId}.bam > ${sampleId}_stats_report.txt
    """
}
