process methylDackel {
    tag "$sampleId"
    label 'medium'
    container 'quay.io/biocontainers/methyldackel:0.6.1--h577a1d6_9'
    input:
    tuple val(sampleId), path(bam), path(bai)
    path genome_fasta

    output:
    //tuple val(sampleId), path("*_CpG.txt"), emit: methylation
    tuple val(sampleId), path("*.bedGraph"), emit: bedgraph

    script:
    """
    MethylDackel extract -o ${sampleId} ${params.referenceGenome}.fa ${bam}
    """
}
