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
    tuple val(sampleId), path("*.cytosine_report_filter.txt"), emit: cytosine_report

    script:
    """
    MethylDackel extract --mergeContext -o ${sampleId} -@ ${task.cpus} ${params.referenceGenome}.fa ${bam}
    MethylDackel extract  --cytosine_report -@ ${task.cpus} -o ${sampleId} ${params.referenceGenome}.fa ${bam}
    awk '\$4+\$5 > 0' ${sampleId}.cytosine_report.txt > ${sampleId}.cytosine_report_filter.txt 
    """
}
