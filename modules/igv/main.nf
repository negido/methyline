process igv {
    tag "$sampleId"
    label 'medium'
    container 'negido/igv-reports:v1.0'
    publishDir "${params.outdir}${results_dir}/html", mode: 'copy', tags: [Type: 'html']
    input:
    tuple val(sampleId), path(bam), path(bai), path(methylation_bedgraph)
    path genome_fasta
    val referenceGenome
    path refGene
    output:
    tuple val(sampleId), path("${sampleId}_IGV.html"), emit: IGV_HTML
    script:
        max_regions = params.igv_max_regions ?: 500
    results_dir = params.analysisName ? "" : "/${sampleId}/${workflow.sessionId}"
    """
    zcat ${refGene} > refgene.refgene
        awk 'BEGIN {FS=OFS="\t"} \$0 !~ /^track/ && NF >= 3 {print \$1, \$2, \$3}' ${methylation_bedgraph} | head -n ${max_regions} > igv_regions.bed
        awk 'BEGIN {FS=OFS="\t"} \$0 !~ /^track/ && NF >= 4 {print \$1, \$2, \$3, \$4}' ${methylation_bedgraph} | head -n ${max_regions} > igv_track.bedGraph

        if [ ! -s igv_regions.bed ]; then
            awk 'BEGIN {OFS="\t"; print "chr1", 1, 1000}' > igv_regions.bed
        fi

        if [ ! -s igv_track.bedGraph ]; then
            cp ${methylation_bedgraph} igv_track.bedGraph
        fi

        create_report igv_regions.bed ${params.referenceGenome}.fa --samples ${sampleId} --tracks igv_track.bedGraph ${bam} refgene.refgene --output ${sampleId}_IGV.html
    """
}
