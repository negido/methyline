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
    tuple val(dmr), path(dmr_bed)
    output:
    tuple val(sampleId), path("${sampleId}_IGV.html"), emit: IGV_HTML
    script:
    def max_regions = params.igv_max_regions ?: 500
    def results_dir = params.analysisName ?: "${sampleId}/${workflow.sessionId}"
    def has_dmrs = dmr_bed != null && dmr_bed.name != 'NO_FILE'
    """
    # Convertir refGene (UCSC) a BED estándar: chrom txStart txEnd name
    # Columnas refGene: bin name chrom strand txStart txEnd ...
    zcat ${refGene} \
    | awk 'BEGIN {FS=OFS="\\t"} !/^#/ && NF >= 6 {print \$3, \$5, \$6, \$2}' \
    > refgene.bed

    # Si hay DMRs (modo multisample), usar como regiones de interés
    # Si no, muestrear aleatoriamente del bedGraph (modo single-sample)
    if ${has_dmrs}; then
        cp ${dmr_bed} igv_regions.bed
    else
        grep -v "^track" ${methylation_bedgraph} \
            | shuf -n ${max_regions} \
            | sort -k1,1 -k2,2n \
            | awk 'BEGIN {FS=OFS="\\t"} NF >= 3 {print \$1, \$2, \$3}' \
            > igv_regions.bed
    fi

    grep -v "^track" ${methylation_bedgraph} \
        | awk 'BEGIN {FS=OFS="\\t"} NF >= 4 {print \$1, \$2, \$3, \$4}' \
        > igv_track.bedGraph

    if [ ! -s igv_regions.bed ]; then
        grep -v "^track" ${methylation_bedgraph} \
            | awk 'BEGIN {FS=OFS="\\t"} NF >= 3 {print \$1, \$2, \$3}' \
            | head -n ${max_regions} > igv_regions.bed
    fi

    create_report igv_regions.bed ${params.referenceGenome}.fa \
        --flanking 1000 \
        --tracks igv_track.bedGraph ${bam} refgene.bed \
        --standalone \
        --output ${sampleId}_IGV.html
    """
}