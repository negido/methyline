process igv {
    tag "$sampleId"
    label 'high'
    container 'negido/igv-reports:v1.0'

    publishDir "${params.analysisName}/html", mode: 'copy', tags: [Type: 'html']

    input:
    tuple val(sampleId), path(bam), path(bai), path(methylation_bedgraph)
    path genome_fasta
    val referenceGenome
    path refGene
    tuple val(sample),path(dmr_bed)

    output:
    tuple val(sampleId), path("${sampleId}_IGV.html"), emit: IGV_HTML

    script:

    def max_regions = params.igv_max_regions ?: 500
    def has_dmrs = dmr_bed != null && dmr_bed.name != 'NO_FILE'

    """
    set -euo pipefail

    echo ">>> Preparando anotaciones"
    zcat ${refGene} > refgene.refgene

    echo ">>> Generando regiones IGV"

    # Si existen DMRs, usarlas
    if ${has_dmrs}; then
        cp ${dmr_bed} igv_regions.bed
    else
        grep -v "^track" ${methylation_bedgraph} \
            | shuf -n ${max_regions} \
            | sort -k1,1 -k2,2n \
            | awk 'BEGIN {FS=OFS="\\t"} NF >= 3 {print \$1, \$2, \$3}' \
            > igv_regions.bed
    fi

    # Fallback por si el BED queda vacío
    if [ ! -s igv_regions.bed ]; then
        grep -v "^track" ${methylation_bedgraph} \
            | awk 'BEGIN {FS=OFS="\\t"} NF >= 3 {print \$1, \$2, \$3}' \
            | head -n ${max_regions} \
            > igv_regions.bed
    fi

    echo ">>> Preparando señal de metilación"

    grep -v "^track" ${methylation_bedgraph} \
        | awk 'BEGIN {FS=OFS="\\t"} NF >= 4 {print \$1, \$2, \$3, \$4}' \
        | sort -k1,1 -k2,2n \
        > igv_track.sorted.bedGraph

    echo ">>> Detectando rango de metilación"

    max_val=\$(awk 'BEGIN{max=0} {if(\$4>max) max=\$4} END{print max}' igv_track.sorted.bedGraph)

    MIN_VAL=0

    awk -v m="\$max_val" 'BEGIN {
        if (m <= 1.0)
            print "MAX_VAL=1"
        else
            print "MAX_VAL=100"
    }' > scale.tmp

    source scale.tmp

    echo ">>> Creando track-config.json"

    cat <<EOF > track_config.json
[
  {
    "name": "${sampleId} methylation",
    "url": "igv_track.sorted.bedGraph",
    "format": "bedgraph",
    "type": "wig",
    "color": "blue",
    "autoscale": false,
    "min": \$MIN_VAL,
    "max": \$MAX_VAL
  },
  {
    "name": "${sampleId} alignments",
    "url": "${bam}",
    "indexURL": "${bai}",
    "format": "bam",
    "type": "alignment",
    "height": 120
  },
  {
    "name": "Genes",
    "url": "refgene.refgene",
    "format": "refgene",
    "type": "annotation"
  }
EOF

    # Añadir DMRs si existen
    if ${has_dmrs}; then
        cat <<EOF >> track_config.json
,
  {
    "name": "DMRs",
    "url": "igv_regions.bed",
    "format": "bed",
    "type": "annotation",
    "color": "purple"
  }
EOF
    fi

    echo "]" >> track_config.json

    echo ">>> Generando IGV report"

    create_report igv_regions.bed ${referenceGenome}.fa \
        --flanking 1000 \
        --track-config track_config.json \
        --standalone \
        --output ${sampleId}_IGV.html

    echo ">>> DONE"
    """
}
