nextflow.enable.dsl=2

params.reads           = params.reads ?: 'test-data/*_{1,2}.fastq.gz'
params.genome_fasta  = params.referenceGenome == 'hg38' ? 'ref/hg38.fa*' : 'ref/hg19.fa*'
params.design_matrix   = params.design_matrix ?: params.design ?: ''
params.referenceGenome = params.referenceGenome ?: 'hg19'
params.refGene         = params.referenceGenome == 'hg38' ? 'ref/refGene_hg38.txt.gz' : 'ref/refGene_hg19.txt.gz'
params.empty_design    = params.empty_design ?: 'test-data/empty_design_matrix.tsv'
params.single_sample   = params.single_sample ?: false
params.analysisName    = params.analysisName ?: ''

include { fastqc }      from './modules/fastqc/main'
include { trimgalore }  from './modules/trimgalore/main'
include { multiQC }     from './modules/multiqc/main'
include { bwamethAlign } from './modules/bwameth/main'
include { markDuplicates } from './modules/markduplicates/main'
include { methylDackel } from './modules/methyldackel/main'
include { methylseekr }  from './modules/methylseekr/main'
include { dss }          from './modules/dss/main'
include { annotatr }     from './modules/annotatr/main'
include { rgreat }       from './modules/rgreat/main'
include { igv }          from './modules/igv/main'

workflow {
    // Entrada de lecturas
    reads = Channel.fromFilePairs(params.reads, checkIfExists: true)

    // Preprocesamiento y control de calidad
    fastqc(reads)
    trimgalore(reads)

    // Alineamiento y marcado de duplicados
    bwamethAlign(trimgalore.out.reads, file(params.genome_fasta), params.referenceGenome)
    markDuplicates(bwamethAlign.out.bam, bwamethAlign.out.bai)

    // Extraccion de metilacion y segmentacion de regiones
    methylDackel(markDuplicates.out.MARKED_DUP, file(params.genome_fasta))
    methylseekr(methylDackel.out.bedgraph, file(params.genome_fasta))

    // Rama comparativa vs single-sample
    if (!params.single_sample) {
        // DMR con DSS (modo comparativo)
        dss(
            methylDackel.out.bedgraph
                .map { sampleId, bedgraph -> bedgraph }
                .collect(),
            params.design_matrix ? file(params.design_matrix) : file(params.empty_design),
            methylseekr.out.pmds
                .map { sampleId, pmd -> pmd }
                .collect()
        )

            // Anotacion funcional de DMR
        annotatr(dss.out.dmr_bed,params.referenceGenome)
        rgreat(dss.out.dmr_bed,params.referenceGenome)

            // Informe calidad
        multiQC(
            fastqc.out.html
                .mix(fastqc.out.zip)
                .mix(trimgalore.out.report)
                .mix(trimgalore.out.fastqc_html)
                .mix(trimgalore.out.fastqc_zip)
                .mix(markDuplicates.out.MD_METRICS)
                .mix(dss.out.status)
                .collect()
        )
    } else {
        // Anotacion funcional de regiones LMR/UMR/PMD
        annotatr(methylseekr.out.lmrs
                    .mix(methylseekr.out.umrs)
                    .mix(methylseekr.out.pmds),
                params.referenceGenome)
        rgreat(
            methylseekr.out.lmrs
                .mix(methylseekr.out.umrs)
                .mix(methylseekr.out.pmds),
            params.referenceGenome
        )

        // Informe calidad
        multiQC(
            fastqc.out.html
                .mix(fastqc.out.zip)
                .mix(trimgalore.out.report)
                .mix(trimgalore.out.fastqc_html)
                .mix(trimgalore.out.fastqc_zip)
                .mix(markDuplicates.out.MD_METRICS)
                .collect()
        )
    }

    // Visualizacion IGV
    igv(
        markDuplicates.out.MARKED_DUP
            .join(methylDackel.out.bedgraph),
        file(params.genome_fasta),
        params.referenceGenome,
        file(params.refGene)
    )
}
