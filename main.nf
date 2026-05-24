nextflow.enable.dsl=2

params.reads           = params.reads ?: 'test-data/*_{1,2}.fastq.gz'
params.genome_fasta  = params.referenceGenome == 'hg38' ? '/data/ref/hg38.fa*' : '/data/ref/hg19.fa*'
params.design_matrix   = params.design_matrix ?: params.design ?: ''
params.referenceGenome = params.referenceGenome ?: 'hg19'
params.refGene         = params.referenceGenome == 'hg38' ? '/data/ref/refGene_hg38.txt.gz' : '/data/ref/refGene_hg19.txt.gz'
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
    methylseekr(methylDackel.out.bedgraph)


    // Rama comparativa vs single-sample
    if (!params.single_sample) {
        // DMR con DSS (modo comparativo)
        dss(
            methylDackel.out.bedgraph
                .map { sampleId, bedgraph -> bedgraph }
                .collect(),
            methylseekr.out.pmds
                .map { sampleId, pmd -> pmd }
                .collect(),
            params.design_matrix ? file(params.design_matrix) : file(params.empty_design)
        )
        pairwise_dmrs = dss.out.dmr_full
        .flatMap { tag, beds ->
             beds.collect { bed ->
                 tuple(bed.baseName, bed)
            }
        }
        pairwise_dmrs_rgreat = dss.out.dmr_bed
        .flatMap { tag, beds ->
             beds.collect { bed ->
                 tuple(bed.baseName, bed)
            }
        }
            // Anotacion funcional de DMR
        annotatr(pairwise_dmrs,params.referenceGenome)
        rgreat(pairwise_dmrs_rgreat)

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

        // Visualizacion IGV con DMRs como regiones de interés
        igv(
            markDuplicates.out.MARKED_DUP
                .join(methylDackel.out.bedgraph),
            file(params.genome_fasta),
            params.referenceGenome,
            file(params.refGene),
            dss.out.consensus_bed.ifEmpty { file("NO_FILE") }
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
                .mix(methylseekr.out.pmds)
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

        // Visualizacion IGV sin DMRs (muestreo de bedGraph)
        igv(
            markDuplicates.out.MARKED_DUP
                .join(methylDackel.out.bedgraph),
            file(params.genome_fasta),
            params.referenceGenome,
            file(params.refGene),
            file("NO_FILE")
        )
    }
}
