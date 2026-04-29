process annotatr {
	tag "$sampleId"
	label 'medium'
    container 'quay.io/biocontainers/bioconductor-annotatr:1.36.0--r45hdfd78af_0'

	input:
	tuple val(sampleId), path(dmr_bed)
    val referenceGenome

    output:
    tuple val(sampleId), path("${sampleId}_annotated.tsv"), emit: annotations
    tuple val(sampleId), path("${sampleId}_annotation_summary.tsv"), emit: summary
    tuple val(sampleId), path("${sampleId}_annotation_bar.pdf"), emit: bar_plot


    script:
    def txdb_pkg = (referenceGenome == "hg38")
        ? "TxDb.Hsapiens.UCSC.hg38.knownGene"
        : "TxDb.Hsapiens.UCSC.hg19.knownGene"
    """
    export TXDB_PKG="${txdb_pkg}"
    export GENOME="${referenceGenome}"
    export SID="${sampleId}"

    Rscript - <<'EOF'
    BiocManager::install(c(Sys.getenv("TXDB_PKG"), "org.Hs.eg.db"), ask = FALSE, update = FALSE)

    library(annotatr)
    library(GenomicRanges)

    genome <- Sys.getenv("GENOME")
    sid    <- Sys.getenv("SID")

    annots <- build_annotations(
        genome      = genome,
        annotations = c(
            paste0(genome, "_cpg_islands"),
            paste0(genome, "_cpg_shores"),
            paste0(genome, "_cpg_shelves"),
            paste0(genome, "_cpg_inter"),
            paste0(genome, "_genes_promoters"),
            paste0(genome, "_genes_exons"),
            paste0(genome, "_genes_introns"),
            paste0(genome, "_genes_intergenic")
        )
    )

    regions <- read_regions(
        con          = "${dmr_bed}",
        genome       = genome,
        format       = "bed",
        rename_score = "mean_meth_diff"
    )

    if (length(regions) == 0) {
        write.table(data.frame(), file = paste0(sid, "_annotated.tsv"),
                    sep = "\\t", quote = FALSE, row.names = FALSE)
        write.table(data.frame(), file = paste0(sid, "_annotation_summary.tsv"),
                    sep = "\\t", quote = FALSE, row.names = FALSE)
        pdf(paste0(sid, "_annotation_bar.pdf"))
        plot.new(); text(0.5, 0.5, "No regions to annotate")
        dev.off()
        quit(save = "no", status = 0)
    }

    annotated <- annotate_regions(
        regions       = regions,
        annotations   = annots,
        ignore.strand = TRUE,
        quiet         = FALSE
    )

    write.table(as.data.frame(annotated),
                file = paste0(sid, "_annotated.tsv"),
                sep = "\\t", quote = FALSE, row.names = FALSE)

    annot_summary <- summarize_annotations(annotated_regions = annotated, quiet = FALSE)
    write.table(as.data.frame(annot_summary),
                file = paste0(sid, "_annotation_summary.tsv"),
                sep = "\\t", quote = FALSE, row.names = FALSE)

    pdf(paste0(sid, "_annotation_bar.pdf"))
    plot_annotation(
        annotated_regions = annotated,
        annotation_order  = NULL,
        plot_title        = paste0("DMR annotation — ", sid),
        x_label           = "Annotation type",
        y_label           = "Count"
    )
    dev.off()
    EOF
    """
}