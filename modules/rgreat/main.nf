process rgreat {
    tag "$sampleId"
    label 'medium'
    container  'docker.io/negido/rgreat-local:latest'
    input:
    tuple val(sampleId), path(dmr_bed)


    output:
    tuple val(sampleId), path("GREAT_GO_BP.tsv"),                    emit: go_bp
    tuple val(sampleId), path("GREAT_GO_CC.tsv"),                    emit: go_cc
    tuple val(sampleId), path("GREAT_GO_MF.tsv"),                    emit: go_mf
    tuple val(sampleId), path("GREAT_volcano.pdf"),                  emit: volcano
    tuple val(sampleId), path("GREAT_region_gene_associations.pdf"), emit: region_gene_plot
    tuple val(sampleId), path("GREAT_region_gene_associations.tsv"), emit: region_gene_tsv

    script:
    def tss_source = (params.referenceGenome == "hg19") ? "txdb:hg19" : "txdb:hg38"
    """
    Rscript - <<'EOF'
    library(rGREAT)
    library(rtracklayer)
    library(GenomicRanges)

    write_empty <- function(path) {
        write.table(
            data.frame(
                id                   = character(),
                description          = character(),
                genome_fraction      = numeric(),
                observed_region_hits = integer(),
                fold_enrichment      = numeric(),
                p_value              = numeric(),
                p_adjust             = numeric()
            ),
            file = path, sep = "\\t", quote = FALSE, row.names = FALSE
        )
    }

    dmr_gr <- rtracklayer::import.bed("${dmr_bed}")

    if (length(dmr_gr) == 0) {
        write_empty("GREAT_GO_BP.tsv")
        write_empty("GREAT_GO_CC.tsv")
        write_empty("GREAT_GO_MF.tsv")
        pdf("GREAT_volcano.pdf")
        plot.new(); text(0.5, 0.5, "No DMR regions to analyse")
        dev.off()
        pdf("GREAT_region_gene_associations.pdf")
        plot.new(); text(0.5, 0.5, "No DMR regions to analyse")
        dev.off()
        write_empty("GREAT_region_gene_associations.tsv")
        quit(save = "no", status = 0)
    }

    tss_source <- "${tss_source}"

    for (ont in c("GO:BP", "GO:CC", "GO:MF")) {
        res  <- great(dmr_gr, ont, tss_source)
        tb   <- getEnrichmentTable(res)
        safe <- gsub(":", "_", ont)
        write.table(tb, file = paste0("GREAT_", safe, ".tsv"),
                    sep = "\\t", quote = FALSE, row.names = FALSE)
    }

    res_bp <- great(dmr_gr, "GO:BP", tss_source)

    pdf("GREAT_volcano.pdf")
    plotVolcano(res_bp)
    dev.off()

    pdf("GREAT_region_gene_associations.pdf")
    plotRegionGeneAssociations(res_bp)
    dev.off()

    assoc    <- getRegionGeneAssociations(res_bp)
    n        <- length(assoc)
    assoc_df <- data.frame(
        chr             = as.character(seqnames(assoc)),
        start           = start(assoc),
        end             = end(assoc),
        annotated_genes = sapply(seq_len(n), function(i)
                              paste(as.character(mcols(assoc)\$annotated_genes[[i]]), collapse = ";")),
        dist_to_TSS     = sapply(seq_len(n), function(i)
                              paste(as.character(mcols(assoc)\$dist_to_TSS[[i]]), collapse = ";")),
        stringsAsFactors = FALSE
    )
    write.table(assoc_df, file = "GREAT_region_gene_associations.tsv",
                sep = "\\t", quote = FALSE, row.names = FALSE)

    cat("rGREAT local analysis completed\\n")
    EOF
    """
}