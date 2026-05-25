process annotatr {
    tag "$sampleId"
    label 'medium'
    container 'quay.io/biocontainers/bioconductor-annotatr:1.36.0--r45hdfd78af_0'

    publishDir "${params.analysisName}/tsv", mode: 'copy', pattern: "*_annotated.tsv"
    publishDir "${params.analysisName}/tsv", mode: 'copy', pattern: "*_annotation_summary.tsv"
    publishDir "${params.analysisName}/tsv", mode: 'copy', pattern: "*_annotation_enrichment.tsv"
    publishDir "${params.analysisName}/pdf", mode: 'copy', pattern: "*_mqc.pdf"

    input:
    tuple val(sampleId), path(dmr_tsv)
    val referenceGenome

    output:
    tuple val(sampleId), path("${sampleId}_annotated.tsv"),             emit: annotations
    tuple val(sampleId), path("${sampleId}_annotation_summary.tsv"),    emit: summary
    tuple val(sampleId), path("${sampleId}_annotation_enrichment.tsv"), emit: enrichment
    tuple val(sampleId), path("${sampleId}_annotation_bar_mqc.pdf"),    emit: bar_plot
    tuple val(sampleId), path("${sampleId}_annotation_enrich_mqc.pdf"), emit: enrich_plot

    script:
    def txdb_pkg = (referenceGenome == "hg38")
        ? "TxDb.Hsapiens.UCSC.hg38.knownGene"
        : "TxDb.Hsapiens.UCSC.hg19.knownGene"
    """
    export TXDB_PKG="${txdb_pkg}"
    export GENOME="${referenceGenome}"
    export SID="${sampleId}"
    export DMR_TSV="${dmr_tsv}"

Rscript - <<'EOF'
    BiocManager::install(c(Sys.getenv("TXDB_PKG"), "org.Hs.eg.db"),
                         ask = FALSE, update = FALSE)

    library(annotatr)
    library(GenomicRanges)

    genome <- Sys.getenv("GENOME")
    sid    <- Sys.getenv("SID")

    # ── Anotaciones ───────────────────────────────────────────────────────────
    cpg_annots  <- paste0(genome, c("_cpg_islands","_cpg_shores",
                                    "_cpg_shelves","_cpg_inter"))
    gene_annots <- paste0(genome, c("_genes_promoters","_genes_exons",
                                    "_genes_introns","_genes_intergenic"))
    annots <- build_annotations(genome      = genome,
                                annotations = c(cpg_annots, gene_annots))

    # ── Leer DMRs y construir BED ─────────────────────────────────────────────
    dm_data <- read.table(Sys.getenv("DMR_TSV"), header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
    tmp_bed <- tempfile(fileext = ".bed")
    write.table(data.frame(
        chr        = dm_data[["chr"]],
        start      = dm_data[["start"]],
        end        = dm_data[["end"]],
        name       = ifelse(dm_data[["meth_diff"]] > 0, "hyper", "hypo"),
        score      = dm_data[["areaStat"]],
        strand     = "*",
        meth_diff  = dm_data[["meth_diff"]],
        meanMethy1 = dm_data[["meanMethy1"]],
        meanMethy2 = dm_data[["meanMethy2"]]
    ), tmp_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

    regions <- tryCatch(
        read_regions(
            con          = tmp_bed,
            genome       = genome,
            extraCols    = c(meth_diff  = "numeric",
                             meanMethy1 = "numeric",
                             meanMethy2 = "numeric"),
            format       = "bed",
            rename_name  = "direction",
            rename_score = "areaStat"
        ),
        error = function(e) NULL
    )

    # ── Placeholder si no hay regiones ────────────────────────────────────────
    if (is.null(regions) || length(regions) == 0) {
        for (f in c("_annotated.tsv","_annotation_summary.tsv",
                    "_annotation_enrichment.tsv"))
            write.table(data.frame(), paste0(sid, f),
                        sep = "\t", quote = FALSE, row.names = FALSE)
        for (f in c("_annotation_bar_mqc.pdf","_annotation_enrich_mqc.pdf")) {
            pdf(paste0(sid, f))
            plot.new(); text(0.5, 0.5, "No regions to annotate")
            dev.off()
        }
        quit(save = "no", status = 0)
    }

    # ── Anotación y TSVs ──────────────────────────────────────────────────────
    annotated <- annotate_regions(regions = regions, annotations = annots,
                                  ignore.strand = TRUE, quiet = FALSE)
    annot_df  <- as.data.frame(annotated)

    write.table(annot_df,
                paste0(sid, "_annotated.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(as.data.frame(summarize_annotations(annotated, quiet = FALSE)),
                paste0(sid, "_annotation_summary.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    # ── Enriquecimiento (Fisher + BH) ─────────────────────────────────────────
    bg_df    <- as.data.frame(build_annotations(genome      = genome,
                                                annotations = c(cpg_annots,
                                                                gene_annots)))
    all_cats <- sort(unique(annot_df[["annot.type"]]))

    enrich_df <- do.call(rbind, lapply(all_cats, function(cat) {
        ft <- fisher.test(matrix(
            c(sum(annot_df[["annot.type"]] == cat),
              nrow(annot_df) - sum(annot_df[["annot.type"]] == cat),
              sum(bg_df[["type"]] == cat),
              nrow(bg_df)    - sum(bg_df[["type"]] == cat)),
            nrow = 2))
        data.frame(annotation = cat,
                   n_dmr      = sum(annot_df[["annot.type"]] == cat),
                   odds_ratio = as.numeric(ft[["estimate"]]),
                   pval       = ft[["p.value"]])
    }))
    enrich_df[["padj"]] <- p.adjust(enrich_df[["pval"]], method = "BH")
    enrich_df[["enriched"]] <- ifelse(
        enrich_df[["padj"]] < 0.05 & enrich_df[["odds_ratio"]] > 1, "enriched",
        ifelse(enrich_df[["padj"]] < 0.05 & enrich_df[["odds_ratio"]] < 1,
               "depleted", "not_significant"))

    write.table(enrich_df, paste0(sid, "_annotation_enrichment.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    # ══════════════════════════════════════════════════════════════════════════
    # PLOT 1 — Barplot annotatr
    # ══════════════════════════════════════════════════════════════════════════
    pdf(paste0(sid, "_annotation_bar_mqc.pdf"))
    plot_annotation(
        annotated_regions = annotated,
        plot_title        = paste0("DMR annotation — ", sid),
        x_label           = "Annotation type",
        y_label           = "Count"
    )
    dev.off()

    # ══════════════════════════════════════════════════════════════════════════
    # PLOT 2 — Lollipop horizontal con leyenda en margen inferior
    # ══════════════════════════════════════════════════════════════════════════
    ef       <- enrich_df[order(enrich_df[["odds_ratio"]]), ]
    log2or   <- log2(pmax(ef[["odds_ratio"]], 1e-9))
    cat_labs <- gsub(paste0(genome, "_"), "", ef[["annotation"]])
    pt_col   <- ifelse(ef[["enriched"]] == "enriched",  "#01696f",
                ifelse(ef[["enriched"]] == "depleted",  "#c94040", "#c0c0c0"))
    cex_pts  <- 1.2 + 2.5 * (ef[["n_dmr"]] / max(ef[["n_dmr"]]))
    padj_lab <- ifelse(ef[["padj"]] < 0.001,
                       formatC(ef[["padj"]], format = "e", digits = 2),
                       round(ef[["padj"]], 3))

    n      <- nrow(ef)
    y_pos  <- seq_len(n)
    xlim_l <- c(min(log2or) - 0.5, max(log2or) + 2.5)

    # oma bottom reserva espacio para la leyenda fuera del área de plot
    pdf(paste0(sid, "_annotation_enrich_mqc.pdf"), width = 9, height = 6)
    par(mar = c(2, 9, 4, 2), oma = c(6, 0, 0, 0))

    plot(NA, xlim = xlim_l, ylim = c(0.5, n + 0.5),
         xlab = expression(log[2]~"Odds Ratio"),
         ylab = "", yaxt = "n",
         main = paste0(sid, " — Annotation enrichment"))
    axis(2, at = y_pos, labels = cat_labs, las = 1, cex.axis = 0.85)
    abline(v = 0, col = "#aaaaaa", lty = 2)
    abline(h = y_pos, col = "#f0f0f0", lty = 1)

    segments(x0 = 0, x1 = log2or, y0 = y_pos, y1 = y_pos,
             col = "#bbbbbb", lwd = 1.5)
    points(log2or, y_pos,
           pch = 21, col = "#444444",
           bg  = pt_col, cex = cex_pts, lwd = 0.8)
    text(max(xlim_l) - 0.1, y_pos,
         labels = paste0("padj=", padj_lab),
         cex = 0.65, adj = 1, col = "#444444")

    # Leyenda en el margen inferior (fuera del área de plot, xpd = NA)
    legend(x      = grconvertX(0.5, from = "npc", to = "user"),
           y      = grconvertY(-0.18, from = "npc", to = "user"),
           legend = c("Enriched (padj<0.05)", "Depleted (padj<0.05)",
                      "Not significant"),
           pch    = 21, col = "#444444",
           pt.bg  = c("#01696f","#c94040","#c0c0c0"),
           horiz  = TRUE, bty = "n", cex = 0.8,
           xjust  = 0.5, xpd = NA)

    dev.off()

EOF
    """
}