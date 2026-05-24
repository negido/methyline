process dss {

    tag 'dss'
    label 'high'

    container 'quay.io/biocontainers/bioconductor-dss:2.58.0--r45h01b2380_0'

    publishDir "${params.analysisName}/bed", mode: 'copy', pattern: "*.bed"
    publishDir "${params.analysisName}/tsv", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.analysisName}/pdf", mode: 'copy', pattern: "*.pdf"

    input:
    path bedgraph_files
    path pmd_files
    path design_matrix

    output:
    tuple val("dmr"), path("*vs*_DMRs.bed"),    emit: dmr_bed
    tuple val("dmr"), path("consensus_DMRs.bed"), emit: consensus_bed
    tuple val("dmr"), path("*_DMRs_full.tsv"),  emit: dmr_full
    path "*_hyper.tsv",                         emit: hyper
    path "*_hypo.tsv",                          emit: hypo
    path "*_DMLtest.tsv",                       emit: smoothed
    path "*_DMR_stats.pdf",                     emit: dmr_stats
    path "cpg_qc_mqc.pdf",                      emit: qc_plot
    path "cpg_stats.tsv",                       emit: qc_stats

    script:
    """
    export DESIGN="${design_matrix}"

Rscript - <<'EOF'

    library(DSS)

    # ══════════════════════════════════════════════════════════════════════════
    # 1. Leer bedGraphs
    # ══════════════════════════════════════════════════════════════════════════

    files <- Sys.glob("*.bedGraph")

    if (length(files) < 2) stop("Se necesitan al menos 2 muestras.")

    read_bg <- function(f) {
        df  <- read.table(f, header = FALSE, sep = "", skip = 1,
                          fill = TRUE, quote = "")
        cov <- as.integer(df[[5]]) + as.integer(df[[6]])
        keep <- !is.na(cov) & cov >= 5
        df   <- df[keep, , drop = FALSE]
        cov  <- cov[keep]
        data.frame(
            chr = df[[1]],
            pos = as.integer(df[[3]]),
            N   = cov,
            X   = as.integer(df[[5]])
        )
    }

    sample_names <- sub(".bedGraph", "", basename(files), fixed = TRUE)

    dat_list <- lapply(seq_along(files), function(i) {
        d <- read_bg(files[i])
        if (nrow(d) == 0) stop(paste("Muestra sin CpGs con cov>=5:", files[i]))
        d
    })

    BSobj <- makeBSseqData(dat_list, sampleNames = sample_names)

    # ══════════════════════════════════════════════════════════════════════════
    # 2. Coverage filtering
    # ══════════════════════════════════════════════════════════════════════════

    keep  <- rowSums(getCoverage(BSobj, type = "Cov") >= 5) >= 2
    BSobj <- BSobj[keep, ]

    # ══════════════════════════════════════════════════════════════════════════
    # 3. Design matrix
    # ══════════════════════════════════════════════════════════════════════════

    design <- read.table(Sys.getenv("DESIGN"), header = TRUE, sep = "",
                         stringsAsFactors = FALSE)

    if (!"group" %in% colnames(design)) stop("design_matrix must contain column 'group'")
    design[["group"]] <- trimws(design[["group"]])
    groups <- unique(design[["group"]])
    message("Groups: ", paste(groups, collapse = ", "))
    if (length(groups) < 2) stop("At least 2 groups are required")

    # ── Detectar columna de nombres de muestra ────────────────────────────
    sample_col <- colnames(design)[which(sapply(colnames(design), function(col)
        any(colnames(BSobj) %in% design[[col]])))[1]]
    message("Columna de muestras detectada: ", sample_col)

    # ── Asignar grupo por muestra — con fallback robusto ──────────────────
    sample_group <- design[["group"]][match(colnames(BSobj), design[[sample_col]])]
    if (any(is.na(sample_group))) {
        warning("NAs en sample_group — asignando grupos en orden")
        sample_group <- rep(groups, each = ceiling(ncol(BSobj) / length(groups)))[
            seq_len(ncol(BSobj))
        ]
    }

    # ══════════════════════════════════════════════════════════════════════════
    # 4. QC PLOTS
    # ══════════════════════════════════════════════════════════════════════════

    cov_matrix  <- getCoverage(BSobj, type = "Cov")
    meth_matrix <- getMeth(BSobj, type = "raw")

    # ── Paleta de colores por grupo ───────────────────────────────────────
    group_base_cols <- c("#4f98a3", "#e07b39", "#7b5ea7", "#6aab4f",
                         "#c94040", "#d4a017", "#3a7ebf", "#888888")
    group_colors <- setNames(group_base_cols[seq_along(groups)], groups)
    sample_colors <- as.character(group_colors[sample_group])

    # ── Métricas por muestra ──────────────────────────────────────────────
    n_covered_per_sample    <- colSums(cov_matrix > 0)
    mean_depth_per_sample   <- colMeans(cov_matrix, na.rm = TRUE)
    median_depth_per_sample <- apply(cov_matrix, 2, median, na.rm = TRUE)
    global_meth_per_sample  <- colMeans(meth_matrix * 100, na.rm = TRUE)

    # ── Solapamiento ──────────────────────────────────────────────────────
    n_samples_per_cpg <- rowSums(cov_matrix > 0)
    n_shared_all      <- sum(n_samples_per_cpg == ncol(BSobj))
    overlap_counts    <- table(factor(n_samples_per_cpg,
                                      levels = seq_len(ncol(BSobj))))

    # ── CpGs compartidos para boxplot ─────────────────────────────────────
    shared_mask <- n_samples_per_cpg == ncol(BSobj)
    meth_shared <- meth_matrix[shared_mask, ] * 100

    # ── Colores y lty ECDF: por grupo + variación de tono intragrupo ──────
    ecdf_colors <- character(ncol(BSobj))
    ecdf_lty    <- integer(ncol(BSobj))
    group_lty   <- setNames(seq_along(groups), groups)

    for (g in groups) {
        idx      <- which(sample_group == g)
        base_col <- group_colors[[g]]
        n_in     <- length(idx)
        if (n_in == 0) next
        shades <- if (n_in == 1) {
            base_col
        } else {
            colorRampPalette(c(
                adjustcolor(base_col, red.f = 1.3, green.f = 1.3, blue.f = 1.3),
                base_col,
                adjustcolor(base_col, red.f = 0.6, green.f = 0.6, blue.f = 0.6)
            ))(n_in)
        }
        ecdf_colors[idx] <- shades
        ecdf_lty[idx]    <- as.integer(group_lty[[g]])
    }

    pdf("cpg_qc_mqc.pdf", width = 14, height = 10)
    par(mfrow = c(2, 3), mar = c(6, 5, 4, 2))

    # ── Plot 1: nº CpGs cubiertos por muestra ─────────────────────────────
    barplot(as.numeric(n_covered_per_sample) / 1e6,
            names.arg = colnames(BSobj),
            col       = sample_colors,
            border    = "white",
            main      = "CpGs covered per sample",
            ylab      = "Million CpGs",
            las       = 2,
            cex.names = 0.7)
    legend("topright",
           legend = names(group_colors),
           fill   = as.character(group_colors),
           bty    = "n", cex = 0.8)

    # ── Plot 2: average read depth por muestra ────────────────────────────
    barplot(as.numeric(mean_depth_per_sample),
            names.arg = colnames(BSobj),
            col       = sample_colors,
            border    = "white",
            main      = "Average read depth per sample",
            ylab      = "Mean depth",
            las       = 2,
            cex.names = 0.7)
    abline(h = mean(mean_depth_per_sample), col = "red", lty = 2, lwd = 1.5)
    text(x      = 0.5,
         y      = mean(mean_depth_per_sample) * 1.05,
         labels = paste0("mean = ", round(mean(mean_depth_per_sample), 1)),
         col    = "red", cex = 0.7, adj = 0)

    # ── Plot 3: solapamiento entre muestras ───────────────────────────────
    barplot(as.numeric(overlap_counts),
            names.arg = seq_len(ncol(BSobj)),
            col       = "#4f98a3",
            border    = "white",
            main      = "CpG overlap across samples",
            xlab      = "Nº samples with coverage",
            ylab      = "Nº CpGs",
            sub       = paste0("Shared all (", ncol(BSobj), " samples): ",
                               format(n_shared_all, big.mark = ",")))

    # ── Plot 4: ECDF de cobertura — colores y lty distinguibles ───────────
    thresholds <- c(1, 2, 3, 5, 8, 10, 15, 20, 30, 50)

    plot(NA,
         xlim = range(thresholds), ylim = c(0, 100),
         main = "CpG cumulative coverage",
         xlab = "Minimum read depth",
         ylab = "% CpGs covered",
         xaxt = "n")
    axis(1, at = thresholds)
    abline(h = c(25, 50, 75), col = "grey85", lty = 2)

    for (i in seq_len(ncol(BSobj))) {
        pct <- sapply(thresholds, function(t)
            mean(cov_matrix[, i] >= t, na.rm = TRUE) * 100)
        lines(thresholds, pct,
              col = ecdf_colors[i],
              lty = ecdf_lty[i],
              lwd = 2)
        points(thresholds, pct,
               col = ecdf_colors[i],
               pch = 19, cex = 0.6)
    }
    legend("topright",
           legend = c(colnames(BSobj),
                      paste0(names(group_lty), " (lty=", group_lty, ")")),
           col    = c(ecdf_colors, as.character(group_colors)),
           lty    = c(ecdf_lty,    as.integer(group_lty)),
           lwd    = 2, cex = 0.55, bty = "n")

    # ── Plot 5: distribución % metilación — boxplot por muestra ──────────
    boxplot(as.data.frame(meth_shared),
            col      = sample_colors,
            border   = "#444444",
            outline  = FALSE,
            main     = paste0("Methylation distribution\n(shared CpGs, n = ",
                              format(sum(shared_mask), big.mark = ","), ")"),
            ylab     = "% methylation",
            las      = 2,
            cex.axis = 0.7)

    # ── Plot 6: global methylation por muestra ────────────────────────────
    barplot(as.numeric(global_meth_per_sample),
            names.arg = colnames(BSobj),
            col       = sample_colors,
            border    = "white",
            ylim      = c(0, 100),
            main      = "Global CpG methylation per sample",
            ylab      = "Mean % methylation",
            las       = 2,
            cex.names = 0.7)

    dev.off()

    # ── Stats TSV ─────────────────────────────────────────────────────────
    stats <- data.frame(
        sample        = colnames(BSobj),
        group         = sample_group,
        n_CpG_covered = as.integer(n_covered_per_sample),
        n_CpG_shared  = n_shared_all,
        mean_depth    = round(mean_depth_per_sample, 2),
        median_depth  = round(median_depth_per_sample, 2),
        global_meth   = round(global_meth_per_sample, 2)
    )
    write.table(stats, "cpg_stats.tsv",
                sep = "\t", quote = FALSE, row.names = FALSE)

    # ══════════════════════════════════════════════════════════════════════════
    # 5. Pairwise comparisons
    # ══════════════════════════════════════════════════════════════════════════

    comparisons <- combn(groups, 2, simplify = FALSE)

    for (cmp in comparisons) {

        group1          <- cmp[1]
        group2          <- cmp[2]
        comparison_name <- paste0(group2, "_vs_", group1)
        message("Running comparison: ", comparison_name)

        samples1 <- colnames(BSobj)[sample_group == group1]
        samples2 <- colnames(BSobj)[sample_group == group2]

        # ── DML test ──────────────────────────────────────────────────────
        DMLtest_res <- DMLtest(BSobj,
                               group1    = samples1,
                               group2    = samples2,
                               smoothing = TRUE)

        write.table(DMLtest_res,
                    paste0(comparison_name, "_DMLtest.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)

        # ── Call DMRs ─────────────────────────────────────────────────────
        DMRs <- callDMR(DMLtest_res,
                        delta       = 0.08,
                        p.threshold = 0.01,
                        minCG       = 3,
                        minlen      = 50,
                        dis.merge   = 100,
                        pct.sig     = 0.5)

        if (!is.null(DMRs) && nrow(DMRs) > 0) {

            # BED
            write.table(DMRs[, c("chr", "start", "end")],
                        paste0(comparison_name, "_DMRs.bed"),
                        sep = "\t", quote = FALSE,
                        row.names = FALSE, col.names = FALSE)

            # Hyper / Hypo
            DMRs[["meth_diff"]] <- if (all(c("meanMethy1","meanMethy2") %in% colnames(DMRs))) {
                DMRs[["meanMethy2"]] - DMRs[["meanMethy1"]]
            } else {
                DMRs[["areaStat"]]
            }

            write.table(DMRs,
                        paste0(comparison_name, "_DMRs_full.tsv"),
                        sep = "\t", quote = FALSE, row.names = FALSE)
            write.table(DMRs[DMRs[["meth_diff"]] > 0, ],
                        paste0(comparison_name, "_hyper.tsv"),
                        sep = "\t", quote = FALSE, row.names = FALSE)
            write.table(DMRs[DMRs[["meth_diff"]] < 0, ],
                        paste0(comparison_name, "_hypo.tsv"),
                        sep = "\t", quote = FALSE, row.names = FALSE)

            # ── DMR stats plot ─────────────────────────────────────────────
            pdf(paste0(comparison_name, "_DMR_stats.pdf"))
            par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))

            pvals <- if (!is.null(DMLtest_res[["pvals"]])) {
                DMLtest_res[["pvals"]]
            } else {
                DMLtest_res[["pval"]]
            }
            hist(na.omit(as.numeric(pvals)),
                 col = "#4f98a3", border = "white",
                 main = paste0(comparison_name, " - p-value distribution"),
                 xlab = "p-value")

            hist(DMRs[["length"]],
                 col = "#4f98a3", border = "white",
                 main = "DMR length distribution",
                 xlab = "Length (bp)")

            hist(DMRs[["nCG"]],
                 col = "#4f98a3", border = "white",
                 main = "CpGs per DMR",
                 xlab = "Nº CpGs")

            counts <- c(Hyper = sum(DMRs[["meth_diff"]] > 0),
                        Hypo  = sum(DMRs[["meth_diff"]] < 0))
            barplot(counts,
                    col = c("#e07b39", "#4f98a3"), border = "white",
                    main = paste0("DMRs: ", comparison_name),
                    ylab = "Nº DMRs")

            dev.off()

        } else {
            message("No DMRs detected for ", comparison_name)
        }
    }

    # ══════════════════════════════════════════════════════════════════════════
    # 6. Consensus merged BED
    # ══════════════════════════════════════════════════════════════════════════

    bed_files <- Sys.glob("*_DMRs.bed")

    if (length(bed_files) > 0) {

        all_beds <- do.call(rbind, lapply(bed_files, function(f) {
            df <- read.table(f, header = FALSE, sep = "\t",
                             stringsAsFactors = FALSE)
            colnames(df) <- c("chr", "start", "end")
            df
        }))

        all_beds <- all_beds[order(all_beds[["chr"]], all_beds[["start"]]), ]
        merged   <- all_beds[1, , drop = FALSE]

        for (i in 2:nrow(all_beds)) {
            current  <- all_beds[i, ]
            last_idx <- nrow(merged)
            if (current[["chr"]] == merged[last_idx, "chr"] &&
                current[["start"]] <= merged[last_idx, "end"]) {
                merged[last_idx, "end"] <- max(merged[last_idx, "end"],
                                               current[["end"]])
            } else {
                merged <- rbind(merged, current)
            }
        }

        write.table(merged, "consensus_DMRs.bed",
                    sep = "\t", quote = FALSE,
                    row.names = FALSE, col.names = FALSE)
    }

EOF
    """
}
