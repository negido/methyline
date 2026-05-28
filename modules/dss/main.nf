process dss {

    tag 'dss'
    label 'high'

    container 'quay.io/biocontainers/bioconductor-dss:2.58.0--r45h01b2380_0'

    publishDir "${params.analysisName}/bed", mode: 'copy', pattern: "*.bed"
    publishDir "${params.analysisName}/tsv", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.analysisName}/pdf", mode: 'copy', pattern: "*.pdf"

    input:
    path bedgraph_files
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
        meth_pct <- as.numeric(df[[4]])
        cov <- as.integer(df[[5]]) + as.integer(df[[6]])
        keep <- !is.na(cov) & cov >= 1
        df   <- df[keep, , drop = FALSE]
        cov  <- cov[keep]
        data.frame(
            chr = df[[1]],
            start = as.integer(df[[2]]),
            pos = as.integer(df[[3]]),
            end = as.integer(df[[3]]),
            N   = cov,
            X   = as.integer(round(as.numeric(df[[4]]) / 100 * cov)),
            meth_pct = meth_pct[keep],
            coord = paste(df[[1]], df[[2]], df[[3]], sep=":"),
            stringsAsFactors = FALSE
        )
    }

    sample_names <- sub(".bedGraph", "", basename(files), fixed = TRUE)

    dat_list <- lapply(seq_along(files), function(i) {
        d <- read_bg(files[i])
        if (nrow(d) == 0) stop(paste("Muestra sin CpGs", files[i]))
        d
    })
    BSobj <- makeBSseqData(dat_list, sampleNames = sample_names)

    # ══════════════════════════════════════════════════════════════════════════
    # 2. Design matrix
    # ══════════════════════════════════════════════════════════════════════════

    design <- read.table(Sys.getenv("DESIGN"), header = TRUE, sep = "",
                        stringsAsFactors = FALSE)

    if (!"group" %in% colnames(design)) {
        stop("design_matrix must contain column 'group'")
    }

    design[["group"]] <- trimws(design[["group"]])

    groups <- unique(design[["group"]])

    message("Groups: ", paste(groups, collapse = ", "))

    if (length(groups) < 2) {
        stop("At least 2 groups are required")
    }

    # ── Normalizar nombres de muestras ─────────────────────────────────────
    bs_samples <- sub("_.*", "", colnames(BSobj))

    # ── Detectar automáticamente columna compatible ────────────────────────
    sample_col <- colnames(design)[
        which(
            sapply(colnames(design), function(col)
                any(bs_samples %in% design[[col]])
            )
        )[1]
    ]

    message("Columna de muestras detectada: ", sample_col)

    if (is.na(sample_col)) {
        stop("No sample column matches BSobj sample names")
    }

    design_sample_names <- design[[sample_col]]
    plot_order <- match(design_sample_names, bs_samples)

    if (any(is.na(plot_order))) {
        missing_samples <- design_sample_names[is.na(plot_order)]
        stop(
            paste0(
                "No sample in BSobj matches design matrix order: ",
                paste(missing_samples, collapse = ", ")
            )
        )
    }

    # ── Asignar grupo por muestra ──────────────────────────────────────────
    sample_group <- design[["group"]][
        match(bs_samples, design[[sample_col]])
    ]

    if (any(is.na(sample_group))) {

        missing_samples <- bs_samples[is.na(sample_group)]

        stop(
            paste0(
                "No matching group for samples: ",
                paste(missing_samples, collapse = ", ")
            )
        )
    }
    # ══════════════════════════════════════════════════════════════════════════
    # 3. QC PLOTS
    # ══════════════════════════════════════════════════════════════════════════

    # ── Asegurar orden consistente: usar el orden de la design matrix ────────
    dat_list_ordered <- dat_list[plot_order]
    names(dat_list_ordered) <- design_sample_names
    sample_group_ordered <- design[["group"]]

    # ── Paleta de colores por grupo ───────────────────────────────────────────
    group_base_cols <- c("#4f98a3", "#e07b39", "#7b5ea7", "#6aab4f",
                         "#c94040", "#d4a017", "#3a7ebf", "#888888")
    group_colors  <- setNames(group_base_cols[seq_along(groups)], groups)
    sample_colors <- as.character(group_colors[sample_group_ordered])

    # ── Métricas por muestra — calculadas desde dat_list_ordered (raw) ────────
    n_covered_per_sample    <- sapply(dat_list_ordered, nrow)
    mean_depth_per_sample   <- sapply(dat_list_ordered,
                                      function(x) mean(x[["N"]], na.rm = TRUE))
    median_depth_per_sample <- sapply(dat_list_ordered,
                                      function(x) median(x[["N"]], na.rm = TRUE))
    # Metilación global: X/N por CpG, luego media — calculado sobre raw
    global_meth_per_sample  <- sapply(dat_list_ordered,
                                      function(x) mean(x[["meth_pct"]],
                                                       na.rm = TRUE))

    # ── Solapamiento entre muestras (sobre dat_list_ordered, sin BSobj) ───────
    # Clave chr:start:end para cada muestra (sin duplicados internos)
    key_lists <- lapply(dat_list_ordered, function(x)
        unique(paste(as.character(x[["chr"]]), x[["start"]], x[["end"]], sep = ":")))

    # Union de todos los CpGs únicos
    all_keys <- unique(unlist(key_lists))
    message("[DSS] Total CpGs union (cov>=1): ", format(length(all_keys), big.mark=","))

    # Para cada CpG: en cuántas muestras aparece
    n_samples_per_cpg <- Reduce("+", lapply(key_lists, function(k)
        as.integer(all_keys %in% k)))

    repeat_levels <- 2:min(9, length(dat_list_ordered))
    repeat_counts <- vapply(repeat_levels, function(k)
        sum(n_samples_per_cpg >= k), integer(1))
    n_shared_all   <- sum(n_samples_per_cpg == length(dat_list_ordered))
    message("[DSS] CpGs compartidos (todas las muestras): ",
            format(n_shared_all, big.mark = ","))
        message("[DSS] Repetidos acumulados por nº de muestras: ",
            paste(sprintf("%s=%s", repeat_levels, format(repeat_counts, big.mark = ",")), collapse = " | "))

        # ── CpGs compartidos para boxplot — calculado desde dat_list_ordered ──────
        shared_keys <- all_keys[n_samples_per_cpg == length(dat_list_ordered)]
        message("[DSS] CpGs compartidos para boxplot: ",
            format(length(shared_keys), big.mark = ","))

    # Para cada muestra: extraer metilación en porcentaje en los CpGs compartidos
    meth_shared_list <- lapply(dat_list_ordered, function(x) {
        key <- paste(as.character(x[["chr"]]), x[["start"]], x[["end"]], sep = ":")
        idx <- match(shared_keys, key)          # NA si no está (no debería pasar)
        pct <- x[["meth_pct"]][idx]
        pct[is.finite(pct)]                     # eliminar NaN/Inf residuales
    })
    # Alinear longitudes para el boxplot (todos deberían tener misma longitud)
    meth_shared_df <- as.data.frame(do.call(cbind, lapply(
        seq_along(meth_shared_list), function(i) {
            v <- meth_shared_list[[i]]
            length(v) <- length(shared_keys)    # padding con NA si hay mismatch
            v
        }
    )))
    colnames(meth_shared_df) <- design_sample_names
    message("[DSS] Rango plot 5 por muestra (min/median/max):")
    for (i in seq_along(dat_list_ordered)) {
        vals <- meth_shared_df[[i]]
        message("[DSS]   ", design_sample_names[i], ": ",
                paste(round(range(vals, na.rm = TRUE), 1), collapse = "-"),
                " | median=", round(median(vals, na.rm = TRUE), 1),
                " | mean=", round(mean(vals, na.rm = TRUE), 1),
                " | n=", sum(!is.na(vals)))
    }
    message("[DSS] Primeras filas plot 5:")
    print(utils::head(meth_shared_df, 3))

    # ── Colores y lty ECDF: por grupo + variación de tono intragrupo ──────────
    ecdf_colors <- character(ncol(BSobj))
    ecdf_lty    <- integer(ncol(BSobj))
    group_lty   <- setNames(seq_along(groups), groups)

    for (g in groups) {
        idx      <- which(sample_group_ordered == g)
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

    # ── Plot 1: nº CpGs cubiertos por muestra ─────────────────────────────────
    message("[DSS] Plot 1/6: CpGs cubiertos por muestra")
    barplot(as.numeric(n_covered_per_sample) / 1e6,
            ylim      = c(0, max(n_covered_per_sample / 1e6) * 1.25),
            names.arg = design_sample_names,
            col       = sample_colors,
            border    = "white",
            main      = "CpGs covered per sample (cov\u22651)",
            ylab      = "Million CpGs",
            las       = 2,
            cex.names = 0.7)
    legend("topright",
           legend = names(group_colors),
           fill   = as.character(group_colors),
           bty    = "n", cex = 0.8)

    # ── Plot 2: average read depth por muestra ────────────────────────────────
        message("[DSS] Plot 2/6: profundidad media por muestra")
    barplot(as.numeric(mean_depth_per_sample),
            ylim      = c(0, max(mean_depth_per_sample) * 1.25),
            names.arg = design_sample_names,
            col       = sample_colors,
            border    = "white",
            main      = "Average read depth per sample (cov\u22651)",
            ylab      = "Mean depth",
            las       = 2,
            cex.names = 0.7)
    abline(h = mean(mean_depth_per_sample), col = "red", lty = 2, lwd = 1.5)
    text(x      = 0.5,
         y      = mean(mean_depth_per_sample) * 1.08,
         labels = paste0("global mean = ", round(mean(mean_depth_per_sample), 1)),
         col    = "red", cex = 0.7, adj = 0)
    legend("topright",
           legend = names(group_colors),
           fill   = as.character(group_colors),
           bty    = "n", cex = 0.8)

    # ── Plot 3: solapamiento entre muestras ───────────────────────────────────
    message("[DSS] Plot 3/6: CpGs repetidos en 2-9 muestras")
    barplot(as.numeric(repeat_counts),
            names.arg = repeat_levels,
            col       = "#4f98a3",
            border    = "white",
            main      = "CpGs repeated across samples (cov\u22651)",
            xlab      = "N\u00ba samples",
            ylab      = "N\u00ba CpGs",
            sub       = paste0("Shared all (", length(dat_list_ordered),
                               " samples): ", format(n_shared_all, big.mark = ",")))

    # ── Plot 4: ECDF de cobertura ─────────────────────────────────────────────
    message("[DSS] Plot 4/6: ECDF de cobertura")
    thresholds <- c(1, 2, 3, 5, 8, 10, 15, 20, 30, 50)

    ecdf_pct <- lapply(dat_list_ordered, function(x)
        sapply(thresholds, function(t) mean(x[["N"]] >= t, na.rm = TRUE) * 100))

    plot(NA,
         xlim = range(thresholds), ylim = c(0, 100),
         main = "CpG cumulative coverage (cov\u22651)",
         xlab = "Minimum read depth",
         ylab = "% CpGs covered",
         xaxt = "n")
    axis(1, at = thresholds)
    abline(h = c(25, 50, 75), col = "grey85", lty = 2)
    abline(v = 5, col = "darkred", lty = 2, lwd = 1.5)
    text(x = 5.5, y = 98, labels = "cov=5", col = "darkred", cex = 0.65, adj = 0)

    for (i in seq_along(dat_list_ordered)) {
        lines(thresholds, ecdf_pct[[i]],
              col = ecdf_colors[i], lty = ecdf_lty[i], lwd = 2)
        points(thresholds, ecdf_pct[[i]],
               col = ecdf_colors[i], pch = 19, cex = 0.6)
    }
    legend("topright",
           legend = c(colnames(BSobj),
                      paste0(names(group_lty), " (lty=", group_lty, ")")),
           col    = c(ecdf_colors, as.character(group_colors)),
           lty    = c(ecdf_lty, as.integer(group_lty)),
           lwd    = 2, cex = 0.55, bty = "n")

    # ── Plot 5: distribución % metilación (CpGs compartidos, raw) ────────────
        message("[DSS] Plot 5/6: distribución de metilación")
    boxplot(meth_shared_df,
            col      = sample_colors,
            border   = "#444444",
            outline  = FALSE,
            main     = paste0("Methylation distribution\n(shared CpGs, n = ",
                              format(length(shared_keys), big.mark = ","), ")"),
            ylab     = "% methylation",
            ylim     = c(0, 100),
            las      = 2,
            cex.axis = 0.7)

    # ── Plot 6: metilación global por muestra (raw) ───────────────────────────
            message("[DSS] Plot 6/6: metilación global por muestra")
    barplot(as.numeric(global_meth_per_sample),
            names.arg = design_sample_names,
            col       = sample_colors,
            border    = "white",
            ylim      = c(0, 100),
            main      = "Global CpG methylation per sample (cov\u22651)",
            ylab      = "Mean % methylation",
            las       = 2,
            cex.names = 0.7)
    abline(h = mean(global_meth_per_sample), col = "red", lty = 2, lwd = 1.5)
    text(x      = 0.5,
         y      = mean(global_meth_per_sample) + 3,
         labels = paste0("global mean = ", round(mean(global_meth_per_sample), 1), "%"),
         col    = "red", cex = 0.7, adj = 0)
    legend("topright",
           legend = names(group_colors),
           fill   = as.character(group_colors),
           bty    = "n", cex = 0.8)

    dev.off()

    # ── Stats TSV ─────────────────────────────────────────────────────────────
    stats <- data.frame(
        sample        = design_sample_names,
        group         = sample_group_ordered,
        n_CpG_covered = as.integer(n_covered_per_sample),
        n_CpG_shared  = n_shared_all,
        mean_depth    = round(mean_depth_per_sample, 2),
        median_depth  = round(median_depth_per_sample, 2),
        global_meth   = round(global_meth_per_sample, 2)
    )
    write.table(stats, "cpg_stats.tsv",
                sep = "\t", quote = FALSE, row.names = FALSE)
    # ══════════════════════════════════════════════════════════════════════════
    # 4. Pairwise comparisons
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
                               smoothing = TRUE,
                               ncores = ${task.cpus})

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
    # 5. Consensus merged BED
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
