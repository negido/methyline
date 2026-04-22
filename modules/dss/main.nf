process dss {
    tag "dss"
    label 'high'
    container 'quay.io/biocontainers/bioconductor-dss:2.58.0--r45h01b2380_0'
    input:
    path bedgraph_files
    path design_matrix
    path pmd_files

    output:
    tuple val("dmr"), path("DMRs.bed"), emit: dmr_bed optional true
    path "smoothed_counts.txt", emit: smoothed optional true
    path "dss_status.txt", emit: status
    script:
    """
    Rscript - <<'EOF'
    library(DSS)
    library(BiocParallel)

    writeLines("started", "dss_status.txt")

    min_cov <- 5L
    min_samples_per_locus <- 2L
    max_loci <- 1000000L

    # ----------------------------
    # 2. Cargar datos desde bedGraph
    # ----------------------------
    files <- Sys.glob("*.bedGraph")
    sample_count <- length(files)
    pmd_paths <- Sys.glob("*_PMDs.bed")

    build_pmd_mask <- function(paths) {
        if (length(paths) == 0) {
            return(list())
        }

        all_pmd <- lapply(paths, function(p) {
            if (file.info(p)\$size == 0) {
                return(NULL)
            }
            x <- tryCatch(
                read.table(p, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = ""),
                error = function(e) NULL
            )
            if (is.null(x) || nrow(x) == 0 || ncol(x) < 3) {
                return(NULL)
            }
            data.frame(
                chr = as.character(x[[1]]),
                start = as.integer(x[[2]]) + 1L,
                end = as.integer(x[[3]]),
                stringsAsFactors = FALSE
            )
        })

        all_pmd <- do.call(rbind, Filter(Negate(is.null), all_pmd))
        if (is.null(all_pmd) || nrow(all_pmd) == 0) {
            return(list())
        }

        split_pmd <- split(all_pmd, all_pmd[["chr"]])
        lapply(split_pmd, function(df) {
            df <- df[order(df[["start"]], df[["end"]]), c("start", "end"), drop = FALSE]

            cur_start <- df[["start"]][1]
            cur_end <- df[["end"]][1]
            merged <- vector("list", nrow(df))
            m <- 1L

            if (nrow(df) > 1) {
                for (i in 2:nrow(df)) {
                    s <- df[["start"]][i]
                    e <- df[["end"]][i]
                    if (s <= cur_end + 1L) {
                        cur_end <- max(cur_end, e)
                    } else {
                        merged[[m]] <- c(cur_start, cur_end)
                        m <- m + 1L
                        cur_start <- s
                        cur_end <- e
                    }
                }
            }

            merged[[m]] <- c(cur_start, cur_end)
            merged <- do.call(rbind, merged[seq_len(m)])
            colnames(merged) <- c("start", "end")
            merged
        })
    }

    filter_out_pmd <- function(df, pmd_mask) {
        if (length(pmd_mask) == 0 || nrow(df) == 0) {
            return(df)
        }

        keep <- rep(TRUE, nrow(df))
        chr_values <- unique(as.character(df[["chr"]]))

        for (chr in chr_values) {
            ranges <- pmd_mask[[chr]]
            if (is.null(ranges) || nrow(ranges) == 0) {
                next
            }

            idx_chr <- which(df[["chr"]] == chr)
            if (length(idx_chr) == 0) {
                next
            }

            pos <- df[["pos"]][idx_chr]
            j <- findInterval(pos, ranges[, "start"])
            
            # Initialize in_pmd as FALSE, set to TRUE only for valid intervals
            in_pmd <- logical(length(pos))
            valid_j <- j > 0 & j <= nrow(ranges)
            if (any(valid_j)) {
                in_pmd[valid_j] <- pos[valid_j] <= ranges[pmax(j[valid_j], 1L), "end"]
            }
            
            keep[idx_chr] <- !in_pmd
        }

        df[keep, , drop = FALSE]
    }

    pmd_mask <- build_pmd_mask(pmd_paths)
    if (length(pmd_mask) > 0) {
        message(paste0(">> PMD mask cargada desde ", length(pmd_paths), " fichero(s) de MethylSeekR"))
    } else {
        message(">> No se detectaron PMDs; DSS se ejecuta sin filtrado PMD")
    }

    parse_bedgraph <- function(file) {
        first_line <- readLines(file, n = 1, warn = FALSE)
        skip_header <- length(first_line) > 0 && startsWith(first_line, "track")

        df <- read.table(
            file,
            header = FALSE,
            sep = "",
            stringsAsFactors = FALSE,
            quote = "",
            comment.char = "",
            skip = if (skip_header) 1 else 0,
            fill = TRUE
        )

        if (ncol(df) < 6) {
            stop(paste0("bedGraph malformado en ", file, ": se esperaban >= 6 columnas y llegaron ", ncol(df)))
        }

        # Keep only positions with some minimum coverage to reduce memory usage.
        cov <- as.integer(df[[5]]) + as.integer(df[[6]])
        keep <- !is.na(cov) & cov >= min_cov
        df <- df[keep, , drop = FALSE]
        cov <- cov[keep]

        out <- data.frame(
            chr = df[[1]],
            pos = as.integer(df[[2]]) + 1L,
            N = cov,
            X = as.integer(df[[5]])
        )

        out <- filter_out_pmd(out, pmd_mask)

        rm(df, cov, keep)
        gc(verbose = FALSE)
        out
    }

    write_light_smoothing <- function(df, outfile, window = 25L) {
        if (nrow(df) == 0) {
            write.table(data.frame(), file = outfile, sep = "\\t", quote = FALSE, row.names = FALSE)
            return(invisible(NULL))
        }

        ord <- order(df[["chr"]], df[["pos"]])
        df <- df[ord, , drop = FALSE]

        meth_raw <- ifelse(df[["N"]] > 0, df[["X"]] / df[["N"]], NA_real_)
        meth_smooth <- as.numeric(stats::filter(meth_raw, rep(1 / window, window), sides = 2))

        out <- data.frame(
            chr = df[["chr"]],
            pos = df[["pos"]],
            N = df[["N"]],
            X = df[["X"]],
            meth_raw = meth_raw,
            meth_smooth = meth_smooth
        )

        write.table(out, file = outfile, sep = "\\t", quote = FALSE, row.names = FALSE)
    }

    write_smooth_or_raw <- function(bsobj, outfile) {
        n_loci <- nrow(bsobj)

        if (is.null(n_loci) || n_loci < 2) {
            message(">> Muy pocos loci tras filtros; se exporta metilación sin smoothing")
            write.table(getMeth(bsobj, type = "raw"), file = outfile, sep = "\\t", quote = FALSE)
            return(invisible(NULL))
        }

        # Avoid unstable multicore behavior in some bsseq/DSS container builds.
        BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

        bsobj_smooth <- tryCatch(
            BSmooth(
                bsobj,
                ns = 25,
                h = 500,
                maxGap = 100000000,
                verbose = FALSE
            ),
            error = function(e) {
                message(paste0(">> BSmooth falló (", conditionMessage(e), "); se exporta metilación sin smoothing"))
                NULL
            }
        )

        if (is.null(bsobj_smooth)) {
            write.table(getMeth(bsobj, type = "raw"), file = outfile, sep = "\\t", quote = FALSE)
        } else {
            write.table(getMeth(bsobj_smooth, type = "smooth"), file = outfile, sep = "\\t", quote = FALSE)
        }
    }

    write_fallback_regions <- function(dmltest, outfile, top_n = 20L, flank = 50L) {
        p_values <- NULL
        if (!is.null(dmltest[["pval"]])) {
            p_values <- dmltest[["pval"]]
        } else if (!is.null(dmltest[["pvals"]])) {
            p_values <- dmltest[["pvals"]]
        }

        if (is.null(p_values)) {
            write.table(data.frame(), file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
            return(invisible(NULL))
        }

        p_values <- suppressWarnings(as.numeric(p_values))
        valid_idx <- which(is.finite(p_values))
        if (length(valid_idx) == 0) {
            write.table(data.frame(), file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
            return(invisible(NULL))
        }

        top_idx <- valid_idx[order(p_values[valid_idx])]
        top_idx <- head(top_idx, top_n)

        fallback <- data.frame(
            chr = as.character(dmltest[["chr"]][top_idx]),
            start = pmax(1L, as.integer(dmltest[["pos"]][top_idx]) - flank),
            end = as.integer(dmltest[["pos"]][top_idx]) + flank
        )

        fallback <- fallback[!is.na(fallback[["chr"]]) & !is.na(fallback[["start"]]) & !is.na(fallback[["end"]]), , drop = FALSE]
        write.table(fallback, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        invisible(NULL)
    }

    # ----------------------------
    # 3. Ejecutar smoothing o análisis diferencial DSS
    # ----------------------------
    has_design_file <- file.exists("${design_matrix}") && file.info("${design_matrix}")\$size > 0

    has_valid_groups <- FALSE
    design <- NULL
    group_levels <- character(0)
    if (has_design_file) {
        design <- read.table("${design_matrix}", header = TRUE, sep = "\\t", stringsAsFactors = FALSE)
        has_group_col <- "group" %in% colnames(design)
        if (has_group_col) {
            group_levels <- unique(trimws(design[["group"]][!is.na(design[["group"]])]))
            group_levels <- group_levels[group_levels != ""]
            message(paste0(">> Grupos detectados en diseño: ", paste(group_levels, collapse = ", ")))
            has_valid_groups <- length(group_levels) > 1
        }
    }

    if (sample_count <= 1) {
        message(">> Una sola muestra: se ejecuta smoothing en DSS")

        single_df <- parse_bedgraph(files[[1]])
        write_light_smoothing(single_df, "smoothed_counts.txt")
        writeLines("completed", "dss_status.txt")
        q(status = 0)
    }

    BSobj <- makeBSseqData(
        lapply(files, parse_bedgraph),
        sampleNames = sub(".bedGraph", "", basename(files), fixed = TRUE)
    )

    cov_mat <- getCoverage(BSobj, type = "Cov")
    keep <- rowSums(cov_mat >= min_cov, na.rm = TRUE) >= min_samples_per_locus
    BSobj <- BSobj[keep, ]

    if (nrow(BSobj) > max_loci) {
        message(paste0(">> Se limita a ", max_loci, " loci con mayor cobertura media para reducir RAM"))
        idx <- order(rowMeans(cov_mat[keep, , drop = FALSE], na.rm = TRUE), decreasing = TRUE)[seq_len(max_loci)]
        BSobj <- BSobj[idx, ]
    }

    rm(cov_mat, keep)
    gc(verbose = FALSE)

    if (has_valid_groups) {
        group_counts <- table(factor(trimws(design[["group"]]), levels = group_levels))
        single_replicate <- any(group_counts < 2)

        if (length(group_levels) == 2 && single_replicate) {
            message(">> Modo: dos grupos sin réplicas biológicas; DSS usará DMLtest con smoothing")

            sample_names <- colnames(BSobj)
            if (length(sample_names) < 2) {
                stop("No hay suficientes muestras para comparar dos grupos en DSS")
            }

            DMLtest <- DMLtest(
                BSobj,
                group1 = sample_names[1],
                group2 = sample_names[2],
                smoothing = TRUE
            )
            DMRs <- callDMR(DMLtest)
        } else {
            message(">> Modo: análisis diferencial (DML/DMR)")

            ncores_dss <- if (.Platform\$OS.type == "windows") {
                1L
            } else {
                max(1L, as.integer("${task.cpus}"))
            }

            if (single_replicate) {
                message(">> Detectado grupo sin replicados biológicos: DSS usará estrategia single-replicate")
            }

            DMLfit <- DMLfit.multiFactor(BSobj, design = design, formula = ~group)

            dmltest_args <- list(
                DMLfit,
                coef = "group"
            )

            mf_formals <- names(formals(DMLtest.multiFactor))
            if ("ncores" %in% mf_formals) {
                dmltest_args[["ncores"]] <- ncores_dss
            }
            if (single_replicate && "smoothing" %in% mf_formals) {
                dmltest_args[["smoothing"]] <- TRUE
            }
            if (single_replicate && "equal.disp" %in% mf_formals) {
                dmltest_args[["equal.disp"]] <- TRUE
            }

            DMLtest <- do.call(DMLtest.multiFactor, dmltest_args)
            DMRs <- callDMR(DMLtest)
        }

        if (is.null(DMRs) || nrow(DMRs) == 0) {
            message(">> No DMR found; se generan regiones candidatas de respaldo desde los mejores CpGs")
            write.table(data.frame(), file = "DMRs_annotated.txt", sep = "\t", quote = FALSE, row.names = FALSE)
            write_fallback_regions(DMLtest, "DMRs.bed", top_n = 20L)
        } else {
            write.table(DMRs, file = "DMRs_annotated.txt", sep = "\t", quote = FALSE, row.names = FALSE)

            bed <- data.frame(
                chr = DMRs[["chr"]],
                start = DMRs[["start"]],
                end = DMRs[["end"]]
            )
            write.table(bed, file = "DMRs.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        }

        p_values <- NULL
        if (!is.null(DMLtest[["pval"]])) {
            p_values <- DMLtest[["pval"]]
        } else if (!is.null(DMLtest[["pvals"]])) {
            p_values <- DMLtest[["pvals"]]
        }

        p_values <- as.numeric(p_values)
        p_values <- p_values[is.finite(p_values)]
        if (length(p_values) > 0) {
            pdf("DMR_statistics.pdf")
            hist(p_values, main = "P-value distribution")
            dev.off()
        } else {
            message(">> No hay p-values numéricos para graficar histograma")
        }
    } else {
        message(">> Más de una muestra, pero sin diseño válido con al menos 2 grupos; se ejecuta smoothing en DSS")
        write_smooth_or_raw(BSobj, "smoothed_counts.txt")
    }

    writeLines("completed", "dss_status.txt")
    EOF
    """
}