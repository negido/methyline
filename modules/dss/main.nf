process dss {
    tag 'dss'
    label 'high'
    container 'quay.io/biocontainers/bioconductor-dss:2.58.0--r45h01b2380_0'
    input:
    path bedgraph_files
    path pmd_files
    path design_matrix

    output:
    tuple val("dmr"), path("DMRs.bed"),     emit: dmr_bed   optional true
    path "smoothed_counts.txt",             emit: smoothed  optional true
    path "DMR_statistics.pdf",              emit: dmr_stats optional true
    path "dss_status.txt",                  emit: status


    script:
    """
    export DESIGN="${design_matrix}"

    Rscript - <<'EOF'
    library(DSS)

    writeLines("started", "dss_status.txt")

    # Leer bedGraphs
    files <- Sys.glob("*.bedGraph")
    if (length(files) < 2) stop("Se necesitan al menos 2 muestras.")

    read_bg <- function(f) {
        df  <- read.table(f, header = FALSE, sep = "", skip = 1, fill = TRUE)
        cov <- as.integer(df[[5]]) + as.integer(df[[6]])
        df  <- df[!is.na(cov) & cov >= 5, , drop = FALSE]
        cov <- cov[!is.na(cov) & cov >= 5]
        data.frame(chr = df[[1]], pos = as.integer(df[[2]]) + 1L,
                   N   = cov,     X   = as.integer(df[[5]]))
    }

    BSobj <- makeBSseqData(
        lapply(files, read_bg),
        sampleNames = sub(".bedGraph", "", basename(files), fixed = TRUE)
    )

    # Filtrar loci con cobertura suficiente
    keep  <- rowSums(getCoverage(BSobj, type = "Cov") >= 5) >= 2
    BSobj <- BSobj[keep, ]

    # Diseño y grupos
    design <- read.table(Sys.getenv("DESIGN"), header = TRUE, sep = "", stringsAsFactors = FALSE)
    groups <- unique(trimws(design[["group"]]))
    message(paste0(">> Grupos: ", paste(groups, collapse = ", ")))

    # Test diferencial
    if (length(groups) == 2) {
        sn      <- colnames(BSobj)
        DMLtest <- DMLtest(BSobj, group1 = sn[design\$group == groups[1]],
                                  group2 = sn[design\$group == groups[2]],
                                  smoothing = TRUE)
    } else {
        DMLfit  <- DMLfit.multiFactor(BSobj, design = design, formula = ~group)
        DMLtest <- DMLtest.multiFactor(DMLfit, coef = "group")
    }

    # Llamar DMRs y guardar
    DMRs <- callDMR(DMLtest)
    if (!is.null(DMRs) && nrow(DMRs) > 0) {
        write.table(DMRs[, c("chr","start","end")], "DMRs.bed",
                    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        pdf("DMR_statistics.pdf")
        p <- na.omit(as.numeric(DMLtest[["pval"]] %||% DMLtest[["pvals"]]))
        hist(p, main = "p-values DSS", xlab = "p-value", col = "steelblue", border = "white")
        dev.off()
    } else {
        message(">> No se detectaron DMRs.")
    }

    writeLines("completed", "dss_status.txt")
    EOF
    """
}