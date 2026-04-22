process annotatr {
	tag "$sampleId"
	label 'medium'
    container 'quay.io/biocontainers/bioconductor-annotatr:1.36.0--r45hdfd78af_0'
	input:
	tuple val(sampleId), path(dmr_bed)
    val referenceGenome

	output:
	tuple val(sampleId), path("*_annotated.txt"), emit: annotations

	script:
    """
    Rscript -e '
    library(annotatr)
    library(GenomicRanges)
    library(readr)
        library(dplyr)

        # ----------------------------
        # 1. Leer BED de regiones
        # ----------------------------
        dmr <- tryCatch(
            read_tsv("${dmr_bed}", col_names = c("chr", "start", "end"), show_col_types = FALSE, progress = FALSE),
            error = function(e) tibble::tibble(chr = character(), start = integer(), end = integer())
        )

        dmr <- dmr %>%
            transmute(
                chr = as.character(chr),
                start = suppressWarnings(as.integer(start)),
                end = suppressWarnings(as.integer(end))
            ) %>%
            filter(!is.na(chr), chr != "", !is.na(start), !is.na(end), end >= start)

        if (nrow(dmr) == 0) {
            write.table(
                data.frame(),
                file = "${sampleId}_annotated.txt",
                sep = "\\t",
                quote = FALSE,
                row.names = FALSE,
                col.names = TRUE
            )
            quit(save = "no", status = 0)
        }

        gr <- GRanges(
            seqnames = dmr[["chr"]],
            ranges = IRanges(start = dmr[["start"]], end = dmr[["end"]])
        )

    # ----------------------------
    # 2. Seleccionar anotaciones
    # ----------------------------
    genome <- "${referenceGenome}"

    core_annotations <- c(
        paste0(genome, "_cpg_islands"),
        paste0(genome, "_cpg_shores"),
        paste0(genome, "_cpg_shelves")
    )

    requested_annotations <- c(paste0(genome, "_basicgenes"), core_annotations)

    annots <- tryCatch(
        build_annotations(genome = genome, annotations = requested_annotations),
        error = function(e) {
            message(paste0(">> annotatr: basicgenes no disponible (", conditionMessage(e), "). Se continúa con anotaciones CpG."))
            build_annotations(genome = genome, annotations = core_annotations)
        }
    )

    # ----------------------------
    # 3. Anotar regiones
    # ----------------------------
    annotated <- tryCatch(
        annotate_regions(
            regions = gr,
            annotations = annots,
            ignore.strand = TRUE,
            quiet = FALSE
        ),
        error = function(e) {
            if (grepl("No annotations intersect the regions", conditionMessage(e), fixed = TRUE)) {
                message(">> annotatr: no hay intersecciones con las anotaciones; se exporta tabla vacía")
                NULL
            } else {
                stop(e)
            }
        }
    )

    # ----------------------------
    # 4. Convertir a tabla y exportar
    # ----------------------------
    annotated_df <- if (is.null(annotated)) data.frame() else as.data.frame(annotated)

    write.table(
        annotated_df,
        file = "${sampleId}_annotated.txt",
        sep = "\\t",
        quote = FALSE,
        row.names = FALSE
    )
    '
    """
}