process methylseekr {
    tag "$sampleId"
    label 'medium'
    container 'quay.io/biocontainers/bioconductor-methylseekr:1.50.0--r45hdfd78af_0'

    input:
    tuple val(sampleId), path(methylation_bed)
    path genome_fasta

    output:
    tuple val("${sampleId}_UMR"), path("*_UMRs.bed"), emit: umrs
    tuple val("${sampleId}_LMR"), path("*_LMRs.bed"), emit: lmrs
    tuple val("${sampleId}_PMD"), path("*_PMDs.bed"), emit: pmds

    script:
    """
    Rscript - <<'EOF'
    library(MethylSeekR)
    library(Biostrings)
    library(rtracklayer)

    skip_header <- startsWith(readLines("${methylation_bed}", n = 1, warn = FALSE), "track")

    bed <- read.table(
        "${methylation_bed}",
        header = FALSE,
        sep = "",
        stringsAsFactors = FALSE,
        quote = "",
        comment.char = "",
        skip = if (skip_header) 1 else 0,
        fill = TRUE
    )

    seq_idx <- Biostrings::fasta.index("${params.referenceGenome}.fa")
    len_col <- if ("seqlengths" %in% colnames(seq_idx)) "seqlengths" else "seqlength"
    name_col <- if ("seqnames" %in% colnames(seq_idx)) "seqnames" else "desc"
    seq_names <- vapply(
        strsplit(as.character(seq_idx[[name_col]]), " ", fixed = TRUE),
        `[`,
        character(1),
        1
    )
    seqLengths <- stats::setNames(as.integer(seq_idx[[len_col]]), seq_names)

    methylome_file <- tempfile(fileext = ".tab")
    write.table(
        data.frame(
            chr = bed[[1]],
            pos = as.integer(bed[[2]]) + 1L,
            T = as.integer(bed[[5]]) + as.integer(bed[[6]]),
            M = as.integer(bed[[5]])
        ),
        file = methylome_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )

    meth.gr <- readMethylome(FileName = methylome_file, seqLengths = seqLengths)
    chr.sel <- if ("chr1" %in% seqlevels(meth.gr)) "chr1" else seqlevels(meth.gr)[1]

    pmd_segments <- segmentPMDs(
        m = meth.gr,
        chr.sel = chr.sel,
        seqLengths = seqLengths,
        num.cores = 1,
        nCGbin = 101
    )

    if (is.list(pmd_segments) && !inherits(pmd_segments, "GRanges")) {
        pmd_segments <- Filter(Negate(is.null), pmd_segments)
        if (length(pmd_segments) > 0) {
            pmd_segments <- do.call(c, pmd_segments)
        }
    }

    if (!inherits(pmd_segments, "GRanges")) {
        stop(paste0("segmentPMDs returned unsupported object of class: ", paste(class(pmd_segments), collapse = ", ")))
    }

    pmds <- pmd_segments[pmd_segments\$type == "PMD"]
    export.bed(pmds, con = "${sampleId}_PMDs.bed")

    file.create("${sampleId}_UMRs.bed")
    file.create("${sampleId}_LMRs.bed")

    cat("MethylSeekR completed\\n")
    EOF
    """
}
