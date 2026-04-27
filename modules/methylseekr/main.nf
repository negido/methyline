process methylseekr {
    tag "$sampleId"
    label 'medium'
    container 'negido/methylseekr_hs_bsgenome:v1.0'
    input:
    tuple val(sampleId), path(methylation_bed)
    output:
    tuple val(sampleId), path("${sampleId}_PMDs.bed"), emit: pmds
    tuple val(sampleId), path("${sampleId}_UMRs.bed"), emit: umrs
    tuple val(sampleId), path("${sampleId}_LMRs.bed"), emit: lmrs
    tuple val(sampleId), path("${sampleId}_segmentUMRsLMRs.pdf"), emit: segmentation_plot
    tuple val(sampleId), path("${sampleId}_finalSegmentation.pdf"), emit: final_plot
    script:
    """
    Rscript - <<'EOF'
    library(MethylSeekR)
    library(GenomicRanges)
    library(BSgenome)
    library(rtracklayer)

    # ---- Cargar BSgenome ----
    if ("${params.referenceGenome}" == "hg38") {
        library(BSgenome.Hsapiens.UCSC.hg38)
        myGenomeSeq <- BSgenome.Hsapiens.UCSC.hg38
    } else {
        library(BSgenome.Hsapiens.UCSC.hg19)
        myGenomeSeq <- BSgenome.Hsapiens.UCSC.hg19
    }

    # ---- Seq lengths ----
    seqLengths <- seqlengths(myGenomeSeq)

    # ---- Leer BED ----
    first_line <- readLines("${methylation_bed}", n = 1, warn = FALSE)
    skip_header <- length(first_line) > 0 && startsWith(first_line, "track")

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

    bed <- bed[!is.na(bed[[1]]) & !is.na(bed[[2]]) & !is.na(bed[[5]]) & !is.na(bed[[6]]), ]

    # ---- Crear methylome ----
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

    # ---- PMDs ----
    pmd_segments <- segmentPMDs(
        m = meth.gr,
        chr.sel = chr.sel,
        seqLengths = seqLengths,
        num.cores = 1,
        nCGbin = 101
    )

    pmds <- pmd_segments[pmd_segments$type == "PMD"]
    rtracklayer::export(pmds, "${sampleId}_PMDs.bed", format = "bed")

    # ---- UMR / LMR ----
    umr_lmr <- segmentUMRsLMRs(
        m = meth.gr,
        meth.cutoff = 0.5,
        nCpG.cutoff = 3,
        PMDs = pmd_segments,
        pdfFilename = "${sampleId}_segmentUMRsLMRs.pdf",
        num.cores = 1,
        myGenomeSeq = myGenomeSeq,
        seqLengths = seqLengths,
        nCpG.smoothing = 3,
        minCover = 5
    )

    umrs <- umr_lmr[umr_lmr$type == "UMR"]
    lmrs <- umr_lmr[umr_lmr$type == "LMR"]

    rtracklayer::export(umrs, "${sampleId}_UMRs.bed", format = "bed")
    rtracklayer::export(lmrs, "${sampleId}_LMRs.bed", format = "bed")

    # ---- Plot final ----
    plotFinalSegmentation(
        m = meth.gr,
        segs = umr_lmr,
        PMDs = pmd_segments,
        meth.cutoff = 0.5,
        numRegions = 3,
        pdfFilename = "${sampleId}_finalSegmentation.pdf",
        minCover = 5,
        nCpG.smoothing = 3
    )

    cat("MethylSeekR completed\\n")
    EOF
    """
}
