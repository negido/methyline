process methylseekr {
    tag "$sampleId"
    label 'medium'
    container 'negido/methylseekr_hs_bsgenome:v1.1'

    publishDir "${params.analysisName}/bed", mode: 'copy', pattern: "*_PMDs.bed"
    publishDir "${params.analysisName}/bed", mode: 'copy', pattern: "*_UMRs.bed"
    publishDir "${params.analysisName}/bed", mode: 'copy', pattern: "*_LMRs.bed"
    publishDir "${params.analysisName}/pdf", mode: 'copy', pattern: "*_segmentPMDs.pdf"
    publishDir "${params.analysisName}/pdf", mode: 'copy', pattern: "*_plotPMDSegmentation.pdf"
    publishDir "${params.analysisName}/pdf", mode: 'copy', pattern: "*_segmentUMRsLMRs.pdf"
    publishDir "${params.analysisName}/pdf", mode: 'copy', pattern: "*_finalSegmentation.pdf"
    input:
    tuple val(sampleId), path(methylation_bed)
    output:
    tuple val(sampleId), path("${sampleId}_PMDs.bed"), emit: pmds
    tuple val(sampleId), path("${sampleId}_segmentPMDs.pdf"), emit: pmds_plot, optional: true
    tuple val(sampleId), path("${sampleId}_plotPMDSegmentation.pdf"), emit: pmds_example, optional: true
    tuple val(sampleId), path("${sampleId}_UMRs.bed"), emit: umrs
    tuple val(sampleId), path("${sampleId}_LMRs.bed"), emit: lmrs
    tuple val(sampleId), path("${sampleId}_segmentUMRsLMRs.pdf"), emit: segmentation_plot, optional: true
    tuple val(sampleId), path("${sampleId}_finalSegmentation.pdf"), emit: final_plot, optional: true
    script:
    """
    Rscript - <<'EOF'
    library(MethylSeekR)
    library(GenomicRanges)
    library(rtracklayer)

    if ("${params.referenceGenome}" == "hg38") {
        library(BSgenome.Hsapiens.UCSC.hg38)
    } else {
        library(BSgenome.Hsapiens.UCSC.hg19)
    }

    seqLengths <- seqlengths(Hsapiens)

    first_line <- readLines("${methylation_bed}", n = 1, warn = FALSE)
    skip_n     <- if (startsWith(first_line, "track")) 1L else 0L

    bed <- read.table(
        "${methylation_bed}",
        header           = FALSE,
        sep              = "",
        stringsAsFactors = FALSE,
        quote            = "",
        comment.char     = "",
        skip             = skip_n,
        fill             = TRUE
    )

    bed <- bed[complete.cases(bed[, c(1, 2, 5, 6)]), ]

    methylome_file <- tempfile(fileext = ".tab")
    write.table(
        data.frame(
            chr = bed[[1]],
            pos = as.integer(bed[[2]]) + 1L,
            T   = as.integer(bed[[5]]) + as.integer(bed[[6]]),
            M   = as.integer(bed[[5]])
        ),
        file      = methylome_file,
        sep       = "\\t",
        quote     = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )

    meth.gr <- readMethylome(FileName = methylome_file, seqLengths = seqLengths)

    chr.sel <- if ("chr1" %in% seqlevels(meth.gr)) "chr1" else seqlevels(meth.gr)[1]

    pmd_segments <- segmentPMDs(
        m          = meth.gr,
        chr.sel    = chr.sel,
        pdfFilename = "${sampleId}_segmentPMDs.pdf",
        seqLengths = seqLengths,
        num.cores  = 1,
        nCGbin     = 101
    )

    savePMDSegments(
        PMDs            = pmd_segments,
        GRangesFilename = "${sampleId}_PMDs.rds",
        TableFilename   = "${sampleId}_PMDs.bed"
    )

    plotPMDSegmentation(
        m           = meth.gr,
        segs        = pmd_segments,
        numRegions  = 1,
        pdfFilename = "${sampleId}_plotPMDSegmentation.pdf",
        minCover    = 1
    )

    umr_lmr <- segmentUMRsLMRs(
        m              = meth.gr,
        meth.cutoff    = 0.5,
        nCpG.cutoff    = 3,
        PMDs           = pmd_segments,
        pdfFilename    = "${sampleId}_segmentUMRsLMRs.pdf",
        num.cores      = 1,
        myGenomeSeq    = Hsapiens,
        seqLengths     = seqLengths,
        nCpG.smoothing = 3,
        minCover       = 5
    )

    saveUMRLMRSegments(
        segs            = umr_lmr,
        GRangesFilename = "${sampleId}_UMRLMRs.bed"
    )

    rtracklayer::export(umr_lmr[umr_lmr\$type == "UMR"], "${sampleId}_UMRs.bed", format = "bed")
    rtracklayer::export(umr_lmr[umr_lmr\$type == "LMR"], "${sampleId}_LMRs.bed", format = "bed")

    plotFinalSegmentation(
        m              = meth.gr,
        segs           = umr_lmr,
        PMDs           = pmd_segments,
        meth.cutoff    = 0.5,
        numRegions     = 3,
        pdfFilename    = "${sampleId}_finalSegmentation.pdf",
        minCover       = 1,
        nCpG.smoothing = 3
    )

    cat("MethylSeekR completed\\n")
    EOF
    """
}