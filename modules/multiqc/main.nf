process multiQC {
    label 'single'
    container 'multiqc/multiqc:v1.33'
    publishDir "${params.analysisName}/html", mode: 'copy'
    input:
    path inputs

    output:
    path "multiqc_report.html", emit: report

    script:
    """
    multiqc -f --outdir . $inputs
    """
}
