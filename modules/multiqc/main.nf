process multiQC {
    label 'single'
    container 'multiqc/multiqc:v1.33'
    input:
    path inputs

    output:
    path "multiqc_report.html", emit: report

    script:
    """
    multiqc -f --outdir . $inputs
    """
}
