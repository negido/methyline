process trimgalore {
    tag "$sampleId"
    label 'medium'
    container 'evolbioinfo/trimgalore:v0.6.11'
    containerOptions '--entrypoint ""'

    input:
    tuple val(sampleId), path(reads)

    output:
    tuple val(sampleId), path("*_val_*.fq.gz"), emit: reads
    path("*_report.txt"), emit: report
    path("*_val_*_fastqc.html"), emit: fastqc_html
    path("*_val_*_fastqc.zip"), emit: fastqc_zip
    script:
    def paired = reads.size() == 2 ? '--paired' : ''
    """
    trim_galore --cores ${task.cpus} --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 5 --three_prime_clip_R2 5 --fastqc --gzip $paired $reads
    """
}
