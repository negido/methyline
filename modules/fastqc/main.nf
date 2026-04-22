process fastqc {
	tag "$sampleId"
	label 'medium'
    container 'biocontainers/fastqc:v0.11.9_cv8'
	input:
	tuple val(sampleId), path(reads)

	output:
	path("*_fastqc.html"), emit: html
	path("*_fastqc.zip"), emit: zip

	script:
	"""
	fastqc --threads ${task.cpus} $reads
	"""
}
