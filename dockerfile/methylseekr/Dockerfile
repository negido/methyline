FROM quay.io/biocontainers/bioconductor-methylseekr:1.50.0--r45hdfd78af_0

# Instalar BiocManager primero
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

# Instalar KernSmooth para evitar fallos de segmentacion en MethylSeekR
RUN R -e "install.packages('KernSmooth', repos='https://cloud.r-project.org')"

# Instalar los genomas
RUN R -e "BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38','BSgenome.Hsapiens.UCSC.hg19'), ask=FALSE, update=FALSE)"