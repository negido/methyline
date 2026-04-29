FROM bioconductor/bioconductor_docker:3.21

RUN Rscript -e "\
    BiocManager::install(c(\
        'rGREAT',\
        'TxDb.Hsapiens.UCSC.hg19.knownGene',\
        'TxDb.Hsapiens.UCSC.hg38.knownGene',\
        'org.Hs.eg.db',\
        'GO.db',\
        'rtracklayer',\
        'GenomicRanges'\
    ), ask = FALSE, update = FALSE)\
"