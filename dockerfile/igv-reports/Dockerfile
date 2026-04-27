FROM continuumio/miniconda3

WORKDIR /usr/src/app

# Create the environment:
RUN conda create -n igvreports python=3.7.1
RUN pip install igv-reports

# Make RUN commands use the new environment:
ENV PATH=/root/miniconda/bin:$PATH

RUN conda config --add channels r
RUN conda config --add channels bioconda
RUN conda install -y -c conda-forge  -c bioconda pysam