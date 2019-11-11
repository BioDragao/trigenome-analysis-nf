FROM nfcore/base:1.7
LABEL authors="Emilyn Costa and Abhinav Sharma" \
      description="Docker image containing all requirements for nf-core/trigenome pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-trigenome-1.0dev/bin:$PATH
