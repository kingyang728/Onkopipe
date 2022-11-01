FROM continuumio/miniconda3
WORKDIR /work
VOLUME /work

VOLUME /data


# RUN mkdir ref
# RUN mkdir InputFastqDir

COPY Snakefile config.yaml environment.yaml ./
COPY envs/ ./envs/
COPY scripts/ ./scripts/
#RUN chown -R admin:admin ./
#RUN chmod 775 /
RUN conda env create -f environment.yaml 
## Make RUN commands use the new environment:
# SHELL ["conda", "run", "-n", "pipeline_env", "/bin/bash", "-c"]
ENV PATH /opt/conda/envs/pipeline_env/bin:$PATH
RUN /bin/bash -c "source activate pipeline_env"


# RUN conda config --add channels conda-forge \
#     && conda config --add channels bioconda \
#     && conda install -c conda-forge mamba \
#     && mamba create -c conda-forge -c bioconda -n snakemake snakemake



# RUN snakemake -v
# RUN snakemake -j all --use-conda -n
# RUN snakemake -j all --use-conda

# ENTRYPOINT	["snakemake","-j","all","--use-conda"]
# ENTRYPOINT	["snakemake","-v"]
# ENTRYPOINT ["conda", "run", "-n", "pipeline_env", "snakemake","-j","all","--use-conda"]
# ENTRYPOINT ["conda", "run", "-n", "pipeline_env", "snakemake","-v"]
# ENTRYPOINT	conda activate snakemake
