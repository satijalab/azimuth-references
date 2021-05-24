# Dockerfile for building references
FROM satijalab/seurat:4.0.1

# Install other R dependencies
RUN R --no-echo --no-restore --no-save -e "install.packages('feather')"

# Install Azimuth (0.3.2)
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('immunogenomics/presto')"
RUN R --no-echo --no-restore --no-save -e "install.packages(c('DT', 'googlesheets4', 'shinyBS', 'shinydashboard', 'shinyjs'))"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install('glmGamPoi')"
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('satijalab/azimuth', ref = 'v0.3.0', dependencies = FALSE)"

# Install python dependencies
RUN pip3 install scanpy
RUN pip3 install zarr

