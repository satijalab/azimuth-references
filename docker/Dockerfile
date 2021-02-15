# Dockerfile for building references
FROM satijalab/seurat:4.0.0

# Install other R dependencies
RUN R --no-echo --no-restore --no-save -e "install.packages('feather')"

# Install Azimuth (0.2.0.9002)
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('immunogenomics/presto')"
RUN R --no-echo --no-restore --no-save -e "install.packages(c('DT', 'googlesheets4', 'shinyBS', 'shinydashboard', 'shinyjs'))"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install('glmGamPoi')"
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('satijalab/azimuth', ref = '4eda381f022edbd03a709bd12f3fe6951321cc39', dependencies = FALSE)"

CMD [ "R" ]
