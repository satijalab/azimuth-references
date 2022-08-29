# Dockerfile for building references
FROM satijalab/azimuth:0.4.5

# Install other R dependencies
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
RUN R --no-echo --no-restore --no-save -e "install.packages(c('feather', 'pagoda2', 'arrow', 'SeuratDisk', 'yaml', 'tidyverse', 'stringr'))"
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')"

# Install for downloading google drive files
RUN pip install gdown 

# Install synapse for download
RUN pip3 install synapseclient

# Install python dependencies
RUN pip3 install scanpy

CMD [ "R" ]

