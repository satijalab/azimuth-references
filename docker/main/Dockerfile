# Dockerfile for building references
FROM satijalab/azimuth:0.4.2

# Install other R dependencies
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
RUN R --no-echo --no-restore --no-save -e "install.packages(c('feather', 'pagoda2'))"

# Install synapse for download
RUN pip3 install synapseclient

# Install python dependencies
RUN pip3 install scanpy

CMD [ "R" ]
