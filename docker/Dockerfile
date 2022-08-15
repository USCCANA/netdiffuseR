FROM rocker/drd

RUN apt-get update && apt-get install \
	-y --no-install-recommends \
	libssl-dev
COPY pkgs.R .
COPY DESCRIPTION .
RUN Rscript pkgs.R && \
	Rscript -e 'utils::install.packages(readLines("pkgs.txt"))'

CMD ["bash"]

