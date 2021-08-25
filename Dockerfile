FROM r-base:4.1.1
MAINTAINER Anthony S. Castanza <acastanza@ucsd.edu>

USER root

RUN apt-get update --yes && apt-get install --no-install-recommends --yes \
  build-essential \
  cmake \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev

RUN R -e 'chooseCRANmirror(ind=1); install.packages(c("getopt","optparse","ROCR"))'

RUN R -e "sessionInfo()"
RUN rm -rf /tmp/downloaded_packages/

CMD ["Rscript", "--version"]

# build using this:
# docker build -t acastanza/gsea:1.4 .
