FROM r-base:4.1.1
MAINTAINER Anthony S. Castanza <acastanza@ucsd.edu>

USER root

RUN R -e 'chooseCRANmirror(ind=1); install.packages(c("getopt","optparse","ROCR"))'

RUN R -e "sessionInfo()"
RUN rm -rf /tmp/downloaded_packages/

CMD ["Rscript", "--version"]

# build using this:
# docker build -t acastanza/gsea:1.4 .
