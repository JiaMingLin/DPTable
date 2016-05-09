FROM robinlin/analysis-base
MAINTAINER robin robin890650@gmail.com

# install python dep
RUN apt-get update && \
    apt-get dist-upgrade -y --force-yes && \
    apt-get install git python-numpy python-scipy python-pip libatlas-base-dev gfortran python-dev -y --force-yes && \
    pip install cvxpy

RUN git clone https://github.com/JiaMingLin/DPTable.git

# install required R packages
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com'))" >.Rprofile && \
    R -e 'source("http://bioconductor.org/biocLite.R") ; biocLite(c("graph","RBGL","Rgraphviz"), ask=FALSE) ' && \
    R -e 'install.packages(c("R.oo", "methods", "R.matlab", "data.table", "distr", "entropy", "plyr", "FNN", "e1071", "stringr", "gRim=0.1.17", "gRain=1.2.4", "gRbase=1.7.2"))'

