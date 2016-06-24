FROM robinlin/analysis-base
MAINTAINER robin robin890650@gmail.com

# install python dep
RUN apt-get update && \
    apt-get dist-upgrade -y --force-yes && \
    apt-get install git python-numpy python-scipy python-pip libatlas-base-dev gfortran python-dev -y --force-yes && \
    pip install cvxpy

# install required R packages
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com'))" > /root/.Rprofile && \
    R -e 'source("http://bioconductor.org/biocLite.R") ; biocLite(c("graph","RBGL","Rgraphviz"), ask=FALSE) ' && \
    R -e 'install.packages(c("R.oo", "methods", "R.matlab", "data.table", "distr", "entropy", "plyr", "FNN", "e1071", "stringr", "gRim", "gRain", "gRbase", "caret", "rPython", "crayon"))'

RUN git clone https://github.com/JiaMingLin/DPTable.git

# install MySQL
ADD my.cnf /my.cnf
RUN apt-get install mariadb-server python-mysqldb -y --force-yes && \
    cp /my.cnf /etc/mysql/my.cnf
    
# install Django
RUN pip install Django


WORKDIR /DPTable
VOLUME ["/DPTable/data", "/DPTable/output", "/workspace"]

EXPOSE 8080
