#generated from 0420c98733ce24d9ef58cd966484da85c19a4236
# --platform=linux/amd64 to avoid 'no match for platform in the manifest' on M1
FROM rocker/tidyverse:4

COPY . /spacemarkers
WORKDIR /spacemarkers

RUN sudo apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install libhdf5-dev build-essential patch -y

RUN Rscript -e 'devtools::install_deps()'

#https://github.com/r-lib/devtools/issues/2395
RUN Rscript -e 'devtools::install()'
