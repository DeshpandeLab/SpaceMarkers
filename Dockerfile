FROM rocker/tidyverse:4.5

COPY . /spacemarkers
WORKDIR /spacemarkers

RUN sudo apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install libhdf5-dev build-essential libglpk-dev libuv1-dev patch -y

RUN Rscript -e 'install.packages("pak");\
                pak::pkg_install("CoGAPS");\
                pak::pkg_install("ComplexHeatmap");\
                pak::pkg_install("BiocNeighbors");\
                pak::pkg_install("jinworks/CellChat");\
                pak::pkg_install("scverse/anndataR");'

RUN Rscript -e 'pak::local_install(dependencies = TRUE);'
