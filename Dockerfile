FROM rocker/rstudio:4.3.3

# Set working directory
WORKDIR /home/rstudio/app

# System dependencies
RUN apt-get update && apt-get install -y curl tar sed \
    libgmp-dev \ 
    libgraphviz-dev \
    graphviz \
    libcurl4-openssl-dev \
    libblas-dev \
    liblapack-dev \
    libssl-dev \
    gfortran \
    libhdf5-dev \
    libpng-dev \
    libfreetype6-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libnlopt-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre2-dev \ 
    wget \
    bzip2 \
    ca-certificates \
    libpython3-dev \
    pandoc \
    pandoc-citeproc \
    libmariadb-dev \
    libpq-dev \
    libmagick++-dev \
    libcairo2-dev \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Miniconda
RUN ARCH=$(uname -m) && \
    case "$ARCH" in \
        x86_64)  MINICONDA=Miniconda3-latest-Linux-x86_64.sh ;; \
        aarch64) MINICONDA=Miniconda3-latest-Linux-aarch64.sh ;; \
        *) echo "Unsupported architecture: $ARCH" && exit 1 ;; \
    esac && \
    wget -O /tmp/miniconda.sh "https://repo.anaconda.com/miniconda/$MINICONDA" && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh
    
# Add conda to PATH
ENV PATH="/opt/conda/bin:$PATH"

# Configure conda and create Python environment with packages
RUN conda create -n py-rstudio -c conda-forge --override-channels \
    python=3.10 numpy pandas pycox lifelines scikit-survival tqdm numba -y

# Set default R library path
ENV R_LIBS=/home/rstudio/app/r-lib
RUN mkdir -p /home/rstudio/app/r-lib && \
    echo "R_LIBS=/home/rstudio/app/r-lib" >> /home/rstudio/app/.Renviron

RUN install2.r --error remotes \ 
    BiocManager \
    RcppArmadillo \ 
    rstan \
    StanHeaders && \
    Rscript -e "remotes::install_version('Rcpp', version = '1.0.13', repos = 'https://cloud.r-project.org')" && \
    Rscript -e "stopifnot(packageVersion('Rcpp') >= '1.0.13')"
RUN Rscript -e "print(.libPaths())"

# Install the pysurvival from the tar
COPY pysurvivalR/ pysurvivalR/
RUN R CMD build pysurvivalR && \
    R CMD INSTALL --library=/usr/local/lib/R/site-library pysurvivalR_*.tar.gz

COPY install_packages.R ./
RUN Rscript install_packages.R

# Point reticulate to the right Python binary
RUN echo "RETICULATE_PYTHON=/opt/conda/envs/py-rstudio/bin/python" >> /home/rstudio/app/.Renviron

# Copy rest of the project files
#COPY . .

# In terminal to build image: 
#  docker build -t cindex_multiverse_project .
# In terminal to run the image in the container:
#  docker run -d -e PASSWORD=mypass123 -p 8787:8787 --name new_c-index cindex_multiverse_project
# (run in Rstudio username: rstudio, password: mypass123)

