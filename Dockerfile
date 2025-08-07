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
    
#ENV RETICULATE_PYTHON=/opt/conda/envs/py-rstudio/bin/python
# Add conda to PATH
#ENV PATH="/opt/conda/bin:$PATH"
# Add conda to PATH
ENV PATH="/opt/conda/bin:$PATH"

ENV LD_LIBRARY_PATH="/opt/conda/envs/py-rstudio/lib:$LD_LIBRARY_PATH"

#RUN conda create -n py-rstudio -c conda-forge --override-channels \
#    python=3.9.21 -y
    
# Configure conda and create Python environment with packages
ENV CONDA_ALWAYS_YES=true

# Create base environment
RUN conda create -n py-rstudio python=3.9.21 -c conda-forge --override-channels

# Install remaining packages via `conda run` so the env is properly used
RUN conda run -n py-rstudio conda install -c conda-forge --override-channels \
    numpy=2.0.2 pandas=2.2.3 tqdm=4.67.1 numba=0.60.0 lifelines=0.30.0 \
    scikit-survival=0.23.1 openssl=3.2.0 h5py libcurl

# Install pycox via pip (more stable)
RUN conda run -n py-rstudio pip install torch torchtuples pycox==0.3.0 && \
    mkdir -p /opt/conda/envs/py-rstudio/lib/python3.9/site-packages/pycox/datasets/data && \
    chmod -R a+rwX /opt/conda/envs/py-rstudio/lib/python3.9/site-packages/pycox/datasets


# Set default R library path
ENV R_LIBS=/home/rstudio/app/r-lib
RUN mkdir -p /home/rstudio/app/r-lib && \
    echo "R_LIBS=/home/rstudio/app/r-lib" >> /home/rstudio/app/.Renviron

# Ensure reticulate uses the correct Python
RUN echo "RETICULATE_PYTHON=/opt/conda/envs/py-rstudio/bin/python" >> /home/rstudio/app/.Renviron

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

RUN echo "LD_LIBRARY_PATH=/opt/conda/envs/py-rstudio/lib" >> /home/rstudio/app/.Renviron

# Copy rest of the project files
#COPY . .

# In terminal to build image: 
#  docker build -t cindex_multiverse_project .
# In terminal to run the image in the container:
#  docker run -d -e PASSWORD=mypass123 -p 8787:8787 --name new_c-index cindex_multiverse_project
# (run in Rstudio username: rstudio, password: mypass123)

# # Steps followed for github-docker
# # tag it
# docker tag cindex_multiverse_project ghcr.io/bbolossierra/cindex_multiverse_project:1.0.0
# # create token I did a fine grained, added just the repo I needed and selected: Actions Read and Write
# # authentification:
# echo github_token | docker login ghcr.io -u bbolossierra --password-stdin
# docker push ghcr.io/bbolossierra/cindex_multiverse_project:1.0.0
# docker pull ghcr.io/bbolossierra/cindex_multiverse_project:1.0.0


