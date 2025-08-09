## The C-index Multiverse

The C-index is a discrimination metric used in time-to-event prediction modeling. However, there are multiple implementations and choices that condition performance results. In this work, we study the effects of ties, censoring, time truncation and input transformations in the model ranking. 

Given the multiple libraries used in this work, from R to python, we provide a Docker image with all code dependencies pre-installed.


### Download files: 

```bash

git clone https://github.com/BBolosSierra/CindexMultiverse.git
cd CindexMultiverse

```

### Docker:

Install docker from https://www.docker.com

For Mac computers follow https://docs.docker.com/desktop/setup/install/mac-install/


```bash

docker pull ghcr.io/bbolossierra/cindex_multiverse_project:1.0.4

```

### Interactively run Rstudio server:

```bash
docker container run --mount type=bind,source="$(pwd)",target=/home/rstudio/project -e PASSWORD=yourpassword -p 8787:8787 ghcr.io/bbolossierra/cindex_multiverse_project:1.0.4

```

Open your browser: http://localhost:8787

To log in; Username: rstudio
Password: yourpassword (or whatever you set)

Use the Files pane to navigate to /CindexMultiverse

### To remove the container:

```bash

docker rmi ghcr.io/bbolossierra/cindex_multiverse_project:1.0.4

```

### R libraries


| Package | Version |
|---------|---------|
| rhdf5 | 2.46.1 |
| graph | 1.80.0 |
| Rgraphviz | 2.46.0 |
| parallelly | 1.44.0 |
| survival | 3.8.3 |
| data.table | 1.17.2 |
| stringr | 1.5.1 |
| dplyr | 1.1.4 |
| tidyr | 1.3.1 |
| gridExtra | 2.3 |
| ggalluvial | 0.12.5 |
| cowplot | 1.1.3 |
| patchwork | 1.3.0 |
| arrow | 20.0.0 |
| reticulate | 1.42.0 |
| rmarkdown | 2.29 |
| knitr | 1.50 |
| Hmisc | 5.2.3 |
| rms | 6.8.1 |
| prodlim | 2025.4.28 |
| pec | 2023.4.12 |
| riskRegression | 2025.5.20 |
| randomForestSRC | 3.4.0 |
| caret | 7.0.1 |
| doFuture | 1.0.2 |
| future | 1.40.0 |
| progressr | 0.15.1 |
| foreach | 1.5.2 |
| flexsurv | 2.3.2 |
| furrr | 0.3.1 |
| gtsummary | 2.2.0 |
| kableExtra | 1.4.0 |
| survivalmodels | 0.1.191 |
| tidyverse | 2.0.0 |
| cBioPortalData | 2.14.2 |


### Python packages

| Package           | Version  |
|-------------------|----------|
| numpy             | 2.0.2    |
| pandas            | 2.2.3    |
| pycox             | 0.3.0    |
| lifelines         | 0.30.0   |
| scikit-survival   | 0.23.1   |
| tqdm              | 4.67.1   |
| numba             | 0.60.0   |



