## Concordance index multiverse


In this markdowns we are testing how the increase in censoring impacts the c index under various implementations from python and R. 

Since it requires python to be run, the reticulate library is necesary. 

Run the following in bash and install libraries:

```bash
conda create -n py-rstudio python=3.10
conda activate py-rstudio
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install numpy
conda install pandas
conda install pycox
conda install lifelines
conda install scikit-survival
conda install tqdm
# pysurvival has a problem with dependency on a deprecated package 'sklearn'. PysurvivalR has been created from the original python code by wrapping the c++ code for R.
# The R wrapper can be installed:
install.packages("./pysurvivalR.tar.gz", repos = NULL, type = "source")

```

The following set up is necessary for the python and R to be run in Rstudio. 

Make sure the conda environment created is loaded. If your bash default reticulate environment is set to a specific path, it might be necessary to unset with "Sys.unsetenv("RETICULATE_PYTHON")"

Run the following inside of a setup in R markdown ```{r setup}```

```r 
library(reticulate)

#Sys.unsetenv("RETICULATE_PYTHON") 
Sys.setenv(OMP_NUM_THREADS = "1")       # Limits OpenMP to 1 thread
Sys.setenv(NUMBA_NUM_THREADS = "1")     # Limits Numba to 1 thread
Sys.setenv(MKL_NUM_THREADS = "1")       # Limits Intel MKL to 1 thread
Sys.setenv(KMP_WARNINGS = "0")          # Disables OpenMP warnings
Sys.setenv(OPENBLAS_NUM_THREADS = "1")  # Limits OpenBLAS to 1 thread

# To compile Rcpp Rpysurvivalwrapper:
#Sys.setenv(CC = "clang")
#Sys.setenv(CXX = "clang++")
#Sys.setenv(CPPFLAGS = "-I/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/Rcpp/include")
#Sys.setenv(LDFLAGS = "-L/Library/Frameworks/R.framework/Resources/lib")


use_condaenv("/opt/homebrew/Caskroom/miniforge/base/envs/py-rstudio", required=TRUE)

```

Make sure you have:

```r
library(arrow)
library(caret)
library(riskRegression)
library(prodlim)
library(pec)
library(survival)
library(rhdf5)  #BiocManager::install("rhdf5")
library(randomForestSRC)
library(survAUC)
library(Hmisc)
library(dplyr)
library(gridExtra)
library(survC1)
library(pysurvivalR) # with install.packages("./pysurvivalR.tar.gz")
library(survivalmodels)

```

If trying bootstraping with parallelization of R libraries

```r
## load libraries
library(progressr) ## use progressr for procession updates
library(doFuture)  ## attaches also foreach and future
library(progressr) ## use progressr for procession updates

```

