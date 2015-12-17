## DPTable
A differentially private tabular data publishing algorithm using sampling and junction tree inference.

####Installation:
Depends:
R library:  R.oo, methods, R.matlab, data.table, distr, entropy, plyr, FNN, e1071, stringr, gRim=0.1.17,  gRain=1.2.4,  gRbase=1.7.2

Matlab & [cvx matlab software](http://cvxr.com/cvx/)

MAKE SURE you have MATLAB installed, otherwise the script will keep waiting to connect matlab.

To install R graphical modeling packages(gRain, gRim, gRbase),
please see the instructions on the website: [http://people.math.aau.dk/~sorenh/software/gR/](http://people.math.aau.dk/~sorenh/software/gR/)
The packages listed above use the graph, RBGL and Rgraphviz packages. These packages are NOT on CRAN but on bioconductor. Please execute the following commands before installing these packages:
```sh
>source("http://bioconductor.org/biocLite.R")
>biocLite(c("graph","RBGL","Rgraphviz"))
```

####Instructions for running the program:
`run-category.R` under folder `exp` is an example script to run junction tree algorithm, generate synthetic data and process k-way marginal queries(be can random query or all k-way marginals) under differential privacy. 

Input files:
The dataset files should be under `data` folder. For each dataset, it needs one binary/categorical dataset file and one domain files that describes all attributes. Numberical values should be coarsened into category first. The dataset file should be named with format [dataset name]+"-coarse.dat". The domain file should be named with [dataset name]+"-coarse.domain".

Output files:
All output files would be put under `output` folder. Each folder is named with dataset name, parameter setting for CV and timestamp. All junction tree files will be saved. The statistics for errors of k-way marginal will also be saved. For synthetic dataset, it will be saved and named with parameter settings(ending with "-sim.dat")

To run the example script `run-category.R`:

In Rstudio, 
first, choose the working directory to the source code directory,

then run
```sh
> system("Rscript exp/run-category.R -f Data3 -e1 0.2 -e2 0.8")
```
```sh
> system("Rscript exp/run-category.R -f Data3 -e1 0.2 -e2 0.8 -nrun 1 -CV 0.3")
```

#####Arguments:
* -f: dataset name
* -e1: epsilon value for constructing noisy junction tree
* -e2: epsilon value for injecting noise into the marginal tables
* -nrun: num of runs, 10 as default
* -CV: numeric value in (0, 1),  threshold  Cramer’s V value for picking correlated attributes pairs, 0.2 as default,typically choose 0.2 for weakly correlated datasets; 0.3 for highly correlated datasets
* -q: boolean value indicates whether to process query on the releasing tables, True as default
* -sim: boolean value indicates whether to simulate and generate synthetic data, False as default. If setting as True, the program shall generate one copy of synthetic data with raw data's size in the output directory.

The parameter `k` for k-way marginal can be reset inside the script by resetting variables as following:
>random_kways = c(4, 6)
>all_kways = c(2, 3)
