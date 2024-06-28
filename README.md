# Code repository for `The $\ell$-test: leveraging sparsity in the Gaussian linear model for improved inference' by Sengupta and Janson (2024).
This repository contains source codes for implementing the $\ell$-test-based procedures and for reproducing the figures from the paper. These procedures are alternatives to the traditional $t$-based inference
procedures in the Gaussian linear model, having exactly the same guarantees showing significant power gains when the coefficient vector is not too dense. In fact, our alternative to the $t$-test, which we
call the $\ell$-test achieves power close to the one-sided $t$-test for sparse coefficient vectors without any knowledge about the true sign of the alternative coefficient. Check out our
paper at [http://arxiv.org/abs/2406.18390](http://arxiv.org/abs/2406.18390) for more details. 

## Usage
The file `l_testing.R` contains the source code for implementing the $\ell$-test-based single testing methods and the unadjusted confidence interval (that doesn't condition on LASSO selection),
while the file `adjusted_l_testing.R` contains the codes for implementing the adjusted $\ell$-test-based confidence interval. Here are some basic examples of how to call these
functions.

```R
source('l_testing.R')
source('adjusted_l_testing.R')
set.seed(1)

n = 100
p = 50
s =  5
A = 2.3

X = matrix(rnorm(n*p), nrow = n)
X = apply(X,2,g)
#Supplying an X with normalized columns is recommended

beta = rep(0,p)
rand_ind = sample(1:p, size = s, replace = FALSE)
j = rand_ind[1]	#index to test
beta[rand_ind] = (1-2*rbinom(s,1,0.5))*A
y = as.numeric(X%*%beta + rnorm(n)) #data

pval_l = l.test(y,X,j)	#l-test for H_j:\beta_j = 0

pval_l = l.test(y-2.3*X[,j],X,j) #for testing H_j(2.3):\beta_j = 2.3

pval_l = l.test(y,X,j, lambda_cv = 0.01) #l-test for H_j:\beta_j = 0 with a supplied lambda for cross-validation

pval_l_adjusted = l.test(y,X,j, adjusted = TRUE, lambda = 0.01) #adjusted l-test for H_j:\beta_j = 0 valid conditionally on LASSO selection using penalty 0.01, and the penalty for the test statistic chosen using cross-validation	

gamma_range = seq(from = beta[j]-10, to = beta[j]+10, length.out = 100) #the grid of \gamma values to test on

ci_l = l.ci(y,X,j, gamma_range = gamma_range, coverage = 0.95) #l-CI

ci_l_adjusted = l.ci_adjusted(y,X,j, gamma_range = gamma_range, coverage = 0.95, lambda = 0.01) #post-selection l-CI for \beta_j valid conditionally on LASSO with penalty 0.01 selecting the coefficient

```

## Reproducing results from the paper
The folder `reproducibility` contains the codes for reproducing the simulation results from the paper. Some of these were directly run on a single computer while for others we used a computing cluster with a Slurm manager.
We provide scripts for submitting multiple jobs to the cluster in `reproducibility/scripts`, however note that the user needs to update them with the file locations, the partition names, etc., before running them.

We first discuss generation of the output files for creating the plots in the paper.
1. _For Figures 2 and 8_: For generating the files for the figures in Figure 8 and the left and center panels of Figure 2, run the file `single_test.R` by un-commenting Blocks 1,2 or 3 of paremeter specifications in the file.
   For generating the figure on the right panel of Figure 2, un-comment Block 4, set the value of the variable `cluster` in line 150 to TRUE and run the file 500 times (in the computing cluster) followed by running `combine_results_ltest_sparsity.R` to produce the desired
   file. You may use `scripts/run_test.sh` to submit jobs to the cluster.
2. _For Figures 3, 9, 10, 11 and 12_: Run the file `reproducibility/robustness_simulations.R`.
   For Figures 4, 5, 13, 14, 15, 16 and 17: Submit multiple jobs to the cluster: `sbatch --array=1-times reproducibility/scripts/run_ci.sh reproducibility/ci_parameter_files/filename.RData`, where
   for the left panel of Figure 4 and Figures 13 and 14 `filename` is `par_unadj_ci_amplitude` and `times` is 50, for the center panel of Figure 4 and Figure 15, `filename` is `par_unadj_ci_sparsity` and `times` is 500,
   for the right panel of Figure 4 and Figure 16, `filename` is `par_unad_ci_correlation` and `times` is 50, for Figures 5 and 17, execute two runs with `filename` taking values `par_adj_ci_amplitude` and
   `par_adj_ci_amplitude_known_sigma` respectively and `times` equal 1000 in both the cases. After this, run `reproducibility/scripts/combine_ci_results.sh` to get the desired files.
3. _For Figure 6_: Run the file `reproducibility/ltest_variability.R` followed by `reproducibility/combine_ltest_variability.R`.
4. _For Figure 7_: Run the file `reproducibility/three_lambda_choices.R`

Finally, to generate the plots, execute `repoducibility/plot_functions.R` followed by running the appropriate code-block from `reproducibility/plot_commands.R` to get the desired figure.

5. _For the claims in Section 5.5 of the paper_: This section involves an analysis on the HIV drug resistance dataset. The data is available in the folder `data`. Submit 313 jobs to the cluster:  `sbatch --array=1-313 ~/reproducibility/scripts/run_HIV_ci.sh` followed by executing `reproducibility/combine_HIV_ci.R` to get a data frame of the widths of the confidence intervals for each of the possible columns. Submit 16 jobs: `sbatch --array=1-16 ~/reproducibility/scripts/run_HIV_pval.sh` followed by `reproducibility/combine_HIV_pval.R` to get two data frames summarizing the number of raw discoveries and the number of gene level discoveries, respectively. In writing these functions, we have borrowed codes from the `R` implementation of [Luo et al. (2022)](https://arxiv.org/pdf/2208.09542).

## Reference
```
@article{SS-LJ:2024,
  title={The $\ell$-test: leveraging sparsity in the Gaussian linear model for improved inference},
  author={Sengupta, Souhardya and Janson, Lucas},
  journal={arXiv preprint arXiv:2406.18390},
  year={2024}
}
```
