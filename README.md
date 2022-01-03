# GxE-with-bayesian-process
This is a project on improving the GxE evaluation method by applying bayesian estimation to the original model.

## Introduction
The gene-environment interaction (GxE) is a common effect that alters normal genetic function in human body. To identify the GxE effect in large sample size of genome-wide association studies (GWAS) data, an [article](https://www.nature.com/articles/s41467-020-15107-0) (Sulc, 2020) presented a maximum likelihood method to estimate the contribution of GxE to continuous traits taking into account all interacting environmental variables, without the need to measure any.

Here we reproduce the method and result of this article, and generate our own simulation data to evaluate the efficacy of the method. Then we improved the method by applying bayesian estimation, Gaussian process, and Metropolis-Hasting Partially Collapsed Gibbs. The final comparison of two methods are also included.

## Files
* [code repeat](`./code repeat`) directory includes all method and figure generating codes to reproduce the result of the original article.
* [analysis](./analysis) directory includes the data simulation and analysis process to evaluate the method effiicacy.
* [calculation](./calculation) directory includes the improved method and the comparison process. The two step Bayesian method is adopted from the [Single Index Model](./reference/zjuthesis.pdf) established by Sun (2020).

