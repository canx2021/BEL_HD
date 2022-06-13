# BEL_HD
Bayesian variable selection based on empirical likelihood for ultra-high dimensional data
## Description
In the semi-parametric domain, under the ultra-high dimensional setting, we propose a 
Bayesian variable selection method based on empirical likelihood, which requires only 
estimating equations and no distributional assumptions. Motivated by Chang et al. 
(2018) on doubly penalized empirical likelihood (EL), we introduce priors to regularize 
both regression parameters and Lagrange multipliers associated with the estimating 
equations. We further develop an efficient Markov chain Monte Carlo sampling 
algorithm based on the active set idea, first proposed by Friedman et al. (2007). 
The proposed method not only inherits merits from the Bayesian paradigm and shares 
the attractiveness of EL approaches, but also has superior performance in both the 
prediction and variable selection, as shown in our numerical studies.

### Depends
 R, MASS, statmod, mvtnorm, picasso, Rcpp

## Implementation
* Run BEL_HD_cpp.R 
* Run EL_DP_demo.R
* Run any code in Simulation/Empirical_Studies folder

## References
<a id="1">[1]</a> 
C. Xu, Y. Cheng, X. Wang, 
Bayesian variable selection based on empirical likelihood for ultra-high dimensional data,
Ready to submit (2021).

<a id="2">[2]</a> 
J. Chang, C. Y. Tang, T. T. Wu, et al., A new scope of penalized empirical likelihood with
high-dimensional estimating equations, Annals of Statistics 46 (2018) 3185–3216.

<a id="3">[3]</a> 
J. Friedman, T. Hastie, H. Höfling, R. Tibshirani, et al., Pathwise coordinate optimization,
Annals of applied statistics 1 (2007) 302–332.

