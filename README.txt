README

Source code for manuscript "Estimating the correlation between semi-competing risk survival endpoints"
by Lexy Sorrell, Yinghui Wei, Malgorzara Wojtys and Peter Rowe.

For questions, comments or remarks about the code please contact L Sorrell (lexy.sorrell@plymouth.ac.uk).

The code has been written using R-3.6.1 (platform: Windows10, 64-bit) with 
package versions copula-0.999-19.1, mvtnorm-1.0-11, mstate-0.2.11.

To reproduce results in this paper, run the Functions.R file. This contains the likelihood functions for each of the copula models.
For the Clayton, Frank and Gumbel copulas, please see the supplementary material for derivation. For the Normal copula see:
Fu H et al. (2013) and Meyer C (2013).

To reproduce the analysis of the Amsterdam Cohort Study data set presented in the main text (Table 2, Table 5 and Figure 2) run the 'Functions.R' 
file, then 'Table 2 and Table 5.R'. The resulting tables will be saved with names: 'Table 2' and 'Table 5'.

To reproduce the simulations in the main paper (Tables 6-11) run the 'Functions.R' file and then 'Table X.R' file, where X is the Table to be 
recreated. The resulting table will be saved as a data frame with name: 'Table X'. To reduce the running time of the simulations, 'runs', describing 
the number of times each simulation will be run, can be reduced from 1000. 

References:

Fu H, et al. 2013. Joint modeling of progression‐free survival and overall survival by a Bayesian normal induced copula 
estimation model. Statistics in medicine, 32(2), pp.240-254.
Meyer C. 2013. The bivariate normal copula. Communications in Statistics-Theory and Methods, 42(13), pp.2402-2422.
