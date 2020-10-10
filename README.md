
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Support.CCA

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/nilanjanalaha/Support.CCA.svg?branch=master)](https://travis-ci.org/nilanjanalaha/Support.CCA)
<!-- badges: end -->

Support.CCA estimates the support of canonical directions using the
methods developed in Laha and Mukherjee (2020).

## Installation

You can install the released version of Support.CCA from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("Support.CCA")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("nilanjanalaha/Support.CCA")
```

## Example

We will generate a toy dataset. First let us generate \(\alpha\) and
\(\beta\).

``` r
library(mvtnorm)
#> Warning: package 'mvtnorm' was built under R version 3.6.2
#Simulate  standard normal data matrix: first generate alpha and beta
p <- 500; q <- 200; al <- rep(0, p); be <- rep(0,q)
al[1:5] <- 1
be[20:25] <- 1

#Standardize alpha and beta
al <- al/sqrt(sum(al^2))
be <- be/sqrt(sum(be^2))
n <- 800; rho <- 0.5
```

Note that the support of \(\alpha\) and \(\beta\) are
\[\text{Sup}(\alpha)=\{1,\ldots,5\},\quad \text{Sup}(\beta)=\{20,\ldots,25\}.\]

Now we will generate two gaussian data matrices X and Y with rank one
covariance matrix, and canonical covariates are \(\alpha\) and
\(\beta\).

``` r

#Creating the  covariance matrix
Sigma_mat <- function(p,q,al,be, rho)
{
  Sx <- diag(rep(1,p), p, p)
  Sy <- diag(rep(1,q), q, q)
  Sxy <- tcrossprod(crossprod(rho*Sx, outer(al, be)), Sy)
  Syx <- t(Sxy)
  rbind(cbind(Sx, Sxy), cbind(Syx, Sy))
}
truesigma <-  Sigma_mat(p,q,al,be, rho)

#Simulating the data
set.seed(4)
Z <- mvtnorm::rmvnorm(n, sigma = truesigma)
x <- Z[,1:p]
y <- Z[,(p+1):(p+q)]
```

The function c\_support provides supports of the leading pair of
canonical covariates \(\alpha\) and \(\beta\). It uses a coordinate
thresholding algorithm (Laha and Mukherjee., 2020).

``` r
library(Support.CCA)
temp <- c_support(x,y)$sup.x

#Proportion of recovered support
which(temp[1:5]==1)/5
#> [1] 0.4
```

If we can provide some estimattors of the number of non-zero elements
(sparsity) of \(\alpha\) and \(\beta\), c\_support will yield better
estimates of the supports.

``` r
Sx <- diag(rep(1,p))
Sy <- diag(rep(1,q))
temp <- c_support(x,y, sv=c(7, 7))$sup.x

#Proportion of recovered support
length(which(temp[1:5]==1))/5
#> [1] 0.4
```

However, c\_support has asymptotic garuantees only if the covariance
matrices are also provided.

``` r
Sx <- diag(rep(1,p))
Sy <- diag(rep(1,q))
temp <- c_support(x,y, sv=c(7, 7), Sx=Sx, Sy=Sy)$sup.x

#Proportion of recovered support
length(which(temp[1:5]==1))/5
#> [1] 0.6
```
