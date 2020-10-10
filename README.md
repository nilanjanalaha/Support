
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Support.CCA

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/nilanjanalaha/Support.CCA.svg?branch=master)](https://travis-ci.org/nilanjanalaha/Support.CCA)
<!-- badges: end -->

Support.CCA estimates the support of canonical directions using the
methods developed in Laha and Mukherjee (2020).

## Installation

Install from [GitHub](https://github.com/) with:

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
length(which(temp[1:5]==1))/5
#> [1] 0.6

#recovered support
which(temp==1)
#>  [1]   1   2   4  29  43 131 159 164 206 232 247 348 388 491
```

We see a lot of elements which are not actually in the support, are
estimated as part of the support. If we put regularization on the number
of non-zero elements (sparsity) of \(\alpha\) and \(\beta\), c\_support
will yield better estimates of the supports. We will use 7 as an
estimator of the sparsities, and pass it as the parameter sv=c(7, 7).

``` r
Sx <- diag(rep(1,p))
Sy <- diag(rep(1,q))
temp <- c_support(x,y, sv=c(7, 7))$sup.x

#Proportion of recovered support
length(which(temp[1:5]==1))/5
#> [1] 0.6

#estimated support
which(temp==1)
#> [1]   1   2   4 131 164 206 247 388
```

Thus, we see providing a regulrization parameter facilitates sparser
support recovery. However, c\_support has asymptotic garuantees only if
the covariance matrices are also provided. Still, as we see, c\_support
includes some elements, which are not in the support of \(\alpha\), to
be in the support.

``` r
Sx <- diag(rep(1,p))
Sy <- diag(rep(1,q))
temp <- c_support(x,y, sv=c(7, 7), Sx=Sx, Sy=Sy)$sup.x

#Proportion of recovered support
length(which(temp[1:5]==1))/5
#> [1] 0.6

#Esimated support
which(temp==1)
#> [1]   1   2   4 131 164 206 247 388
```

The other method g\_support also estimates the support of \(\alpha\) and
\(\beta\).

``` r
temp <- g_support(x,y)$sup.x

#Proportion of recovered support
length(which(temp[1:5]==1))/5
#> [1] 0.8
```

We see that g\_support recovers 80\(\%\) of the data. Now we will see
which elements are in estimated support.

``` r
#Estimated support 
which(temp==1)
#> [1] 1 2 4 5
```

We see that for this dataset, g\_support does not include any element
which is not in the support.

One can pass a tuning parameter Cg to g\_support. Cg controls the length
of the estimated support. If we increase Cg, then the estimated support
shrinks. A smaller Cg results in an estimated support with more zeros.
The default value is 0.50. See the documentation typing\`\` ?g.support"
in the console for more details.

``` r
temp <- g_support(x,y, Cg=.4)$sup.x

#Proportion of recovered support
length(which(temp[1:5]==1))/5
#> [1] 1

#Estimated support 
which(temp==1)
#> [1] 1 2 3 4 5
```

We see that Cg=0.40 leads to exact support recovery of \(\alpha\). One
can also specify Sx and Sy in g\_support, the variances of X and Y, but
it is not necessary.

``` r
temp <- g_support(x,y, Cg=.4, Sx=Sx, Sy=Sy)$sup.x

#Proportion of recovered support
length(which(temp[1:5]==1))/5
#> [1] 1

#Estimated support 
which(temp==1)
#> [1] 1 2 3 4 5
```
