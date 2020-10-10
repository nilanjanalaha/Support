# Calculating the co-ordinate thresholding estimator

 #-------------Install packages -----------


 #library(MASS)

 #------------ function for Generating observations----------------

 # Sigma_nat _________INPUT:
 # p,q <- dimensions of alpha
 # alpha and beta
 # rho = canonical correlation

 #Sigma_ mat ________ Output
 # Gaussian X,Y with indentity covariance matrix and covariance matrix
 # as \rho\alpha\beta^T

 Sigma_mat <- function(p,q,al,be, rho)
 {
   Sx <- diag(rep(1,p), p, p)
   Sy <- diag(rep(1,q), q, q)
   Sxy <- tcrossprod(crossprod(rho*Sx, outer(al, be)), Sy)
   Syx <- t(Sxy)
   rbind(cbind(Sx, Sxy), cbind(Syx, Sy))
 }

 #-------- Coordinate thresholding of matrix------------------------
 #---------A as in Algorithm 3 of our paper-------------------------

 #ct________ Input:
 # A <- The matrix
 # t <- the threshold level
 # n <- the sample size

 #ct_________ Output
 # The thresholded matrix \eta(A)

 ct <- function(A, t, n)
 {
  apply(A, c(1,2), soft_t , t=t, n=n)
 }

 # The entrywise  soft-threshold function
 soft_t <- function(x, t, n)
 {
   sqn <- sqrt(n)
   if (x>t/sqn)
   {ans <- x-t/sqn}
   else if (x< -t/sqn)
   {ans <- x+t/sqn}
   else ans <- 0
   ans
 }
 #---------------------Estimating alpha and beta-----
 #------------ input of give_svd
 # t <- the threshold level
 # n <- the sample size
 # sxy <- the covarriance matix

 #------------ output of give_svd
# A matrix whose first column is alpha and second column is beta

 give_svd <- function(sxy, t, n)
 {
    #the thresholding step
  etam <- ct(sxy, t, n)
  temp <- svd(etam, 1, 1)
  list(temp$u, temp$v)
 }

 #--------------------- Cleaning --------------------

 #Parameters_______ tpl, see Theorem 4 of paper

 #--------Input of clean:
 # tp: patameter for threshold, see Theorem 4
 # sx, sy: sparsity
 # alpha, beta: estimates of alpha, beta
 #-------Output for clean:
 #A matrix with two coloumns, the first column gives the support for
 # X, where the second column gives the support for Y

 clean <- function(tp, syx, alpha, beta)
 {
 # Threshold for cleaning


 # Right multiplication
 al <- crossprod(syx, beta)
 #Left thresholding
 Dx <- ifelse(abs(al)>tp, 1, 0)
 #Left multiplication
 be <- crossprod(t(syx), alpha)
 # Left thresholding
 Dy <- ifelse(abs(be)>tp, 1, 0)

 #The matrix of the supports
 list(sup.x=Dx, sup.y=Dy)
 }

 #---------- Support recovery given a datset ---------
  #--------Input of c_support:
  # x <- a matrix with n rows
  # y <- a matrix with n rows
 # t:  the threshold level
 # tp: patameter for threshold, see Theorem 4


 #--------- Output
 # A matrix whose left column is the support of alpha and right
 # column is the support of beta
 ############## c_support ################################
 #' Support recovery of the first pair of canonical directions
 #'
 #' This function implements the coordinate thresholding algorithm of \emph{Laha et al. (2020)}. 
 #' to estimate the support of the first pair of canonical directions of X and Y.
 #'
 #'
 #' @details  \code{sv}: sv can be a rough estimate of the sparsity (number of non-zero elements)
 #'             of \eqn{\alpha} and \eqn{\beta}, since only the order is important (see laha et al., 2020) 
 #'             here. The returned supports of \eqn{\alpha} and \eqn{\beta} will not necessarily  match
 #'             the same sparsity levels as provided in sv. In absence of user specified sv,
 #'             it is estimated from the data using Mai and Zhang (2017)'s SCCA. 
 #' @details  \code{c1, c2}:   See Theorem 3 of Laha and Mukherjee (2020) to see how
 #'                         these tuning parameters control the thresholding of
 #'                         the covariance matrix. Larger value of c1 and c2 will
 #'                         result in a higher value of thresholding parameter, which
 #'                         will give a sparser estimator of the covariance matrix. The method is
 #'                         less sensitive to c2 than it is to c1.
 #' @details \code{tau, nl}:   The cleaning step of the co-ordinate threshold
 #'                           uses the cut-off level
 #'                           \deqn{tau \frac{\log(2s)^{0.5+nl}{s}.}
 #'                           Therefore, higher values of tau and nl result in
 #'                           a higher cut-off value, shrinking the estimated support.
 #'                           See Theorem 4 of Laha and Mukherjee (2020) for more details.
 #'
 #' @param  x  A matrix with n rows and p columns;
 #'           corresponds to the first data matrix. 
 #' @param  y A matrix with n rows and q columns; corresponds to the second data matrix.
 #' @param c1  Optional. A positive constant corresponding to the threshold level in the
 #'           co-ordinate thresholding step to estimate the covariance matrix. If 
 #'           not provided, it is calculated from the data.
 #' @param c2  Optional. A positive constant greater than c1 corresponding to the threshold level in the
 #'           co-ordinate thresholding step to estimate the covariance matrix. If 
 #'           not provided, it is calculated from the data.
 #' @param sv  Optional. A vector giving the number of non-zero elements in the canonical
 #'            covariates \eqn{\alpha} and \eqn{\beta}, respectively. See 'details'
 #'            for more information.
 #' @param tau Optional. A positive tuning parameter corresponding to the cut-off level in the
 #'           cleaning step of the co-ordinate thresholding algorithm. If 
 #'           not provided, it is calculated from the data.
 #' @param nl  Optional. A positive tuning parameter taking value in (0,0.50).  Corresponds to the cut-off level in the
 #'           cleaning step of the co-ordinate thresholding algorithm. If 
 #'           not provided, it is calculated from the data.
 #' @param Sx Optional. A \eqn{p\times p} matrix, the variance of X. Must be positive definite. If missing, estimated from the data.
 #' @param Sy Optional. A \eqn{q\times q} matrix, the variance of Y. Must be positive definite.  If missing, estimated from the data.
 #' @param is.standardize Can be either ``TRUE" or ``FALSE". Indicats hether the variables will be standardized to have mean zero.
 #'                    The default is TRUE.
 #'@return  A list of two arrays.
 #'\itemize{
 #'\item sup.x - An array with length p with binary entries. If the i-th element
 #'             is 0, it means that i is not in the support of \eqn{\alpha}.
 #'             Conversely, if the i-th element is 1, it means i is in the support
 #'             of \eqn{\alpha}.
 #'\item sup.y - An array with length p with binary entries. If the i-th element
 #'             is 0, it means that i is not in the support of \eqn{\beta}.
 #'             Conversely, if the i-th element is 1, it means i is in the support
 #'             of \eqn{\beta}.
 #' }
 #'@references Laha, N., Mukherjee, R. (2020) 
 #'\emph{Support recovery of canonical correlation analysis}. Submitted
 #' @references Mai, Q., Zhang, X. (2019) \emph{An iterative penalized least squares 
 #' approach to sparse canonical correlation analysis}, Biometrics, 75, 734-744.
 #'@author \href{https://connects.catalyst.harvard.edu/Profiles/display/Person/184207}{Nilanjana Laha}
 #' (maintainer), \email{nlaha@@hsph.harvard.edu},
 #' Rajarshi Mukherjee, \email{ram521@@mail.harvard.edu}.
 #' @seealso \code{\link{g_support}}
 #' @examples
 #' library(mvtnorm)
#' #Simulate  standard normal data matrix: first generate alpha and beta
#' p <- 500; q <- 200; al <- c(rep(1, 10), rep(0, 490));
#' be <- c(rep(0,150), rnorm(50,1))
#' 
#' #Normalize alpha and beta
#' al <- al/sqrt(sum(al^2))
#'  be <- be/sqrt(sum(be^2))
#' n <- 300; rho <- 0.5
#' 
#'  #Creating the  covariance matrix
#' Sigma_mat <- function(p,q,al,be, rho)
 #'{
 #' Sx <- diag(rep(1,p), p, p)
 #' Sy <- diag(rep(1,q), q, q)
 #' Sxy <- tcrossprod(crossprod(rho*Sx, outer(al, be)), Sy)
 #' Syx <- t(Sxy)
 #' rbind(cbind(Sx, Sxy), cbind(Syx, Sy))
 #' }
#' truesigma <-  Sigma_mat(p,q,al,be, rho)
#' 
#' #Simulating the data
 #'Z <- mvtnorm::rmvnorm(n, sigma = truesigma)
#' x <- Z[,1:p]
#' y <- Z[,(p+1):(p+q)]
#' 
#' #Support of alpha
#' which(c_support(x=x, y=y, sv=c(10, 50))$sup.x==1)
 #' @export
 c_support <- function(x, y, sv, c1, c2,  tau, nl, Sx, Sy, is.standardize)
 {
   p <- ncol(x)
   q <- ncol(y)
   
   if(missing(is.standardize)) is.standardize <- "FALSE"
   if (is.standardize=="FALSE")
   {
     x <-  scale(x, scale= FALSE)
     y <-  scale(y, scale= FALSE)
   }
   #setting default values and standardizing the data
   if(missing(c1)) c1 <- 1.2
   if(missing(c2)) c2 <- 100*c1
   if(missing(tau)) tau <- 1
   if(missing(nl)) nl <- 0.25
   
   
   # Sx and Sy
   if(missing(Sx)) 
     {x <- scale(x, center = FALSE)}
   else   x <- crossprod(t(x), solve(expm::sqrtm(Sx)))
   
   if(missing(Sy)) 
   {y <- scale(y, center = FALSE)}
   else  y <- crossprod(t(y), solve(expm::sqrtm(Sy)))
      
   
   #sparsity estimation
   if(missing(sv))
   {
   temp <- SCCA(x,y, lambda.alpha = sqrt(log(p+q)/n),
                lambda.beta = sqrt(log(p+q)/n))
   al.mai <- temp$beta
   be.mai <- temp$alpha
   
   sxy <- cov(x,y)
   tb <- crossprod(sxy, al.mai)
   ta <- crossprod(t(sxy), be.mai)
   #general method's tau
    sv <- c(length(which(al.mai!=0)), length(which(be.mai!=0)))
   }
   
   # Algorithm starts here
   s <- max(sv)
    syx <- cov(y, x)
    sxy <- t(syx)
    n <- nrow(x)

    t <- givet(s, p, q, c1, c2)
    # calculating alpha, beta
    temp <- give_svd(sxy, t, n)

    tp <- givetp(tau, s, nl)
    #cleaning alpha, beta gives the support
     clean(tp, syx, temp[[1]], temp[[2]])

 }

 #-------- function for t ---------------
 # This usees the definition of t from Theorem 3
 #----------- Input
 # s: sparsity
 # p,q
 # We take C_1=1, C_2=3 in Theorem 3
 #------------ Output
 # A value of t
 givet <- function(s, p, q, C1, C2)
 {
    hz <- log((p)/(4*s^2))
    tstar <- C1*sqrt(hz)
    if(hz>0)
    {
       if(tstar< sqrt(log(2*p))/2)
        {  return(tstar)}
       else { return(C2* tstar)}
    }
    return(0)
 suppressWarnings()
 }
 
  givetp <- function(tau, s, nl)  return(tau* (log(2*s))^(0.5+nl)/(2*s))
