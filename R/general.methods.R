#' Support recovery of the first pair of canonical directions: unknown sparsity
#'
#' Suppose \eqn{\alpha} and \eqn{\beta} are the first pair of canonical covariates
#' corresponding to random vectors X and Y.
#' This function uses Mai and Zhang (2017)'s SCCA to estimate \eqn{\alpha} and
#' \eqn{\beta}. The estimation step is followed by a cleaning step, which
#' results in  refined estimates of the supports of \eqn{\alpha} and \eqn{\beta}.
#'
#'
#' @details  \code{sv}: sv can be a rough estimate of the sparsity (number of non-zero elements)
#'             of \eqn{\alpha} and \eqn{\beta}. The returned supports of \eqn{\alpha} and \eqn{\beta} will not necessarily  match
#'             the same sparsity levels as provided in sv. The sparsity is only required for
#'             setting a good cut-off in the cleaning step (see Laha and Mukherjee, 2020).
#'              In absence of user specified sv, it is estimated from the data using Mai and Zhang (2017)'s SCCA. 
#' @details  \code{Cg}:   The cut-off in the cleaning step is set by
#'                       \deqn{C = Cg\log(2s)/s.}
#'                        See Theorem 2 of Laha and Mukherjee (2020) for more details. Therefore, a
#'                        higher value of Cg will lead to  a smaller estimated support.
#' @details \code{Sx, Sy}:   The method requires an estimator of the inverse covariance
#'                           matrices, i.e. the precision matrices of X and Y. Unless provided by
#'                           the user, estimated using \code{\link[clime]{clime}} method of Cai et al. (2011).
#'
#' @param  x  A matrix with n rows and p columns;
#'           corresponds to the first data matrix. 
#' @param  y A matrix with n rows and q columns; corresponds to the second data matrix.
#' @param Cg  Optional. A positive constant corresponding to the threshold level in the
#'           co-ordinate thresholding step to estimate the covariance matrix. If 
#'           not provided, it is calculated from the data.
#' @param sv  Optional. A vector giving the  number of non-zero elements in the canonical
#'            covariates \eqn{\alpha} and \eqn{\beta}, respectively. See 'details'
#'            for more information.
#' @param Sx Optional. A \eqn{p\times p} positive definite matrix, the variance of X. Must be positive definite. If missing, estimated from the data.
#' @param Sy Optional. A \eqn{q\times q} positive definite matrix, the variance of Y. Must be positive definite.  If missing, estimated from the data.
#' @param is.standardize Can be either ``TRUE" or ``FALSE". Indicates whether the variables will be standardized to have mean zero.
#'                    The default is FALSE.
#'@return  A list of two arrays
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
#'\emph{Support recovery of canonical correlation analysis}. Submitted.
#' @references Mai, Q., Zhang, X. (2019) \emph{An iterative penalized least squares 
#' approach to sparse canonical correlation analysis}, Biometrics, 75, 734-744.
#' @references Cai, T., Liu, W., and Luo, X. (2012) \emph{A Constrained \eqn{l_1}
#'  Minimization Approach to Sparse Precision Matrix Estimation}, JASA, 106, 594-607.
#'@author \href{https://connects.catalyst.harvard.edu/Profiles/display/Person/184207}{Nilanjana Laha}
#'(maintainer), \email{nlaha@@hsph.harvard.edu},
#' Rajarshi Mukherjee, \email{ram521@@mail.harvard.edu}.
#' @seealso \code{\link{c_support}}
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
#' #Support of beta
#' which(g_support(x,y)$sup.y==1)
#' @export
 g_support <- function(x, y, Cg, sv, Sx, Sy, is.standardize)
 {
   p <- ncol(x)
   q <- ncol(y)
   
   n <- nrow(x)
   
   #setting default values and standardizing the data
   if(missing(Cg)) Cg <- 1
   if(missing(is.standardize)) is.standardize <- "FALSE"
  rho <- sqrt(log(p+q)/n)
     #CLIME
   if(missing(Sx)) 
   {Phix <- glasso::glasso(var(x), rho)$wi} else  { Phix <- solve(Sx)}
   
   if(missing(Sy)) 
   {Phiy <- glasso::glasso(var(y), rho)$wi} else  { Phiy <- solve(Sy)}
   
  
   #Standardize
   if (is.standardize=="FALSE")
   {
     x <-  scale(x, scale=FALSE)
     y <-  scale(y, scale=FALSE)
   }
     
  
  
   # Finding the canonical covariates
   temp <- SCCA(x,y, lambda.alpha = sqrt(log(p+q)/n),
                lambda.beta = sqrt(log(p+q)/n))
   al.mai <- temp$beta
   be.mai <- temp$alpha
   
   sxy <- cov(x,y)
   tb <- crossprod(Phix, sxy)
   tb <- crossprod(tb, al.mai)
   
   syx <- t(sxy)
   ta <- crossprod(Phiy, syx)
   ta <- crossprod(ta, be.mai)
   
   #general method's tau
   if(missing(sv)) sv <- c(length(which(al.mai!=0)), length(which(be.mai!=0)))
   s <- max(sv)
   
   #Calculating cut-off
   gen.tau <- Cg*log(2*s)/s #if s=1, this becomes 0
   
   al.est <- ifelse(abs(ta)>gen.tau, 1, 0)
   be.est <- ifelse(abs(tb)>gen.tau, 1, 0)
   
   list(sup.x = al.est, sup.y = be.est)
 } 
 
 #' Support recovery of the first pair of canonical directions using Mai and Zhang (2017)
 #'
 #' Suppose \eqn{\alpha} and \eqn{\beta} are the first pair of canonical covariates
 #' corresponding to random vectors X and Y.
 #' This function uses Mai and Zhang (2017)'s SCCA to estimate \eqn{\alpha} and
 #' \eqn{\beta}, and outputs their supports as estimates of the supports of \eqn{\alpha}
 #' and \eqn{\beta}. The code in the authors' website is used to implement the method.
 #' This gives an initial estimator of the support, which can be refined 
 #' by a cleaning step (see Laha and Mukherjee, 2020). The refined estimator
 #' of the support is given by the function \code{\link{g_support}}.
 #'
 #' @param  x  A matrix with n rows and p columns;
 #'           corresponds to the first data matrix. 
 #' @param  y A matrix with n rows and q columns; corresponds to the second data matrix.
  #'@return  A list of two arrays
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
 #'\emph{Support recovery of canonical correlation analysis}. Submitted.
 #' @references Mai, Q., Zhang, X. (2019) \emph{An iterative penalized least squares 
 #' approach to sparse canonical correlation analysis}, Biometrics, 75, 734-744.
 #'@author \href{https://connects.catalyst.harvard.edu/Profiles/display/Person/184207}{Nilanjana Laha}
 #'(maintainer), \email{nlaha@@hsph.harvard.edu},
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
 #' #Support of beta
 #' which(m_support(x,y)$sup.y==1)
 #' @export
 
 m_support <- function(x, y )
 {
    p <- ncol(x)
    q <- ncol(y)
    n <- nrow(x)
    
      # Finding the canonical covariates
    temp <- SCCA(x,y, lambda.alpha = sqrt(log(p+q)/n),
                 lambda.beta = sqrt(log(p+q)/n))
   ta <- temp$beta
   tb <- temp$alpha
    
    
    al.est <- ifelse(abs(ta)>0, 1, 0)
    be.est <- ifelse(abs(tb)>0, 1, 0)
    
    list(sup.x = al.est, sup.y = be.est)
 } 
