#' Bootstrap Krippendorff's Alpha
#' 
#' This function implements Prof. Klaus Krippendorff's algorithm for bootstrapping 
#' the Krippendorff's alpha coefficient. It computes confidence values (reliability 
#' estimates) for the given probabilities. 
#'  
#' @usage 
#' kripp.boot(x, iter = 2000, probs = c(.025, .975), 
#'            method = c("nominal", "ordinal", "interval", "ratio"))
#'
#' @param x is a matrix with rows (in R: observations) corresponding to judges and
#' columns (in R: variables) corresponding to rated objects. Should be numeric, with NAs for missing data.
#' 
#' @param iter the number of iterations for bootstrapping.
#' 
#' @param probs a vector of probabilities for which confidence values are computed.
#' 
#' @param method the metric used to calculate the difference function.
#' "nominal", "ordinal", "interval", and "ratio" are currently implemented.
#'
#' @return A list containing the following components:
#' @return \item{$mean.alpha}{the mean value of all bootstrapped alpha replicates}
#' @return \item{$alpha}{a vector of bootstrapped alphas}
#' @return \item{$upper}{upper alpha value for given probabilities}
#' @return \item{$lower}{lower alpha value for given probabilities}
#' @return \item{$raters}{number of raters used in calculating alpha}
#' @return \item{$iter}{number of bootstrap replications}
#' @return \item{$probs}{vector of probabilities used}
#' @return \item{$size}{number of items used in calculating alpha}
#' 
#' @author Polina Proutskova (proutskova@googlemail.com) 
#' @author Mike Gruszczynski (mikewgruz@gmail.com) 
#' @references Krippendorff, K. (2011). Computing Krippendorff's Alpha-Reliability. Retrieved from http://repository.upenn.edu/asc_papers/43
#' @seealso Algorithm for bootstrapping a distribution of alpha (http://www.afhayes.com/public/alphaboot.pdf)
#' @seealso Andrew F. Hayes's SPSS code (http://www.afhayes.com/public/kalpha.sps)
#' @seealso Krippendorff, K. (2012). Content analysis: An introduction to its methodology. Sage.
#' 
#' @examples 
#' # Krippendorff's "C" data (2011, 2)
#' nmm<-matrix(c(1,1,NA,1,2,2,3,2,3,3,3,3,3,3,3,3,2,2,2,2,1,2,3,4,4,4,4,4,
#' 1,1,2,1,2,2,2,2,NA,5,5,5,NA,NA,1,1,NA,NA,3,NA),nrow=4)
#'
#' # assume default nominal classification with 2000 replicates
#' kripp.boot(nmm)
#'
#' # nominal classification with 5000 replicates
#' kripp.boot(nmm, iter=5000)
#'
#' # ordinal classification with 5000 replicates
#' kripp.boot(nmm, iter=5000, method="ordinal")
#' @export

kripp.boot <- function(x, iter=2000, probs=c(.025,.975), method = c("nominal", "ordinal", "interval", "ratio")) {
  if(missing(x))
    stop("missing data: kripp.alpha(x,method=c(\"nominal\",\"ordinal\",\"interval\",\"ratio\"))\n", 
         "\twhere x is a classifier by object matrix of classifications or scores\n")
  method <- match.arg(method)  
  
  # dimensions of x
  dimx<-dim(x)
  njudges<-dimx[1]
  nvars<-dimx[2]
  
  # levels of measurment in x
  levx<-(levels(as.factor(x)))
  nvalx<-length(levx)
  if(match(method[1],"confidence",0)) {
    ratingsx<-as.complex(levx)
  }
  else ratingsx<-as.numeric(levx)
  
  # mu - the number of judgements for each variable (column of x) 
  vn<-function(datavec) sum(!is.na(datavec))
  mu<-apply(x, 2, vn)
  
  # observed coincidence matrix
  coincidence.matrix<-function(x) {
    coincm<-matrix(rep(0, nvalx * nvalx), nrow = nvalx)    
    
    # construct the matrix
    for(col in 1:nvars) {  #for each variable
      for(i1 in 1:(njudges - 1)) {
        for(i2 in (i1 + 1):njudges) { #for each pair of judges
          if(!is.na(x[i1, col]) && !is.na(x[i2, col])) {
            judgement1<-which(levx == x[i1, col])
            judgement2<-which(levx == x[i2, col])
            coincm[judgement1, judgement2]<-coincm[judgement1,judgement2] + (1 + (judgement1 == judgement2))/(mu[col]-1)
            if(judgement1 != judgement2) coincm[judgement2,judgement1]<-coincm[judgement1,judgement2]
          }
        }
      } #for each pair of judges
    }#for each variable
    
    return(coincm)
  } #coincidence.matrix
  
  coinc <- coincidence.matrix(x)
  nc<-apply(coinc,1,sum)
  # ncnc<-sum(nc * (nc - 1))
  n.. <- sum(apply(coinc, 2, sum))
  dim_coinc<-dim(coinc)
  tri_coinc<-as.vector(coinc[upper.tri(coinc, diag=TRUE)])
  
  # expected coincidence matrix
  expected<-matrix(rep(0, nvalx * nvalx), nrow = nvalx) 
  for(k in 1:nvalx) {
    for(c in 1:nvalx) {
      expected[c,k] <- nc[c] * (nc[k] - (c==k)) / (n..-1)
    }
  }
  
  # difference function
  delta<-matrix(rep(0, nvalx * nvalx), nrow = nvalx)
  for(k in 1:nvalx) {
    for(c in 1:nvalx) {
      if(match(method[1],"nominal",0)) delta[c,k] <- as.numeric(c!=k)
      if(match(method[1],"ordinal",0)) {
        delta[c,k]<-0
        if (c<=k) for(g in (c):(k)) delta[c,k]<-delta[c,k] + nc[g]
        else for(g in (k):(c)) delta[c,k]<-delta[c,k] + nc[g]      
        delta[c,k]<-delta[c,k]-(nc[k]+nc[c])/2
        delta[c,k]<-delta[c,k]^2
      }
      if(match(method[1],"interval",0)) delta[c,k]<-(ratingsx[c]-ratingsx[k])^2
      if(match(method[1],"ratio",0)) delta[c,k]<-(ratingsx[c]-ratingsx[k])^2/(ratingsx[c]+ratingsx[k])^2
      if(match(method[1],"confidence",0)) {
        delta[c,k]<-((Re(ratingsx[c])-Re(ratingsx[k]))^2)*Im(ratingsx[c])*Im(ratingsx[k])
      }
    }
  } # delta
  
  # observed disagreement
  observed_dis <- sum(coinc*delta)
  # expected disagreement
  expected_dis <- sum(expected*delta)  # ??? /n..
  # alpha
  alpha <- 1 - observed_dis/expected_dis
  
  # probability function for bootstrapping pmat
  pcoinc <- 2*(coinc/n..)
  diag(pcoinc)<-diag(coinc)/n..
  psum<-0
  ck<-1
  pmat <- matrix(rep(0, 2*length(tri_coinc)), nrow = 2)
  for(k in 1:nvalx) {
    for(c in k:nvalx) {
      psum <- psum+pcoinc[c,k]
      pmat[1,ck] <- psum
      pmat[2,ck] <- delta[c,k]
      ck <- ck+1
    }
  }
  pmat <- matrix(c(c(0,0), pmat), nrow=2)  # add a (0,0) column at the beginning
  
  # bootstrapping
  q <- sum(coinc!=0)
  bootM <- min(25*q, (((njudges-1)*n..)/2))
  
  btalpha = rep(NA, iter)
  numone <- 0
  
  for (z in 1:iter) {   
    # constract distribution
    rand <- stats::runif(bootM)
    numsum <- 0
    for (i in 1:length(rand)) 
      for (j in 2:ncol(pmat)) 
        if (rand[i] <= pmat[1,j])
          if (rand[i] >= pmat[1,j-1])
            numsum <- numsum + pmat[2,j]
    
    # compute alpha
    alpha <- 1 - (numsum/((expected_dis/n..)*bootM))  # 
    ## not quite sure why expected_dis should be divided by n..
    ## Krippendorff's algorithm says: 1 - SUM / (M * De)  ---- no norming
    ## Hayes's code says: compute expdis=csum(rsum((expect&*delta)))/n ----- normed
    ## Unnormed version produces apparently wrong results
    if (alpha < -1) alpha <- -1
    if (alpha == 1 &  sum(diag(coinc) != 0)==1 )
      alpha <- 0
    if (alpha == 1 &  sum(diag(coinc) != 0) >1 )
      numone <- numone+1
    
    btalpha[z]<-alpha
    
    #debug
    #     if(alpha<0.2) {
    #       cat('alpha=',alpha,', c=',c,', k=',k,'\n')
    #     }
    
  } #for (z in 1:iter+1)
  
  # Correct the distribution for situations in which the lack of variation 
  # should cause alpha to be indeterminate ( alpha = 1 – 0/0 )
  nx <- round(iter*sum((diag(coinc)/n..)^bootM))
  lim <- min(nx, numone)
  chk <- 0
  for (i in 1:length(btalpha)) {
    if (btalpha[i] == 1 & chk < lim) {
      btalpha[i] <- 0
      chk <- chk + 1
    }
  }
  
  # sort the bootstrap estimates
  btalpha <- sort(btalpha)
  
  # results
  lowlim <- max(trunc(min(probs) * iter),1)
  highlim <- min(trunc(max(probs) * iter) + 1,length(btalpha))
  boot.stats <- list(mean.alpha=mean(btalpha, na.rm=TRUE), 
                     upper=btalpha[highlim], 
                     lower=btalpha[lowlim], 
                     alphas=btalpha,
                     raters=njudges,
                     iter=iter,
                     probs=probs,
                     size=ncol(x))
  class(boot.stats) <- 'kripp.boot'
  return(boot.stats)
}

#' @export
print.kripp.boot <- function(x, ...) {
  cat('\n',"Bootstrapped Krippendorff's Alpha","\n")
  cat('Alpha Levels for', x$raters,'raters,', 
      x$iter, ' iterations', '\n')
  cat('     ', x$probs[1],'= ', x$lower, '\n')
  cat(' Mean alpha = ', x$mean.alpha, '\n')
  cat('     ', x$probs[2],'= ', x$upper, '\n')
}
