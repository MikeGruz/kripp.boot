
kripp.boot <- function(x, raters='rows', probs=c(.025,.975), iter=100, method='nominal') {
    
    if (!is.numeric(as.matrix(x))) stop(x, ' contains non-numeric cells')

    alphas <- numeric(iter)
    
    if (raters == 'cols') x <- t(x)
        
    for (i in 1:iter) {
        alphas[i] <- kripp.alpha(x[,sample(ncol(x),
                                           size=ncol(x), 
                                           replace=TRUE)], 
                                 method=method)$value
    }
    
    kripp.ci <- quantile(alphas, probs=probs, na.rm=TRUE)
    boot.stats <- list(mean.alpha=mean(alphas, na.rm=TRUE), 
                       upper=kripp.ci[2], 
                       lower=kripp.ci[1], 
                       alphas=alphas,
                       raters=nrow(x),
                       iter=iter,
                       probs=probs,
                       size=ncol(x))
    class(boot.stats) <- 'kripp.boot'
    return(boot.stats)
}
