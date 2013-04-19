print.kripp.boot <- function(boot.stats) {
    cat("Bootstrapped Krippendorff's Alpha","\n\n")
    cat('Alpha Levels for', boot.stats$raters,'raters,', 
        boot.stats$iter, ' iterations', '\n')
    cat('     ', boot.stats$probs[1],'= ', boot.stats$lower, '\n')
    cat(' Mean alpha = ', boot.stats$mean.alpha, '\n')
    cat('     ', boot.stats$probs[2],'= ', boot.stats$upper, '\n')
}
