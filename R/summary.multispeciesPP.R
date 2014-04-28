summary.multispeciesPP <-
function(object,standardized=TRUE,...) {
    if(standardized) {
        coefs <- object$normalized.all.coef
        se <- object$normalized.std.errs
    } else {
        coefs <- object$all.coef
        se <- object$std.errs
    }
    coef.table <- data.frame(coefs,se,coefs/se,2*pnorm(-abs(coefs/se)))
    colnames(coef.table) <- c("Estimate","Std. Error","z value","Pr(>|z|)")
    tr <- list(coef.table=coef.table)
    class(tr) <- "summary.multispeciesPP"
    tr
}
