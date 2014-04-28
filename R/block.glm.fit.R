block.glm.fit <-
function(x, z, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
                          mustart = NULL, offset = rep(0, nobs), families=list(gaussian()),
                          row.families = rep(1, nobs), control = list())
{
    control <- do.call("glm.control", control)
    stopifnot(is.list(families))
    if(is.null(start)) stop("block.glm.fit only implemented for non-null start")
    if(!is.matrix(x)) x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    if(!is.matrix(z)) z <- as.matrix(z)
    znames <- dimnames(z)[[2L]]
    ynames <- if (is.matrix(y))
        rownames(y)
    else names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nblocks <- nobs / nrow(x)
    nvars <- ncol(x) * nblocks + ncol(z)
    EMPTY <- nvars == 0
    if (is.null(weights))
        weights <- rep.int(1, nobs)
    if (is.null(offset))
        offset <- rep.int(0, nobs)
    apply.by.index <- function(args,funs,indices) {
        if(is.matrix(args)) {
            stop("args must be a vector", call.=FALSE)
        }
        if(length(args) != length(indices)) {
            stop("args and indices must have same length", call.=FALSE)
        }
        if(length(funs) < max(indices) || min(indices) < 1) {
            stop("indices must be between 1 and length(funs)", call.=FALSE)
        }
        ## initialize outputs to have right mode
        outputs <- rep( funs[[ indices[1] ]](args[1]),
                       NROW(args))
        for(k in unique(indices)) {
            outputs[indices == k] <- funs[[k]](args[indices == k])
        }
        outputs
    }
    variance <- function(mu) {
        apply.by.index(mu,lapply(families, function(fam) fam$variance),row.families)
    }
    linkinv <- function(eta) {
        apply.by.index(eta,lapply(families, function(fam) fam$linkinv),row.families)
    }
    linkfun <- function(mu) {
        apply.by.index(mu,lapply(families, function(fam) fam$linkfun),row.families)
    }
    dev.resids <- function(y, mu, wt) {
        outputs <- numeric(length(y))
        for(k in unique(row.families)) {
            rowsk <- row.families==k
            outputs[rowsk] <-
                families[[k]]$dev.resids(y[rowsk],mu[rowsk],wt[rowsk])
        }
        outputs
    }
    ##apply.by.index(cbind(y,mu,wt),
    ##                   lapply(families, function(fam) {function(args) fam$dev.resids(args[,1],args[,2],args[,3])}),
    ##                   row.families)
    mu.eta <- function(eta) {
        apply.by.index(eta,lapply(families, function(fam) fam$mu.eta),row.families)
    }
    unless.null <- function(x, if.null) if (is.null(x))
        if.null
    else x

    valideta <- function(eta) {
        all(apply.by.index(eta,
                           lapply(families, function(fam) unless.null(fam$valideta, function(et) TRUE)),
                           row.families))
    }
    validmu <- function(mu) {
        all(apply.by.index(mu,
                           lapply(families, function(fam) unless.null(fam$validmu, function(m) TRUE)),
                           row.families))
    }
    ######### NO INITIALIZATION YET
    eta.from.coef <- function(coef) {
        eta <- numeric(nrow(x)*nblocks)

        jz <- nblocks*ncol(x) + 1:ncol(z)
        for(k in 1:nblocks) {
            ik <- (k-1)*nrow(x) + 1:nrow(x)
            jk <- (k-1)*ncol(x) + 1:ncol(x)
            eta[ik] <- x %*% coef[jk] + z %*% coef[jz]
        }
        eta + offset
    }
    if (EMPTY) {
        stop("you passed in an empty model", call.=FALSE)
    }
    else {
        good.resp <- !is.na(y)
        coefold <- NULL
        eta <- if (!is.null(etastart))
            etastart
        else if (!is.null(start))
            if (length(start) != nvars)
                stop(gettextf("length of 'start' should equal %d and correspond to initial coefs",
                  nvars),
                  domain = NA)
            else {
                coefold <- start
                eta.from.coef(coefold)
            }
        else linkfun(mustart)
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta)))
            stop("cannot find valid starting values: please specify some",
                call. = FALSE)
        devold <- sum(dev.resids(y, mu, weights)[good.resp])
        boundary <- conv <- FALSE
        for (iter in 1L:control$maxit) {
            ##good <- weights > 0           # not doing this because need to retain structure of y
            varmu <- variance(mu)#[good]
            ##if (any(is.na(varmu)))
            ##    stop("NAs in V(mu)")
            ##if (any(varmu == 0))
            ##    stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            ##if (any(is.na(mu.eta.val[good])))
            ##    stop("NAs in d(mu)/d(eta)")
            ##good <- (weights > 0) & (mu.eta.val != 0)
            ##if (all(!good)) {
            ##    conv <- FALSE
            ##    warning("no observations informative at iteration ",
            ##      iter)
            ##    break
            ##}
            u <- (eta - offset) + (y - mu)/mu.eta.val
            w <- sqrt((weights * mu.eta.val^2)/variance(mu))
            #ngoodobs <- as.integer(nobs - sum(!good))
            fit <- block.projection(x,z,w,u,inverse.hessian=FALSE,
                                    wt.tol=min(1e-7,control$epsilon/1000))
            if (any(!is.finite(fit$coefficients))) {
                conv <- FALSE
                warning(gettextf("non-finite coefficients at iteration %d",
                  iter), domain = NA)
                break
            }
            #if (nobs < fit$rank)
            #    stop(gettextf("X matrix has rank %d, but only %d observations",
            #      fit$rank, nobs), domain = NA)
            start <- fit$coefficients
            eta <- eta.from.coef(start)
            mu <- linkinv(eta)
            dev <- sum(dev.resids(y, mu, weights)[good.resp])
            if (control$trace)
                cat("Deviance =", dev, "Iterations -", iter,
                  "\n")
            boundary <- FALSE
            if (!is.finite(dev) || dev > devold) {
                if (is.null(coefold))
                  stop("no valid set of coefficients has been found: please supply starting values",
                    call. = FALSE)
                ##warning("step size truncated due to divergence",
                ##  call. = FALSE)
                ii <- 1
                while (!is.finite(dev) || dev > devold) {
                  if (ii > control$maxit)
                    stop("inner loop 1; cannot correct step size",
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- eta.from.coef(start)
                  mu <- linkinv(eta)
                  dev <- sum(dev.resids(y, mu, weights)[good.resp])
                }
                boundary <- TRUE
                if (control$trace)
                  cat("Step halved: new deviance =", dev, "\n")
            }
            if (!(valideta(eta) && validmu(mu))) {
                if (is.null(coefold))
                  stop("no valid set of coefficients has been found: please supply starting values",
                    call. = FALSE)
                warning("step size truncated: out of bounds",
                  call. = FALSE)
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                  if (ii > control$maxit)
                    stop("inner loop 2; cannot correct step size",
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- eta.from.coef(start)
                  mu <- linkinv(eta)
                 }
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights)[good.resp])
                if (control$trace)
                  cat("Step halved: new deviance =", dev, "\n")
            }
            if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
                conv <- TRUE
                coef <- start
                break
            }
            else {
                devold <- dev
                coef <- coefold <- start
            }
        }
        if (!conv)
            warning("glm.fit: algorithm did not converge", call. = FALSE)
        if (boundary)
            warning("glm.fit: algorithm stopped at boundary value",
                call. = FALSE)
        eps <- 10 * .Machine$double.eps
        ##if (family$family == "binomial") {
        ##    if (any(mu > 1 - eps) || any(mu < eps))
        ##        warning("glm.fit: fitted probabilities numerically 0 or 1 occurred",
        ##          call. = FALSE)
        ##}
        ##if (family$family == "poisson") {
        ##    if (any(mu < eps))
        ##        warning("glm.fit: fitted rates numerically 0 occurred",
        ##          call. = FALSE)
        ##}
        ##if (fit$rank < nvars)
        ##    coef[fit$pivot][seq.int(fit$rank + 1, nvars)] <- NA
        ##xxnames <- xnames[fit$pivot]
        residuals <- (y - mu)/mu.eta(eta)
        ##fit$qr <- as.matrix(fit$qr)
        ##nr <- min(sum(good), nvars)
        ##if (nr < nvars) {
        ##    Rmat <- diag(nvars)
        ##    Rmat[1L:nr, 1L:nvars] <- fit$qr[1L:nr, 1L:nvars]
        ##}
        ##else Rmat <- fit$qr[1L:nvars, 1L:nvars]
        ##Rmat <- as.matrix(Rmat)
        ##Rmat[row(Rmat) > col(Rmat)] <- 0
        ##names(coef) <- xnames
        ##colnames(fit$qr) <- xxnames
        ##dimnames(Rmat) <- list(xxnames, xxnames)
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    ##wt <- rep.int(0, nobs)
    wt <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    ##if (!EMPTY)
    ##    names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int("",
    ##        sum(good) - fit$rank))
    ##wtdmu <- if (intercept)
    ##    sum(weights * y)/sum(weights)
    ##else linkinv(offset)
    ##nulldev <- sum(dev.resids(y, wtdmu, weights))
    ##n.ok <- nobs - sum(weights == 0)
    ##nulldf <- n.ok - as.integer(intercept)
    ##rank <- if (EMPTY)
    ##    0
    ##else fit$rank
    ##resdf <- n.ok - rank
    ##aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
    list(coefficients = coef, residuals = residuals, fitted.values = mu,
        ##effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat,
        ##rank = rank, qr = if (!EMPTY) structure(fit[c("qr", "rank",
        ##    "qraux", "pivot", "tol")], class = "qr"), family = family,
         linear.predictors = eta, deviance = dev, #aic = aic.model,
        ##null.deviance = nulldev,
         iter = iter, weights = wt, prior.weights = weights,
         fit = fit,
        ##df.residual = resdf, df.null = nulldf,
         y = y, converged = conv,
         boundary = boundary)
}
