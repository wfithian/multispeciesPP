\name{block.glm.fit}
\alias{block.glm.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to fit a generalized linear model with a particular block
  structure for the design matrix.
}
\description{
  This function fits a GLM where each row has a different family and the
  block structure is like that described in the Proportional Observer
  Bias paper.
}
\usage{
block.glm.fit(x, z, y, weights = rep(1, nobs), start = NULL, etastart = NULL, mustart = NULL, offset = rep(0, nobs), families = list(gaussian()), row.families = rep(1, nobs), control = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    The X matrix lying on the diagonal blocks of the design.
}
  \item{z}{
    The Z matrix lying on the blocks in the right column of
    the design.
}
  \item{y}{
    Responses
}
  \item{weights}{
    Case weights
}
  \item{start}{
    Starting coefficients, which must be given
}
  \item{etastart}{
    Starting linear predictors (not implemented)
}
  \item{mustart}{
    Starting fits for the mean of y
}
  \item{offset}{
    Offsets for the generalized linear model
}
  \item{families}{
    List of all the families that occur for any observation.
}
  \item{row.families}{
    List of integers whose length is \code{length(y)} between 1 and
    \code{length(families)}.  If \code{row.families[i]} is 3, then row
    \code{i} uses the third family in the list of families.
}
  \item{control}{
    Similar to \code{control} for \code{glm.fit}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  William Fithian
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, z, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
    mustart = NULL, offset = rep(0, nobs), families = list(gaussian()),
    row.families = rep(1, nobs), control = list())
{
    control <- do.call("glm.control", control)
    stopifnot(is.list(families))
    if (is.null(start))
        stop("block.glm.fit only implemented for non-null start")
    if (!is.matrix(x))
        x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    if (!is.matrix(z))
        z <- as.matrix(z)
    znames <- dimnames(z)[[2L]]
    ynames <- if (is.matrix(y))
        rownames(y)
    else names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nblocks <- nobs/nrow(x)
    nvars <- ncol(x) * nblocks + ncol(z)
    EMPTY <- nvars == 0
    if (is.null(weights))
        weights <- rep.int(1, nobs)
    if (is.null(offset))
        offset <- rep.int(0, nobs)
    apply.by.index <- function(args, funs, indices) {
        if (is.matrix(args)) {
            stop("args must be a vector", call. = FALSE)
        }
        if (length(args) != length(indices)) {
            stop("args and indices must have same length", call. = FALSE)
        }
        if (length(funs) < max(indices) || min(indices) < 1) {
            stop("indices must be between 1 and length(funs)",
                call. = FALSE)
        }
        outputs <- rep(funs[[indices[1]]](args[1]), NROW(args))
        for (k in unique(indices)) {
            outputs[indices == k] <- funs[[k]](args[indices ==
                k])
        }
        outputs
    }
    variance <- function(mu) {
        apply.by.index(mu, lapply(families, function(fam) fam$variance),
            row.families)
    }
    linkinv <- function(eta) {
        apply.by.index(eta, lapply(families, function(fam) fam$linkinv),
            row.families)
    }
    linkfun <- function(mu) {
        apply.by.index(mu, lapply(families, function(fam) fam$linkfun),
            row.families)
    }
    dev.resids <- function(y, mu, wt) {
        outputs <- numeric(length(y))
        for (k in unique(row.families)) {
            rowsk <- row.families == k
            outputs[rowsk] <- families[[k]]$dev.resids(y[rowsk],
                mu[rowsk], wt[rowsk])
        }
        outputs
    }
    mu.eta <- function(eta) {
        apply.by.index(eta, lapply(families, function(fam) fam$mu.eta),
            row.families)
    }
    unless.null <- function(x, if.null) if (is.null(x))
        if.null
    else x
    valideta <- function(eta) {
        all(apply.by.index(eta, lapply(families, function(fam) unless.null(fam$valideta,
            function(et) TRUE)), row.families))
    }
    validmu <- function(mu) {
        all(apply.by.index(mu, lapply(families, function(fam) unless.null(fam$validmu,
            function(m) TRUE)), row.families))
    }
    eta.from.coef <- function(coef) {
        eta <- numeric(nrow(x) * nblocks)
        jz <- nblocks * ncol(x) + 1:ncol(z)
        for (k in 1:nblocks) {
            ik <- (k - 1) * nrow(x) + 1:nrow(x)
            jk <- (k - 1) * ncol(x) + 1:ncol(x)
            eta[ik] <- x \%*\% coef[jk] + z \%*\% coef[jz]
        }
        eta + offset
    }
    if (EMPTY) {
        stop("you passed in an empty model", call. = FALSE)
    }
    else {
        good.resp <- !is.na(y)
        coefold <- NULL
        eta <- if (!is.null(etastart))
            etastart
        else if (!is.null(start))
            if (length(start) != nvars)
                stop(gettextf("length of 'start' should equal \%d and correspond to initial coefs",
                  nvars), domain = NA)
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
            varmu <- variance(mu)
            mu.eta.val <- mu.eta(eta)
            u <- (eta - offset) + (y - mu)/mu.eta.val
            w <- sqrt((weights * mu.eta.val^2)/variance(mu))
            fit <- block.projection(x, z, w, u, inverse.hessian = FALSE,
                wt.tol = min(1e-07, control$epsilon/1000))
            if (any(!is.finite(fit$coefficients))) {
                conv <- FALSE
                warning(gettextf("non-finite coefficients at iteration \%d",
                  iter), domain = NA)
                break
            }
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
        residuals <- (y - mu)/mu.eta(eta)
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    list(coefficients = coef, residuals = residuals, fitted.values = mu,
        linear.predictors = eta, deviance = dev, iter = iter,
        weights = wt, prior.weights = weights, fit = fit, y = y,
        converged = conv, boundary = boundary)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
