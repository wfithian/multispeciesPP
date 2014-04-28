\name{block.projection}
\alias{block.projection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Weighted least-squares solve with a certain block structure
}
\description{
  The function is called by \code{block.glm.fit} at each iteration to do
  a computationally efficient step of iterative reweighted least squares.
}
\usage{
block.projection(x, z, w, y, inverse.hessian = FALSE, wt.tol = 1e-30)
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
  \item{w}{
    Vector of weights.
}
  \item{y}{
    Vector to be projected
}
\item{inverse.hessian}{
  Should be FALSE
}
\item{wt.tol}{
  Used in internal computations.  Should be small
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
function (x, z, w, y, inverse.hessian = FALSE, wt.tol = 1e-30)
{
    n.blocks <- length(w)/nrow(x)
    block.size <- nrow(x)
    resid <- numeric(length(y))
    fitted <- numeric(length(y))
    coef.yz.on.x <- matrix(NA, n.blocks * ncol(x), ncol(z) +
        1)
    zTyz <- 0
    sq.std.errs <- rep(NA, n.blocks * ncol(x) + ncol(z))
    if (inverse.hessian) {
        inv.hess <- matrix(0, n.blocks * ncol(x) + ncol(z), n.blocks *
            ncol(x) + ncol(z))
    }
    for (b in 1:n.blocks) {
        ib <- (b - 1) * block.size + (1:block.size)
        which.good <- (w[ib] > wt.tol) & is.finite(y[ib])
        ib <- ib[which.good]
        jb <- (b - 1) * ncol(x) + (1:ncol(x))
        wb <- w[ib]
        qrb <- qr(x[which.good, , drop = FALSE] * wb)
        coef.yz.on.x[jb, ] <- qr.coef(qrb, cbind(y[ib] * wb,
            z[which.good, , drop = FALSE] * wb))
        yz.x <- qr.resid(qrb, cbind(y[ib] * wb, z[which.good,
            , drop = FALSE] * wb))
        zTyz <- zTyz + t(yz.x[, -1, drop = FALSE]) \%*\% yz.x
        if (ncol(x) < length(ib)) {
            if (inverse.hessian) {
                inv.hess[jb, jb] <- chol2inv(qr.R(qrb))
            }
            sq.std.errs[jb] <- diag(chol2inv(qr.R(qrb)))
        }
    }
    qr.zTz <- qr(zTyz[, -1])
    coefZ <- qr.solve(qr.zTz, zTyz[, 1, drop = FALSE])
    coefX <- coef.yz.on.x[, 1, drop = FALSE] - coef.yz.on.x[,
        -1, drop = FALSE] \%*\% coefZ
    for (b in 1:n.blocks) {
        ib <- (b - 1) * block.size + (1:block.size)
        jb <- (b - 1) * ncol(x) + (1:ncol(x))
        resid[ib] <- y[ib] - x \%*\% coefX[jb] - z \%*\% coefZ
        fitted[ib] <- x \%*\% coefX[jb] + z \%*\% coefZ
    }
    sq.std.errs[n.blocks * ncol(x) + (1:ncol(z))] <- diag(qr.solve(qr.zTz))
    sq.std.errs[1:(n.blocks * ncol(x))] <- sq.std.errs[1:(n.blocks *
        ncol(x))] + colSums(t(coef.yz.on.x[, -1, drop = FALSE]) *
        qr.solve(qr.zTz, t(coef.yz.on.x[, -1, drop = FALSE])))
    if (inverse.hessian) {
        inv.hess <- inv.hess + rbind(-coef.yz.on.x[, -1, drop = FALSE],
            diag(ncol(z))) \%*\% qr.solve(qr.zTz, cbind(-t(coef.yz.on.x[,
            -1, drop = FALSE]), diag(ncol(z))))
        return(list(resid = resid, fitted = y - resid, coefficients = c(coefX,
            coefZ), std.errs = sqrt(sq.std.errs), inverse.hessian = inv.hess))
    }
    else {
        return(list(resid = resid, fitted = fitted, coefficients = c(coefX,
            coefZ), std.errs = sqrt(sq.std.errs)))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
