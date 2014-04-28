block.projection <-
function(x,z,w,y,inverse.hessian=FALSE,wt.tol=1E-30) {
    n.blocks <- length(w)/nrow(x)
    block.size <- nrow(x)
    resid <- numeric(length(y))
    fitted <- numeric(length(y))
    coef.yz.on.x <- matrix(NA,n.blocks*ncol(x),ncol(z)+1)
    zTyz <- 0
    sq.std.errs <- rep(NA,n.blocks*ncol(x)+ncol(z))

    if(inverse.hessian) {
        inv.hess <- matrix(0,n.blocks * ncol(x) + ncol(z),n.blocks * ncol(x) + ncol(z))
    }
    for(b in 1:n.blocks) {
        ib <- (b-1)*block.size + (1:block.size)
        which.good <- (w[ib] > wt.tol) & is.finite(y[ib])
        #xb <- x[which.good,,drop=FALSE]
        #zb <- z[which.good,,drop=FALSE]
        ib <- ib[which.good]
        jb <- (b-1)*ncol(x) + (1:ncol(x))
        wb <- w[ib]
        qrb <- qr(x[which.good,,drop=FALSE]* wb)
        coef.yz.on.x[jb,] <- qr.coef(qrb,cbind(y[ib] * wb,z[which.good,,drop=FALSE] * wb))
        yz.x <- qr.resid(qrb,cbind(y[ib] * wb,z[which.good,,drop=FALSE] * wb))
        ##yz.x <- cbind(y[ib],zb) * wb - (xb * wb) %*% coef.yz.on.x[jb,]
        zTyz <- zTyz + t(yz.x[,-1,drop=FALSE]) %*% yz.x
        if(ncol(x) < length(ib)) {
            if(inverse.hessian) {
                inv.hess[jb,jb] <- chol2inv(qr.R(qrb))
            }
            sq.std.errs[jb] <- diag(chol2inv(qr.R(qrb)))
        }
    }
    ## now, compute coefficients
    qr.zTz <- qr(zTyz[,-1])
    coefZ <- qr.solve(qr.zTz,zTyz[,1,drop=FALSE])
    coefX <- coef.yz.on.x[,1,drop=FALSE] - coef.yz.on.x[,-1,drop=FALSE] %*% coefZ

    ## and residuals
    for(b in 1:n.blocks) {
        ib <- (b-1)*block.size + (1:block.size)
        jb <- (b-1)*ncol(x) + (1:ncol(x))
        resid[ib] <- y[ib] - x %*% coefX[jb] - z %*% coefZ
        fitted[ib] <- x %*% coefX[jb] + z %*% coefZ
    }

    ## and diagonal of inverse Hessian
    sq.std.errs[n.blocks*ncol(x) + (1:ncol(z))] <- diag(qr.solve(qr.zTz))
    sq.std.errs[1:(n.blocks*ncol(x))] <- sq.std.errs[1:(n.blocks*ncol(x))] +
        colSums(t(coef.yz.on.x[,-1,drop=FALSE]) * qr.solve(qr.zTz,t(coef.yz.on.x[,-1,drop=FALSE])))

    if(inverse.hessian) {
        inv.hess <- inv.hess + rbind(-coef.yz.on.x[,-1,drop=FALSE],diag(ncol(z))) %*%
            qr.solve(qr.zTz,cbind(-t(coef.yz.on.x[,-1,drop=FALSE]),diag(ncol(z))))
        return(list(resid=resid,fitted=y-resid,coefficients=c(coefX,coefZ),
                    std.errs=sqrt(sq.std.errs),inverse.hessian=inv.hess))
    }
    else {
        return(list(resid=resid,fitted=fitted,coefficients=c(coefX,coefZ),
                    std.errs=sqrt(sq.std.errs)))
    }
}
