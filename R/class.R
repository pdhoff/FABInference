#' @title Summarizing Linear Model Fits with FAB Inference
#' 
#' @description \code{summary} method for class 'lmFAB'
#' 
#' @param object an object of class \code{lmFAB}
#' @param correlation see \code{summary.lm}
#' @param symbolic.cor see \code{summary.lm}
#' @param ... see \code{summary.lm}
#'
#' @details A mod of \code{summary.lm} that shows FAB p-values in table
#'
#' @import stats
#' 
#' @export 
summary.lmFAB<-function (object, correlation = FALSE, symbolic.cor = FALSE,
    ...)
{
    z <- object
    p <- z$rank
    rdf <- z$df.residual
    if (p == 0) {
        r <- z$residuals
        n <- length(r)
        w <- z$weights
        if (is.null(w)) {
            rss <- sum(r^2)
        }
        else {
            rss <- sum(w * r^2)
            r <- sqrt(w) * r
        }
        resvar <- rss/rdf
        ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
        class(ans) <- "summary.lm"
        ans$aliased <- is.na(coef(object))
        ans$residuals <- r
        ans$df <- c(0L, n, length(ans$aliased))
        ans$coefficients <- matrix(NA_real_, 0L, 4L, dimnames = list(NULL,
            c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
        ans$sigma <- sqrt(resvar)
        ans$r.squared <- ans$adj.r.squared <- 0
        ans$cov.unscaled <- matrix(NA_real_, 0L, 0L)
        if (correlation)
            ans$correlation <- ans$cov.unscaled
        return(ans)
    }
    if (is.null(z$terms))
        stop("invalid 'lm' object:  no 'terms' component")
    if (!inherits(object, "lm"))
        warning("calling summary.lm(<fake-lm-object>) ...")
    Qr <- qr.lmFAB(object)
    n <- NROW(Qr$qr)
    if (is.na(z$df.residual) || n - p != z$df.residual)
        warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
    r <- z$residuals
    f <- z$fitted.values
    w <- z$weights
    if (is.null(w)) {
        mss <- if (attr(z$terms, "intercept"))
            sum((f - mean(f))^2)
        else sum(f^2)
        rss <- sum(r^2)
    }
    else {
        mss <- if (attr(z$terms, "intercept")) {
            m <- sum(w * f/sum(w))
            sum(w * (f - m)^2)
        }
        else sum(w * f^2)
        rss <- sum(w * r^2)
        r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    if (is.finite(resvar) && resvar < (mean(f)^2 + var(c(f))) *
        1e-30)
        warning("essentially perfect fit: summary may be unreliable")
    p1 <- 1L:p
    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    se <- sqrt(diag(R) * resvar)
    est <- z$coefficients[Qr$pivot[p1]]
    tval <- est/se
    ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    ans$residuals <- r 

    ## pvals 
    pU<-2 * pt(abs(tval), rdf,lower.tail = FALSE)
    pF<-z$FABpv
    nU<-length(pU)-length(pF) 
    pF<-c(pU[seq(1,nU,length=nU)],pF) 

    ans$coefficients <- cbind(Estimate = est, `Std. Error` = se,
        `t value` = tval, `Pr(>|t+bfab|)` =pF)

    #ans$coefficients <- cbind(Estimate = est, `Std. Error` = se,
    #    `t value` = tval, `Pr(>|t|)` = 2 * pt(abs(tval), rdf,
    #        lower.tail = FALSE)) 

    ans$aliased <- is.na(z$coefficients)
    ans$sigma <- sqrt(resvar)
    ans$df <- c(p, rdf, NCOL(Qr$qr))
    if (p != attr(z$terms, "intercept")) {
        df.int <- if (attr(z$terms, "intercept"))
            1L
        else 0L
        ans$r.squared <- mss/(mss + rss)
        ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n -
            df.int)/rdf)
        ans$fstatistic <- c(value = (mss/(p - df.int))/resvar,
            numdf = p - df.int, dendf = rdf)
    }
    else ans$r.squared <- ans$adj.r.squared <- 0
    ans$cov.unscaled <- R
    dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1,
        1)]
    if (correlation) {
        ans$correlation <- (R * resvar)/outer(se, se)
        dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
        ans$symbolic.cor <- symbolic.cor
    }
    if (!is.null(z$na.action))
        ans$na.action <- z$na.action
    class(ans) <- c("summary.lmFAB","summary.lm")
    ans
}

#' @title QR decomposition 
#' 
#' @description QR decomposition for lmfFAB objects
#' 
#' @param x \code{lmFAB} object 
#' @param ... see \code{qr.lm}, if you can find it
#' 
#' @export
qr.lmFAB<-function (x, ...)
{
    if (is.null(r <- x$qr)) 
        stop("lm object does not have a proper 'qr' component.\n Rank zero or should not have used lm(.., qr=FALSE).")
    r
}

#' @title Summarizing Generalized Linear Model Fits with FAB Inference
#' 
#' @description \code{summary} method for class 'glmFAB'
#' 
#' @param object an object of class \code{glmFAB}
#' @param dispersion see \code{summary.glm}
#' @param correlation see \code{summary.glm} 
#' @param symbolic.cor see \code{summary.glm} 
#' @param ... see \code{summary.glm}
#'
#' @details A mod of \code{summary.glm} that shows FAB p-values in table
#'
#' @export 
summary.glmFAB<-function(object, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE,
    ...)
{
    est.disp <- FALSE
    df.r <- object$df.residual
    if (is.null(dispersion))
        dispersion <- if (object$family$family %in% c("poisson",
            "binomial"))
            1
        else if (df.r > 0) {
            est.disp <- TRUE
            if (any(object$weights == 0))
                warning("observations with zero weight not used for calculating dispersion")
            sum((object$weights * object$residuals^2)[object$weights >
                0])/df.r
        }
        else {
            est.disp <- TRUE
            NaN
        }
    aliased <- is.na(coef(object))
    p <- object$rank
    if (p > 0) {
        p1 <- 1L:p
        Qr <- qr.lmFAB(object)
        coef.p <- object$coefficients[Qr$pivot[p1]]
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
        covmat <- dispersion * covmat.unscaled
        var.cf <- diag(covmat)
        s.err <- sqrt(var.cf)
        tvalue <- coef.p/s.err
        dn <- c("Estimate", "Std. Error")
        if (!est.disp) {
            pvalue <- 2 * pnorm(-abs(tvalue))
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef.p), c(dn,
                "z value", "Pr(>|z+bfab|)"))
        }
        else if (df.r > 0) {
            pvalue <- 2 * pt(-abs(tvalue), df.r)
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef.p), c(dn,
                "t value", "Pr(>|z+bfab|)"))
        }
        else {
            coef.table <- cbind(coef.p, NaN, NaN, NaN)
            dimnames(coef.table) <- list(names(coef.p), c(dn,
                "t value", "Pr(>|z+bfab|)"))
        }
        df.f <- NCOL(Qr$qr)
    }
    else {
        coef.table <- matrix(, 0L, 4L)
        dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error",
            "t value", "Pr(>|z+bfab|)"))
        covmat.unscaled <- covmat <- matrix(, 0L, 0L)
        df.f <- length(aliased)
    }

    ## substitute in p-vals 
    pW<-coef.table[,4]
    pF<-object$FABpv  
    nW<-length(pW)-length(pF) 
    pF<-c(pW[seq(1,nW,length=nW)],pF)      
    coef.table[,4]<-pF 

    keep <- match(c("call", "terms", "family", "deviance", "aic",
        "contrasts", "df.residual", "null.deviance", "df.null",
        "iter", "na.action"), names(object), 0L)
    ans <- c(object[keep], list(deviance.resid = residuals(object,
        type = "deviance"), coefficients = coef.table, aliased = aliased,
        dispersion = dispersion, df = c(object$rank, df.r, df.f),
        cov.unscaled = covmat.unscaled, cov.scaled = covmat))
    if (correlation && p > 0) {
        dd <- sqrt(diag(covmat.unscaled))
        ans$correlation <- covmat.unscaled/outer(dd, dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- c("summary.glmFAB","summary.glm") 
    return(ans)
}

