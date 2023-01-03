## From prodlim package - finds matching rows between two data frames. "Thus the function returns a vector with the row numbers of (first) matches of its first argument in its second.", https://www.rdocumentation.org/packages/prodlim/versions/2018.04.18/topics/row.match
#' @export
row_match <- function(x, table, nomatch = NA) {
    if (class(table) == "matrix") {
        table <- as.data.frame(table)
    }
    if (is.null(dim(x))) {
        x <- as.data.frame(matrix(x, nrow = 1))
    }
    cx <- do.call("paste", c(x[, , drop = FALSE], sep = "\r"))
    ct <- do.call("paste", c(table[, , drop = FALSE], sep = "\r"))
    match(cx, ct, nomatch = nomatch)
}

#' Sums nth elements
#'
#' Sums every n values of a numeric vector
#' @param v the vector to be summed
#' @param n the number of elements to sum
#' @export
sum_n_vector <- function(v,n){
    nv <- length(v)
    if (nv %% n)
        v[ceiling(nv / n) * n] <- NA
    colSums(matrix(v, n), na.rm=TRUE)
}

#' Prints the average and quantiles of a vector
#'
#' @param x the vector of values to summarize
#' @param average character describing the average version to use (mean or median)
#' @param quantiles vector of length 2 giving the quantiles to calculate
#' @param nsignif number of significant figures to report
#' @return a string with the formatted estimate
#' @export
print_estimate <- function(x,average="mean", quantiles=c(0.025,0.975),nsignif=3){
    if(average=="mean"){
        par1 <- mean(x,na.rm=TRUE)
    } else if(average=="median"){
        par1 <- median(x, na.rm=TRUE)
    } else {
        message("Warning - average should be mean or median")
        par1 <- mean(x,na.rm=TRUE)
    }
    
    if(length(quantiles) > 2){
        message("Warning - quantiles should be of length 2")
    }
    
    par2 <- signif(quantile(x, quantiles,na.rm=TRUE),nsignif)
    label <- paste0(signif(par1,nsignif), "; (", par2[1], "-",par2[2],")")
    return(label)
}