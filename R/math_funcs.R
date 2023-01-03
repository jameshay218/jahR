#' Logit function
#' 
#' Converts a given value from the scale of 0 and 1 to between -Inf and Inf
#' @param p numeric value to be converted
#' @return numeric value of the converted input
#' @seealso \code{\link{logistic}} returns the inverse value (ie. the logistic function)
#' @export
#' @examples
#' logit(0.5)
logit <- function(p){
    return(log(p/(1-p)))
}

#' Logistic function
#' 
#' Converts a given value from the scale of -Inf and +Inf to 0 and 1
#' @param p numeric value to be converted
#' @return numeric value of the converted input
#' @seealso \code{\link{logistic}} returns the inverse value (ie. the logit function)
#' @export
#' @examples
#' logistic(56)
#' logistic(logit(0.5))
logistic <- function(p){
    return(1/(1+exp(-p)))
}

#' Custom logistic transformation
#' 
#' Transforms a given value from between -Inf and Inf to xmin and xmax
#' @param x numeric value to be scaled
#' @param xmin numeric value representing the lower bound of the new scale
#' @param xmax numeric value representing the upper bound of the new scale
#' @return numeric value with the converted input
#' @seealso \code{\link{transform_logit}} inverse function
#' @export
#' @examples
#' transform_logistic(5,3,10)
transform_logistic <- function(x, xmin, xmax){
    y <- (xmax-xmin)/(1 + exp(-x)) + xmin
    return(y)
}
#' Custom logit transform
#' 
#' Transforms a given value from between to xmin and xmax to -Inf and Inf
#' @param x numeric value to be scaled
#' @param xmin numeric value representing the lower bound of the old scale
#' @param xmax numeric value representing the upper bound of the old scale
#' @return numeric value with the converted input
#' @seealso \code{\link{transform_logistic}} inverse function
#' @export
#' @examples
#' transform_logit(0,3,10)
transform_logit <- function(x, xmin, xmax){
    y <- -log(((xmax-xmin)/(x-xmin))-1)
    return(y)
}
#' to_unit_scale
#' 
#' Transforms a given value from the scale of min-max to 0-1, either directly from a linear scale or via a log scale
#' @param y numeric value to be transformed
#' @param min the bottom of the linear scale. Defaults to 1
#' @param max the top of the linear scale. Defaults to 100
#' @param logflag Optional flag to convert to a unit scale via the log scale. Defaults to FALSE
#' @param logbase Optional numeric value to be used as the base for the log scale
#' @return numeric value with the converted input
#' @seealso \code{\link{fromUnitScale}} inverse function
#' @export
#' @examples
#' to_unit_scale(65,20,80)
#' to_unit_scale(15,0,16,TRUE,10)
to_unit_scale <- function(y, min=1,max=100,logflag=FALSE,logbase=10){
    if(logflag){
        rtn <- (log(y,logbase)-log(min,logbase))/(log(max,logbase)-log(min,logbase))
    }
    else{
        rtn <- (y-min)/(max-min)
    }
    rtn
}
#' from_unit_scale
#' 
#' Transforms a given value from the scale of 0-1 to min-max, either directly to a linear scale or via a log scale
#' @param x numeric value to be transformed
#' @param min the bottom of the linear scale. Defaults to 1
#' @param max the top of the linear scale. Defaults to 100
#' @param logflag Optional flag to convert to a unit scale via the log scale. Defaults to FALSE
#' @param logbase Optional numeric value to be used as the base for the log scale
#' @return numeric value with the converted input
#' @seealso \code{\link{toUnitScale}} inverse function
#' @export
#' @examples
#' from_unit_scale(0.5,20,80)
#' from_unit_scale(0.33,0,16,TRUE,10)
from_unit_scale <- function(x,min=1,max=100,logflag=FALSE,logbase=10){
    if(logflag){
        rtn <- min*logbase^(x*(log(max,logbase)-log(min,logbase)))
    } 
    else{
        rtn <- min + (max-min)*x
    }
    rtn
}

