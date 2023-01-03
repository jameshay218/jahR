#' Convert years to text quarters
#' 
#' Converts a vector of times, as years with decimals, into quarters. Appends "QX-" to the start of the number, where X is the quarter number.
#' @param times a vector of years
#' @return a character vector of modified times
#' @export
convert_years_to_quarters <- function(times){
    times <- sapply(time_key, function(x){
        if(x %% 1 == 0){
            paste0("Q1-",x)
        } else if(x %% 1 == 0.25){
            paste0("Q2-",x-0.25)
        } else if(x%%1 == 0.5){
            paste0("Q3-",x-0.5)
        } else {
            paste0("Q4-",x-0.75)
        }
    })
    times
}

#' Get days per month of the year
#'
#' Returns a vector of days per month
#' @param breaks total number of months to consider (defaults to 12 to give raw days per each month)
#' @return vector of days
#' @export
get_days_per_month <- function(breaks=12){
    days <- c(31,28,31,30,31,30,31,31,30,31,30,31)
    days <- colSums(matrix(days,ncol=breaks))
    return(days)
}

#' Convert integer to date
#' 
#' @param x integer to be converted
#' @param origin character of format "1970-01-01" as date origin
#' @return a date
#' @export
convert_number_to_date <- function(x, origin="2013-01-01"){
    y <- as.Date(x,origin=origin)
    return(y)
}

#' Convert date to integer
#' 
#' @param x date to be converted
#' @param origin character of format "1970-01-01" as date origin
#' @return an integer
#' @export
convert_date_to_number <- function(x,origin="2013-01-01"){
    y <- as.numeric(as.Date(x,origin=origin)) - as.numeric(as.Date(origin,origin=origin))
    return(y)
}