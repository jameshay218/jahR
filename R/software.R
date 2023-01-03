#' Protect function
#'
#' Wrapper function to protect calls to the posterior function. If posterior does not compute correctly, returns -100000.
#' @param f the function to be protected
#' @return the protected function
#' @export
protect <- function(f) {
    function(...) {
        tryCatch(f(...), error = function(e) {
            message("caught error: ", e$message)
            -10000000
        })
    }
}

