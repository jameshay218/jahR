#' My overall ggplot2 theme
#'
#' Adds my custom ggplot2 theme to a ggplot2 object, with the intention of standardizing my plots.
#' @keywords ggplot2, theme
#' @return A ggplot2 theme object
#' @export
#' @import ggplot2
#' @examples
#' ggplot(data=data.frame(x=rnorm(100))) + geom_histogram(aes(x=x)) + theme_overall()
theme_overall <- function(){
  theme_minimal() + theme(axis.text=element_text(size=7),
                     axis.title = element_text(size=8),
                     plot.title=element_text(size=10),
                     plot.tag = element_text(size=10,face="bold"),
                     legend.text=element_text(size=7),legend.title=element_text(size=8))
}

#' ggplot2 theme no x axis
#'
#' Removes the x-axis line, ticks and text from a ggplot2 object
#' @keywords ggplot2, theme
#' @return A ggplot2 theme object
#' @export
#' @examples
#' ggplot(data=data.frame(x=rnorm(100))) + geom_histogram(aes(x=x)) + theme_no_x_axis()
theme_no_x_axis <- function(){
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank(),axis.line.x=element_blank())
}

#' ggplot2 theme nice axes
#'
#' Makes some nice axis lines for a ggplot2 object
#' @keywords ggplot2, theme
#' @return A ggplot2 theme object
#' @export
#' @examples
#' ggplot(data=data.frame(x=rnorm(100))) + geom_histogram(aes(x=x)) + theme_nice_axes()
theme_nice_axes <- function(){
  theme_nice_axes <- theme(axis.line.x=element_line(size=0.5,color="black"),axis.line.y=element_line(size=0.5,color="black"),panel.border = element_blank())
}


