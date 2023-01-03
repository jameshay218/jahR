#' Saves a pdf and png of a ggplot2 object
#' 
#' @param p the ggplot2 object to be saved
#' @param wd the full working directory to save to
#' @param save_name the filename, without file extension
#' @param width the width of the plot in inches
#' @param height the height of the plot in inches
#' @return NULL
#' @export
save_plots <- function(p, wd, save_name, width, height){
    message("Saving to ", paste0(wd, save_name,".pdf/png"))
    ggsave(p, filename=paste0(wd, save_name,".pdf"),width=width,height=height)
    ggsave(p, filename=paste0(wd, save_name,".png"),width=width,height=height,units='in',dpi=300)
    message("Success")
    NULL
}