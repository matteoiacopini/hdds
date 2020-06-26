#' Creates a sequence of PNG files and animated plot of an HDDS, from data points
#'
#' This function allows you to create a HDDS for comparing two distributions.
#' @param data            data for distribution
#' @param gif.create      do you want to create a GIF file in output?
#' @param gif.name        name of the GIF file
#' @param ...             additional parrameters to be passed to the function hdds
#'
#' @keywords animation HDDS
#'
#' @examples
#' DATA = rgamma(15,4.0,0.5)
#' colors.from.to = c(hcl.colors(70,"Blues")[70], hcl.colors(70,"Blues")[1])
#' color.data = "black"
#' animated_hdds(data=DATA, gif.create=TRUE, gif.name="animated_HDDS_gamma",
#'               discrete=FALSE, bounds=NA, add.median=TRUE,
#'               colors.from.to=colors.from.to,
#'               add.data.points=TRUE, color.data=color.data)
#'
#' @return a set of PNG files containing the HDDSs as each observation comes in, and an animated HDDS as GIF file (optional)
#'
#' @import stats
#' @import ggplot2
#' @import grDevices
#' @import graphics
#' @import gridExtra
#' @import magick
#' @import purrr
#'
#' @export

animated_hdds = function(data, gif.create=TRUE, gif.name="animated_HDDS", ...){

   Ntot = length(data)  # total number of data points
   N0 = 3               # number of initial data points used for first kernel estimate
   N  = Ntot - N0

   # print frames: one HDDS by adding one observation at time
   for(idx in 1:N){
      # hdds(data[1:(idx+N0)], discrete=F, bounds=bounds, gamma=0.8, colors.from.to=col,
      #      add.data.points=T, jitter.data=F, labels=T)
      hdds(data[1:(idx+N0)], ...)
      ggsave(paste0("frame.",idx,".png"), width=15, height=6, units="cm")

      # convert each png file as magick object
      png.dat = image_read(paste0("frame.",idx,".png"))
      assign(paste0("frame.",idx),png.dat)
   }

   animation = NULL
   if(gif.create){
      # stack the frames
      df = expand.grid(idx = 1:N)
      df$file = paste0("frame.",df$idx,".png")
      images  = map(df$file, image_read)
      images  = image_join(images)
      # create animated GIF
      animation = image_animate(images, fps = 1)
      image_write(animation, paste0(gif.name,".gif"))
   }
   return(animation)
}
