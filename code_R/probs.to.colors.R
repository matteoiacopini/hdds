#' Converts a vector of probabilities into an array defining a set of colors.
#'
#' This function allows you to create a HDDS for comparing two distributions.
#' @param x          probabilty vector
#' @param ramp       matrix of colors (RGB or other encoding)
#' @param from       lower bound color scale, in relative scale (0,1)
#' @param to         upper bound color scale, in relative scale (0,1)
#' @param inv_scale  inverse scaling factor
#'
#' @keywords colors
#'
#' @examples
#' data = rnorm(800,1,2)
#' pdf  = dnorm(data,1,2)
#' colramp = colorRampPalette(c("white","black"))(100)
#' probs.to.colors = function(x=pdf, ramp=colramp, from=0, to=1, inv_scale=1.0)
#'
#' @return matrix referring to color shades for each point in the probability vector
#'
#' @import grDevices
#'
#' @export


probs.to.colors = function(x, ramp, from=0, to=1, inv_scale){
   ramp[findInterval(x, seq(from,to/inv_scale,length.out=length(ramp)+1), all.inside=T)]
}
