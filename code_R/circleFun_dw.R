#' Creates a circular sector, starting from the III quadrant and ending in the IV quadrant,
#' moving in counterclockwise direction.
#'
#' This function allows you to create a HDDS for comparing two distributions.
#' @param center   cener of the disk
#' @param diameter length of diameter of the disk
#' @param npoints  number of discretization points in the disk
#' @param start    lower bound
#' @param end      upper bound
#' @param filled TRUE/FALSE, fill-in the disk?
#'
#' @keywords HDDS
#'
#' @examples
#' circleFun_dw(diameter = 1.0)
#'
#' @return matrix object with values relating to the color shading
#'
#' @export

circleFun_dw = function(center=c(0,-0.02), diameter, npoints=100, start=0, end=2, filled=TRUE){
   # REVERSE order, to making correct density strip (values from left to right)
   tt = seq(end*pi, start*pi, length.out=npoints)  # tt = seq(start*pi, end*pi, length.out=npoints)
   output = data.frame(
      x = center[1] + diameter / 2 * cos(tt),
      y = center[2] - diameter / 2 * sin(tt)
   )
   if(filled==T){
      output = rbind(output, center)
   }else{
      output = rbind(output, output[npoints,])
   }
   return(output)
}
