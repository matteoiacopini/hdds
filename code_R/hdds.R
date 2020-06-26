#' Creates and plots a single HDDS from data points
#'
#' This function allows you to create a HDDS for comparing two distributions.
#' @param data            data for distribution
#' @param discrete        data: discrete (1 bin for each point) <--> continuous (kernel estimation)
#' @param diameter        length of diameter for the disk
#' @param bounds          (vector) of lower and upper bounds where to truncate the distribution for visualisation
#' @param kernel.type     (string) type of kernel for density estimation
#' @param kernel.n        number of evaluation points in density estimation
#' @param add.median        add median line?
#' @param add.data.points   add data points?
#' @param color.data        color of data points
#' @param jitter.data       jitter data points? (otherwise similar data points are superimposed)
#' @param colors.from.to    (vector of 2 strings) bounds of color palette
#' @param color.n           number of colors in the palette
#' @param circle.center   center of the disk
#' @param circle.bins     number of bins in the disk
#' @param inv.scale       inverse rescaling factor of the densities
#' @param gamma           gamma correction: < 1 --> darken the tails;  > 1 shorten the black area around peak
#' @param legend          TRUE/FALSE, add legend ?
#' @param labels          TRUE/FALSE, add label for the bounds ?
#' @param labels.sz       labels size
#' @param ticks           ticks to be displayed as thin lines
#' @param tlen            ticks line length
#' @param tcol            ticks line color
#' @param twd             ticks line width
#'
#' @keywords HDDS
#'
#' @examples
#' # discrete, uni-modal
#' DATA1 = rnbinom(400,10,0.5)
#' colors.from.to = c(hcl.colors(70,"Reds")[70], hcl.colors(70,"Reds")[1]);
#' hdds(data=DATA1, discrete=TRUE, bounds=NA, add.median=TRUE,
#'      add.data.points=FALSE, colors.from.to=colors.from.to, labels=TRUE)
#'
#' # continuous, uni-modal
#' DATA2 = rgamma(200,4.0,0.5)
#' colors.from.to = c(hcl.colors(70,"Blues")[70], hcl.colors(70,"Blues")[1]);
#' color.data = "black"
#' hdds(data=DATA2, discrete=FALSE, bounds=NA, add.median=TRUE,
#'      colors.from.to=colors.from.to, inv.scale=20,
#'      add.data.points=TRUE, color.data=color.data, labels=TRUE)
#'
#' @return ggplot object (list of values) for plotting
#'
#' @references{
#'    \insertRef{paperHDDS}{hdds}
#' }
#'
#' @importFrom Rdpack reprompt
#'
#' @import stats
#' @import ggplot2
#' @import grDevices
#' @import graphics
#'
#' @export

hdds    = function(data,
                   discrete = FALSE,      # data discrete (1 bin for each point) <--> continuous (kernel estimation)
                   diameter = 1.00,       # diameter of the HDDS
                   bounds   = NA,         # bounds of the HDDS
                   kernel.type = c("epanechnikov"),
                   kernel.n    = 512,
                   add.median      = TRUE,        # add median line?
                   add.data.points = FALSE,       # add data points?
                   color.data      = "black",     # color of data points
                   jitter.data     = TRUE,        # jitter data points? (otherwise similar data points are superimposed)
                   colors.from.to  = c(hcl.colors(70,"Grays")[70], hcl.colors(70,"Grays")[1]), #c("white","black"),
                   color.n         = 100,
                   circle.center = c(0,0.02),
                   circle.bins   = 100,
                   inv.scale = 10,     # scaling factor for color intensity
                   gamma     = 1.0,    # gamma correction:  < 1 --> darken the tails;  > 1 shorten the black area around peak
                   legend    = FALSE,
                   labels    = FALSE,
                   labels.sz = 9.0,   # labels size
                   ticks = NULL,      # ticks
                   tlen  = NULL,      # ticks line length
                   tcol  = "black",   # ticks line color
                   twd   = 1.00       # ticks line width
){

   # preliminary checks
   if (!is.numeric(data)){
      stop("\'data\' must be numeric")
   }
   if (any(is.na(data))){
      warning("Data contains NA")
   }
   if (diameter <= 0){
      stop("\'diameter\' must be a positive scalar")
   }

   # determine bounds of the HDDS
   if (any(is.na(bounds))){
      bounds = c(min(data), max(data))
   }else if (is.double(bounds)){
      bounds = sort(bounds, decreasing = FALSE)
   }else if (is.character(bounds)){
      if(bounds == "quantile_real"){
         bounds = quantile(data,c(0.05,0.95))
      }else if (bounds == "quantile_pos"){
         bounds = quantile(data,c(0.0,0.90))
      }else{
         stop("Bounds not correct. Specify a 2-dim vector in ascending order, or \"quantile_real\", or \"quantile_pos\"")
      }
   }

   # determine color palette
   ColorRamp = colorRampPalette(colors.from.to)(color.n)

   # estimate the DENSITY
   if(discrete){
      BINS   = seq(bounds[1], bounds[2], by=1)    # equispaced by 1
      breaks = seq(min(min(data),bounds[1]), max(max(data),bounds[2]), by=1)
      ll     = hist(data, breaks = breaks, plot = FALSE)$density
      colVal = ll[(breaks>=bounds[1]-1) & (breaks<=bounds[2]-1)]
      colVal[is.na(colVal)] = 0
      NBINS  = length(colVal)
   }else{
      NBINS = circle.bins
      BINS  = seq(bounds[1], bounds[2], length.out= NBINS) # by=diff(range(bounds_up))/(NBINS-1))
      if(kernel.type == "empirical"){
         colVal = hist(data, breaks= seq(bounds[1],bounds[2],length.out= color.n), plot=F)$density
      }else{
         DENS = density(data, adjust=1, kernel=kernel.type, n=kernel.n, from=bounds[1], to=bounds[2])
         DF   = splinefun(DENS) #approxfun(DENS)
         colVal = DF(BINS)
         colVal[is.na(colVal)] = 0
      }
   }
   # TRANSFORM: density --> color shading
   dens = (colVal/sum(colVal))^gamma
   BIN.COLORS = probs.to.colors(x=dens, ramp=ColorRamp, from=0, to=1, inv_scale=inv.scale)



   # create HDDS object
   colBin = circleFun_up(center=circle.center, diameter=diameter, start=0, end=1/NBINS, filled=F)
   hc = ggplot() +
      geom_polygon(data=colBin, aes(x,y,colour=colBin[,1]), fill=BIN.COLORS[NBINS]) +
      coord_equal() +
      xlim(-0.58, 0.58) +
      ylim(-0.05, 0.56) +
      theme_bw()

   # fill-in the HDDS
   for(i in 1:(NBINS-0)){
      colBin = circleFun_up(center=circle.center, diameter=diameter, start= 1-i/NBINS, end= 1-(i-1)/NBINS, filled=T)
      hc = hc + geom_polygon(data=colBin, aes(x,y), fill=BIN.COLORS[i])
   }

   # circle border
   colBin = circleFun_up(center=circle.center, diameter=diameter, start=0, end=1, filled=F)
   hc = hc + geom_path(data=colBin, aes(x,y), color=colors.from.to[2], size=0.5)


   # add MEDIAN line
   if (add.median){
      med = median(data)
      # compute coordinates of median along the half circle
      # |med-A| / |B-A|  =  arc(a,med) / arc(a,b)     where arc(a,b)= pi*r  (semicirle length)
      rr = diameter/2
      arc = (med-bounds[1])*pi*rr / (bounds[2]-bounds[1])
      angle = arc/rr
      # NOTE: angle is counterclockwise oriented, from the point (1,0). To obtain the opposite orientation:
      #       (cos(pi-angle), sin(pi-angle)) * rr + circle.center, or equivalently:
      coord_med = c( -cos(angle), sin(angle) ) * rr + circle.center
      # compute starting and ending points of the line passing through the median
      v = coord_med-circle.center
      u = v / sqrt(sum(v^2))
      step = (0.05*sqrt(sum(v^2)))
      coord_start = coord_med - step*u
      coord_end   = coord_med + step*u
      med_data = data.frame(
         x = coord_start[1],
         xend = coord_end[1],
         y = coord_start[2],
         yend = coord_end[2]
      )
      hc = hc + geom_segment(data= med_data, aes(x= x, y= y, xend= xend, yend= yend),
                             color="black", linetype="solid", size=1.10)
   }


   # add TICKS
   if (!is.null(ticks)){
      # preliminary check
      if (any(ticks < bounds[1] | ticks > bounds[2])){
         warning("Some ticks are outside bounds and will not be displayed.")
      }
      # compute coordinates of ticks along the half circle
      cc  = t(replicate(length(ticks), circle.center))
      rr  = diameter/2
      arc = (ticks-bounds[1])*pi*rr / (bounds[2]-bounds[1])
      angle = arc/rr
      # NOTE: angle is counterclockwise oriented, from the point (1,0). To obtain the opposite orientation:
      #       (cos(pi-angle), sin(pi-angle)) * rr + circle.center, or equivalently:
      coord_ticks = cbind( -cos(angle), sin(angle) ) * rr + cc
      # compute starting and ending points of the line passing through the ticks
      v = coord_ticks - cc
      if(length(ticks)==1){
         u = v / sqrt(sum(v^2))
      }else{
         u = diag( 1/sqrt(rowSums(v^2)) ) %*% v  # u = v / sqrt(sum(v^2))
      }
      if(is.null(tlen)){
         tlen = (0.05*sqrt(rowSums(v^2)))
         tlen = replicate(2, tlen)
      }else if(is.atomic(tlen) && length(tlen) == 1L){
         tlen = tlen * matrix(1,length(ticks),2)
      }
      coord_start = coord_ticks - tlen*u
      coord_end   = coord_ticks + tlen*u

      segment_data = data.frame(
         x = coord_start[,1],
         xend = coord_end[,1],
         y = coord_start[,2],
         yend = coord_end[,2]
      )
      hc = hc + geom_segment(data= segment_data, aes(x= x, y= y, xend= xend, yend= yend),
                             color=tcol, linetype="solid", size=twd)
   }


   # add DATA POINTS
   if(add.data.points){
      dots = data[data > bounds[1] & data < bounds[2]]
      # apply map from [bounds[1],bounds[2]] -> [-1,1]
      r_points = diameter/2 - 0.04
      # equispaced points on the semi-circle: same arc length
      xx = vector(mode = "numeric", NBINS+1)
      yy = vector(mode = "numeric", NBINS+1)
      stp = (pi * r_points) / NBINS       # arc length. NOTE: (pi * r) = length of semi-circle with radius r
      for(i in 1:(NBINS+1)){
         theta = (i-1)*stp / r_points     # get angle under arc
         xx[i] = r_points * cos(theta)    # get x-coordinate
         yy[i] = r_points * sin(theta)    # get y-coordinate
      }
      xx = sort(xx)
      yy = yy + circle.center[2]
      # Map dots into (xx,yy) coordinates, such that: min(dots)-> xx[NBINS], max(dots)-> xx[1]
      xend_pts = vector(mode = "numeric", length(dots))
      yend_pts = vector(mode = "numeric", length(dots))
      for(j in 1:length(dots)){
         ll = hist(dots[j], breaks = c(bounds[1]-1,seq(bounds[1], bounds[2], length.out = NBINS)), plot = FALSE)$counts
         xend_pts[j] = xx[as.logical(ll)]
         yend_pts[j] = yy[as.logical(ll)]
      }
      if(jitter.data){
         xend_pts = xend_pts + rnorm(length(xend_pts),0,0.008)
         yend_pts = yend_pts + rnorm(length(yend_pts),0,0.008)
      }
      hc = hc + geom_point(aes(x = xend_pts, y = yend_pts), colour=color.data, size = 0.80)
   }


   # add LABELS
   if(labels == TRUE){
      hc = hc +
         annotate("text", x=circle.center[1]+0.56*diameter, y=circle.center[2]+0.02*diameter,
                  label=round(max(BINS),2), size=labels.sz, color="black") +
         annotate("text", x=circle.center[1]-0.56*diameter, y=circle.center[2]+0.02*diameter,
                  label=round(min(BINS),2), size=labels.sz, color="black")
   }


   # add LEGEND
   if(legend == TRUE){
      hc = hc +
         theme(axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
               axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
               panel.background=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(), plot.background=element_blank()) +
         scale_colour_gradient(low= colors.from.to[1], high= colors.from.to[2], limits=c(0,1), name="Density") +
         theme(legend.direction = "horizontal", legend.position=c(0.47,-0.01)) #c(0.47,0.05) #"bottom") #c(0.98,0.5))
   }else{
      hc = hc +
         theme(axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
               axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
               panel.background=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(), plot.background=element_blank(), legend.position="none")
   }

   return(hc)
}
