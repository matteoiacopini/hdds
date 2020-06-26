#' Creates and plots a two HDDS from two sets of data points, vertically stacking them
#' by turning upside down the second plot, to facilitate comparisons.
#'
#' This function allows you to create a HDDS for comparing two distributions.
#' @param data.up         data for distribution on top
#' @param data.dw         data for distribution on bottom
#' @param discrete.up     data on top: discrete (1 bin for each point) <--> continuous (kernel estimation)
#' @param discrete.dw     data on bottom: discrete (1 bin for each point) <--> continuous (kernel estimation)
#' @param diameter.up     length of diameter for the disk on top
#' @param diameter.dw     length of diameter for the disk on bottom
#' @param bounds.up       (vector) of lower and upper bounds where to truncate the distribution on top for visualisation
#' @param bounds.dw       (vector) of lower and upper bounds where to truncate the distribution on bottom for visualisation
#' @param kernel.type.up  (string) type of kernel for density estimation on top
#' @param kernel.type.dw  (string) type of kernel for density estimation on bottom
#' @param kernel.n        number of evaluation points in density estimation
#' @param add.median.up       add median line on top?
#' @param add.median.dw       add median line on bottom?
#' @param add.data.points.up  add data points on top?
#' @param add.data.points.dw  add data points on bottom?
#' @param color.data.up       color of data points on top
#' @param color.data.dw       color of data points on bottom
#' @param jitter.data         jitter data points? (otherwise similar data points are superimposed)
#' @param colors.from.to.up   (vector of 2 strings) bounds of color palette for top disk
#' @param colors.from.to.dw   (vector of 2 strings) bounds of color palette for bottom disk
#' @param color.n           number of colors in the palette
#' @param circle.center.up  center of the disk on top
#' @param circle.center.dw  center of the disk on bottom
#' @param circle.bins     number of bins in the disk
#' @param inv.scale       inverse rescaling factor of the densities
#' @param gamma.up        gamma correction top: < 1 --> darken the tails;  > 1 shorten the black area around peak
#' @param gamma.dw        gamma correction bottom: < 1 --> darken the tails;  > 1 shorten the black area around peak
#' @param legend          TRUE/FALSE, add legend ?
#' @param labels          TRUE/FALSE, add label for the bounds ?
#' @param labels.sz       labels size
#'
#' @keywords two HDDS, comparison
#'
#' @examples
#' # continuous, multi-modal
#' DATA3 = c(rnorm(floor(0.4*800),-4,1.1), rnorm(floor((1-0.4)*800),4,2))
#' # continuous, uni-modal
#' DATA4 = rnorm(800,0,3.0)
#' colors.from.to.up = c(hcl.colors(70,"Reds")[70], hcl.colors(70,"Reds")[1]);
#' colors.from.to.dw = c(hcl.colors(70,"Blues")[70], hcl.colors(70,"Blues")[1]);
#' hdds_updw(data.up = DATA3, data.dw = DATA4, discrete.up = FALSE, discrete.dw = FALSE,
#'           bounds.up = c(-6,6), bounds.dw = c(-6,6), add.median.up = TRUE, add.median.dw = TRUE,
#'           add.data.points.up = FALSE, add.data.points.dw = FALSE, labels = TRUE, inv.scale=20,
#'           colors.from.to.up = colors.from.to.up, colors.from.to.dw = colors.from.to.dw)
#'
#' @return ggplot object (list of values) for plotting
#'
#' @import stats
#' @import ggplot2
#' @import grDevices
#' @import graphics
#'
#' @export

hdds_updw = function(data.up,
                     data.dw,
                     discrete.up = FALSE,     # data discrete (1 bin for each point) <--> continuous (kernel estimation)
                     discrete.dw = FALSE,
                     diameter.up = 1.0,       # diameter of the HDDS
                     diameter.dw = 1.0,
                     bounds.up = NA,          # bounds of the HDDS
                     bounds.dw = NA,
                     kernel.type.up = c("epanechnikov"),
                     kernel.type.dw = c("epanechnikov"),
                     kernel.n  = 512,
                     add.median.up = TRUE,           # add median line?
                     add.median.dw = TRUE,
                     add.data.points.up = FALSE,     # add data points line?
                     add.data.points.dw = FALSE,
                     color.data.up      = "black",   # color of data points
                     color.data.dw      = "black",
                     jitter.data = TRUE,             # jitter data points? (otherwise similar data points are superimposed)
                     colors.from.to.up = c(hcl.colors(70,"Grays")[70], hcl.colors(70,"Grays")[1]), #c("white","black"),
                     colors.from.to.dw = c(hcl.colors(70,"Grays")[70], hcl.colors(70,"Grays")[1]), #c("white","black"),
                     color.n = 100,
                     circle.center.up = c(0,0.02),
                     circle.center.dw = c(0,-0.0), #c(0,-0.02),
                     circle.bins = 100,
                     inv.scale = 10,     # scaling factor for color intensity
                     gamma.up  = 1.0,    # gamma correction:  < 1 --> darken the tails;  > 1 shorten the black area around peak
                     gamma.dw  = 1.0,    # gamma correction:  < 1 --> darken the tails;  > 1 shorten the black area around peak
                     legend    = FALSE,
                     labels    = FALSE,
                     labels.sz = 9.0     # labels size
                     ){

   # preliminary checks
   if (any(is.na(data.up))){
      warning("Data.up contains NA")
   }
   if (any(is.na(data.dw))){
      warning("Data.dw contains NA")
   }
   if (diameter.up <= 0){
      stop("Diameter.up must be a positive scalar")
   }
   if (diameter.dw <= 0){
      stop("Diameter.dw must be a positive scalar")
   }

   ColorRamp.up = colorRampPalette(colors.from.to.up)(color.n)
   ColorRamp.dw = colorRampPalette(colors.from.to.dw)(color.n)


   # determine BOUNDS of the plot
   if (any(is.na(bounds.up))){
      bounds.up = c(min(data.up), max(data.up))
   }else if (is.double(bounds.up)){
      bounds.up = sort(bounds.up, decreasing = FALSE)
   }else if (is.character(bounds.up)){
      if(bounds.up == "quantile_real"){
         bounds.up = quantile(data.up,c(0.05,0.95))
      }else if (bounds.up == "quantile_pos"){
         bounds.up = quantile(data.up,c(0.0,0.90))
      }else{
         stop("Bounds.up not correct. Specify a 2-dim vector in ascending order, or \"quantile_real\", or \"quantile_pos\"")
      }
   }
   if (any(is.na(bounds.dw))){
      bounds.dw = c(min(data.dw), max(data.dw))
   }else if (is.double(bounds.dw)){
      bounds.dw = sort(bounds.dw, decreasing = FALSE)
   }else if (is.character(bounds.dw)){
      if(bounds.dw == "quantile_real"){
         bounds.dw = quantile(data.dw,c(0.05,0.95))
      }else if (bounds.dw == "quantile_pos"){
         bounds.dw = quantile(data.dw,c(0.0,0.90))
      }else{
         stop("Bounds.dw not correct. Specify a 2-dim vector in ascending order, or \"quantile_real\", or \"quantile_pos\"")
      }
   }



   # estimate the DENSITY
   if(discrete.up){
      BINS.up   = seq(bounds.up[1], bounds.up[2], by=1)    # equispaced by 1
      breaks.up = seq(min(min(data.up),bounds.up[1]), max(max(data.up),bounds.up[2]), by=1)
      ll.up     = hist(data.up, breaks = breaks.up, plot = FALSE)$density
      colVal.up = ll.up[(breaks.up>=bounds.up[1]-1) & (breaks.up<=bounds.up[2]-1)]
      colVal.up[is.na(colVal.up)] = 0
      NBINS.up  = length(colVal.up)
   }else{
      NBINS.up = circle.bins
      BINS.up  = seq(bounds.up[1], bounds.up[2], length.out= NBINS.up) # by=diff(range(bounds.up))/(NBINS-1))
      if(kernel.type.up == "empirical"){
         colVal.up = hist(data.up, breaks= seq(bounds.up[1],bounds.up[2], length.out= color.n), plot=F)$density
      }else{
         DENS.up = density(data.up, adjust=1, kernel=kernel.type.up, n=kernel.n, from=bounds.up[1], to=bounds.up[2])
         DF.up   = splinefun(DENS.up) #approxfun(DENS)
         colVal.up = DF.up(BINS.up)
         colVal.up[is.na(colVal.up)] = 0
      }
   }
   # TRANSFORM: kernel density --> density strip
   dens = (colVal.up/sum(colVal.up))^gamma.up
   BIN.COLORS.up = probs.to.colors(x=dens, ramp=ColorRamp.up, from=0, to=1, inv_scale=inv.scale)


   if(discrete.dw){
      BINS.dw   = seq(bounds.dw[1], bounds.dw[2], by=1)    # equispaced by 1
      breaks.dw = seq(min(min(data.dw),bounds.dw[1]), max(max(data.dw),bounds.dw[2]), by=1)
      ll.dw     = hist(data.dw, breaks.dw = breaks.dw, plot = FALSE)$density
      colVal.dw = ll.dw[(breaks.dw>=bounds.dw[1]-1) & (breaks.dw<=bounds.dw[2]-1)]
      colVal.dw[is.na(colVal.dw)] = 0
      NBINS.dw  = length(colVal.dw)
   }else{
      NBINS.dw = circle.bins
      BINS.dw  = seq(bounds.dw[1], bounds.dw[2], length.out= NBINS.dw) # by=diff(range(bounds.dw))/(NBINS-1))
      if(kernel.type.dw == "empirical"){
         colVal.dw = hist(data.dw, breaks= seq(bounds.dw[1],bounds.dw[2], length.out= color.n), plot=F)$density
      }else{
         DENS.dw = density(data.dw, adjust=1, kernel=kernel.type.dw, n=kernel.n, from=bounds.dw[1], to=bounds.dw[2])
         DF.dw   = splinefun(DENS.dw) #approxfun(DENS)
         colVal.dw = DF.dw(BINS.dw)
         colVal.dw[is.na(colVal.dw)] = 0
      }
   }
   # TRANSFORM: kernel density --> density strip
   dens = (colVal.dw/sum(colVal.dw))^gamma.dw
   BIN.COLORS.dw = probs.to.colors(x=dens, ramp=ColorRamp.dw, from=0, to=1, inv_scale=inv.scale)


   # create HDDS object
   colBin_up = circleFun_up(circle.center.up, diameter=diameter.up, start=0, end=1/NBINS.up, filled=F)
   colBin_dw = circleFun_dw(circle.center.dw, diameter=diameter.dw, start=0, end=1/NBINS.dw, filled=F)
   hc = ggplot() +
      geom_polygon(data=colBin_up, aes(x,y,colour=colBin_up[,1]), fill=BIN.COLORS.up[NBINS.up]) +
      geom_polygon(data=colBin_dw, aes(x,y,colour=colBin_dw[,1]), fill=BIN.COLORS.dw[NBINS.dw]) +
      coord_equal() +
      xlim(-0.58, 0.58) +
      ylim(-0.58, 0.58) +
      theme_bw()

   # fill-in the HDDS
   for(i in 1:(NBINS.up-0)){
      colBin_up = circleFun_up(circle.center.up, diameter=diameter.up, start= 1-i/NBINS.up, end= 1-(i-1)/NBINS.up, filled=T)
      hc = hc + geom_polygon(data=colBin_up, aes(x,y), fill=BIN.COLORS.up[i])
   }
   for(i in 1:(NBINS.dw-0)){
      colBin_dw = circleFun_dw(circle.center.dw, diameter=diameter.dw, start= 1-i/NBINS.dw, end= 1-(i-1)/NBINS.dw, filled=T)
      hc = hc + geom_polygon(data=colBin_dw, aes(x,y), fill=BIN.COLORS.dw[i])
   }

   # circle border
   colBin_up = circleFun_up(circle.center.up, diameter=diameter.up, start=0, end=1, filled=F)
   colBin_dw = circleFun_dw(circle.center.dw, diameter=diameter.dw, start=0, end=1, filled=F)
   hc = hc +
      geom_path(data=colBin_up, aes(x,y), color=colors.from.to.up[2], size=0.5) +
      geom_path(data=colBin_dw, aes(x,y), color=colors.from.to.dw[2], size=0.5)


   # add MEDIAN line
   if (add.median.up){
      med_up = median(data.up)
      # compute coordinates of median along the half circle
      # |med-A| / |B-A|  =  arc(a,med) / arc(a,b)     where arc(a,b)= pi*r  (semicirle length)
      rr_up    = diameter.up/2 #+0.015
      arc_up   = (med_up-bounds.up[1])*pi*rr_up / (bounds.up[2]-bounds.up[1])
      angle_up = arc_up/rr_up
      # NOTE: angle is counterclockwise oriented, from the point (1,0). To obtain the opposite orientation:
      #       (cos(pi-angle), sin(pi-angle)) * rr + circle.center, or equivalently:
      coord_med_up = c( -cos(angle_up), sin(angle_up) ) * rr_up + circle.center.up
      # compute starting and ending points of the line passing through the median
      v_up = coord_med_up - circle.center.up
      u_up = v_up / sqrt(sum(v_up^2))
      step_up = (0.05*sqrt(sum(v_up^2)))
      coord_start_up = coord_med_up - step_up*u_up
      coord_end_up   = coord_med_up + step_up*u_up
      med_data_up = data.frame(
         x = coord_start_up[1],
         xend = coord_end_up[1],
         y = coord_start_up[2],
         yend = coord_end_up[2]
      )
      hc = hc + geom_segment(data= med_data_up, aes(x= x, y= y, xend= xend, yend= yend),
                             color="black", linetype="solid", size=1.10)
   }
   if (add.median.dw){
      med_dw = median(data.dw)
      # compute coordinates of median along the half circle
      rr_dw    = diameter.dw/2 #+0.015
      arc_dw   = (med_dw-bounds.dw[1])*pi*rr_dw / (bounds.dw[2]-bounds.dw[1])
      angle_dw = arc_dw/rr_dw
      # NOTE: angle is counterclockwise oriented, from the point (1,0). To obtain the opposite orientation:
      #       (cos(pi-angle), sin(pi-angle)) * rr + circle.center, or equivalently:
      coord_med_dw = c( -cos(angle_dw), sin(angle_dw) ) * rr_dw + circle.center.dw
      # compute starting and ending points of the line passing through the median
      v_dw = coord_med_dw - circle.center.dw
      u_dw = v_dw / sqrt(sum(v_dw^2))
      step_dw = (0.05*sqrt(sum(v_dw^2)))
      coord_start_dw = coord_med_dw - step_dw*u_dw
      coord_end_dw   = coord_med_dw + step_dw*u_dw
      # NOTE: use opposite sign for y-coordinate since the HDDS is in the 3-4th quadrants
      med_data_dw = data.frame(
         x = coord_start_dw[1],
         xend = coord_end_dw[1],
         y = -coord_start_dw[2],
         yend = -coord_end_dw[2]
      )
      hc = hc + geom_segment(data= med_data_dw, aes(x= x, y= y, xend= xend, yend= yend),
                             color="black", linetype="solid", size=1.10)
   }


   # add DATA points
   if(add.data.points.up){
      dots_up = data.up[data.up > bounds.up[1] & data.up < bounds.up[2]]
      # apply map from [bounds[1],bounds[2]] -> [-1,1]
      r_points_up = diameter.up/2-0.04
      # equispaced points on the semi-circle: same arc length
      xx = vector(mode = "numeric", NBINS.up+1)
      yy = vector(mode = "numeric", NBINS.up+1)
      stp_up = (pi * r_points_up) / NBINS.up       # arc length. NOTE: (pi * r) = length of semi-circle with radius r
      for(i in 1:(NBINS.up+1)){
         theta = (i-1)*stp_up / r_points_up     # get angle under arc
         xx[i] = r_points_up * cos(theta)    # get x-coordinate
         yy[i] = r_points_up * sin(theta)    # get y-coordinate
      }
      xx = sort(xx)
      yy = yy + circle.center.up[2]
      # Map dots into (xx,yy) coordinates, such that: min(dots)-> xx[NBINS], max(dots)-> xx[1]
      xend_pts_up = vector(mode = "numeric", length(dots_up))
      yend_pts_up = vector(mode = "numeric", length(dots_up))
      for(j in 1:length(dots_up)){
         ll = hist(dots_up[j], breaks = c(bounds.up[1]-1,seq(bounds.up[1], bounds.up[2], length.out=NBINS.up)), plot=F)$counts
         xend_pts_up[j] = xx[as.logical(ll)]
         yend_pts_up[j] = yy[as.logical(ll)]
      }
      if(jitter.data){
         xend_pts_up = xend_pts_up + rnorm(length(xend_pts_up),0,0.008)
         yend_pts_up = yend_pts_up + rnorm(length(yend_pts_up),0,0.008)
      }
      hc = hc + geom_point(aes(x = xend_pts_up, y = yend_pts_up), colour=color.data.up, size = 0.80)
   }
   if(add.data.points.dw){
      dots_dw = data.dw[data.dw > bounds.dw[1] & data.dw < bounds.dw[2]]
      # apply map from [bounds[1],bounds[2]] -> [-1,1]
      r_points_dw = diameter.dw/2-0.04
      # equispaced points on the semi-circle: same arc length
      xx = vector(mode = "numeric", NBINS.dw+1)
      yy = vector(mode = "numeric", NBINS.dw+1)
      stp_dw = (pi * r_points_dw) / NBINS.dw       # arc length. NOTE: (pi * r) = length of semi-circle with radius r
      for(i in 1:(NBINS.dw+1)){
         theta = (i-1)*stp_dw / r_points_dw     # get angle under arc
         xx[i] = r_points_dw * cos(theta)    # get x-coordinate
         yy[i] = r_points_dw * sin(theta)    # get y-coordinate
      }
      xx = sort(xx)
      yy = yy + circle.center.dw[2]
      # Map dots into (xx,yy) coordinates, such that: min(dots)-> xx[NBINS], max(dots)-> xx[1]
      xend_pts_dw = vector(mode = "numeric", length(dots_dw))
      yend_pts_dw = vector(mode = "numeric", length(dots_dw))
      for(j in 1:length(dots_up)){
         ll = hist(dots_dw[j], breaks = c(bounds.dw[1]-1,seq(bounds.dw[1], bounds.dw[2], length.out=NBINS.dw)), plot=F)$counts
         xend_pts_dw[j] = xx[as.logical(ll)]
         yend_pts_dw[j] = yy[as.logical(ll)]
      }
      if(jitter.data){
         xend_pts_dw = xend_pts_dw + rnorm(length(xend_pts_dw),0,0.008)
         yend_pts_dw = yend_pts_dw + rnorm(length(yend_pts_dw),0,0.008)
      }
      hc = hc + geom_point(aes(x = xend_pts_dw, y = yend_pts_dw), colour=color.data.dw, size = 0.80)
   }


   if(labels == TRUE){
      hc = hc +
         annotate("text", x=circle.center.up[1]+0.56*diameter.up, y=circle.center.up[2]+0.02*diameter.up,
                  label=round(max(BINS.up),2), size=labels.sz, color="black") +
         annotate("text", x=circle.center.up[1]-0.56*diameter.up, y=circle.center.up[2]+0.02*diameter.up,
                  label=round(min(BINS.up),2), size=labels.sz, color="black")

      # # add labels UP
      # hc = hc +
      #      annotate("text", x=circle.center.up[1]+0.57*diameter.up, y=circle.center.up[2]+0.06*diameter.up,
      #               label=round(max(BINS_up),2), size=5, color=colors.from.to.up[2]) +
      #      annotate("text", x=circle.center.up[1]-0.57*diameter.up, y=circle.center.up[2]+0.06*diameter.up,
      #               label=round(min(BINS_up),2), size=5, color=colors.from.to.up[2]) +
      #      annotate("text", x=circle.center.up[1], y=circle.center.up[2]+0.53*diameter.up,
      #               label=round((max(BINS_up)-min(BINS_up))/2,2), size=5, color=colors.from.to.up[2])
      # # add labels DOWN
      # hc = hc +
      #      annotate("text", x=circle.center.dw[1]+0.57*diameter.dw, y=circle.center.dw[2]-0.06*diameter.dw,
      #               label=round(max(BINS_dw),2), size=5, color=colors.from.to.dw[2]) +
      #      annotate("text", x=circle.center.dw[1]-0.57*diameter.dw, y=circle.center.dw[2]-0.06*diameter.dw,
      #               label=round(min(BINS_dw),2), size=5, color=colors.from.to.dw[2]) +
      #      annotate("text", x=circle.center.dw[1], y=circle.center.dw[2]-0.53*diameter.dw,
      #               label=round((max(BINS_dw)-min(BINS_dw))/2,2), size=5, color=colors.from.to.dw[2])
   }

   if(legend == TRUE){
      hc = hc +
         theme(axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
               axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
               panel.background=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(), plot.background=element_blank()) +
         # scale_colour_gradient(low="white", high="black", limits=c(0,1), name="Density") +
         # theme(legend.direction = "vertical", legend.position=c(0.98,0.5))
         scale_colour_gradient(low= colors.from.to.up[1], high= colors.from.to.up[2], limits=c(0,1), name="Density") +
         # theme(legend.direction = "horizontal", legend.position=c(0.45,0.04)) + #"bottom")
         scale_fill_gradient(low= colors.from.to.dw[1], high= colors.from.to.dw[2], limits=c(0,1), name="Density") #+
      # theme(legend.direction = "vertical", legend.position=c(0.98,0.5))
   }else{
      hc = hc +
         theme(axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
               axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
               panel.background=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(), plot.background=element_blank(), legend.position="none")
   }
   return(hc)
}

