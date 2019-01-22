# Check whether a soap film smoother boundary and knots make sense

# see Readme.md for details on how to use this

library(mgcv)
library(rgeos)
library(sp)
source("autocrunch.R")

soap_check <- function(bnd, knots=NULL, data=NULL, plot=TRUE,
                       tol=sqrt(.Machine$double.eps)){

  ## check that the boundary makes sense
  # check that boundary is a list
  stopifnot(is.list(bnd))

  # check that the boundary (or boundary part) have x and y elements
  lapply(bnd, function(x) stopifnot(c("x","y") %in% names(x)))
  # each boundary part must have at least 4 elements!
  lapply(bnd, function(x) stopifnot(length(x$x)>3, length(x$y)>3))

  # check that the boundary loops are actually loops
  check_ends <- function(x, tol){
    all.equal(c(x$x[1], x$y[1]),
              c(x$x[length(x$y)], x$y[length(x$y)]),tolerance=tol)
  }
  end_check <- unlist(lapply(bnd, check_ends, tol=tol))
  end_check_logical <- is.character(end_check)
  if(any(end_check_logical)){
    stop(paste("Boundary loop(s)",which(end_check_logical),
               "don't have identical start & end points",collapse=" "))
  }

  islands <- FALSE
  # check for intersections
  if(length(bnd)>1){
    inds <- combn(1:length(bnd),2)

    ## make the bnds into polys here
    make_bnd_poly <- function(bnd){

      bnd$x[length(bnd$x)] <-  bnd$x[1]
      bnd$y[length(bnd$y)] <-  bnd$y[1]
      SpatialPolygons(list(Polygons(list(Polygon(bnd)),ID=1)))
    }
    bnd_poly <- lapply(bnd, make_bnd_poly)

    # function to see if two polygons intersect
    intersects <- function(this.ind, bnd){
      poly1 <- bnd[[this.ind[1]]]
      poly2 <- bnd[[this.ind[2]]]

      gIntersects(poly1, poly2)
    }
    # apply over all the combinations
    inter <- apply(inds, 2, intersects, bnd=bnd_poly)

    if(any(inter)){
      # get the index for the prospective "outer" loop
      outer_ind <- which.max(unlist(lapply(bnd_poly, gArea)))
      outer_bnd <- bnd_poly[[outer_ind]]

      other_bnd <- bnd_poly
      other_bnd[[outer_ind]] <- NULL

      # is everything else inside that?
      islands <- unlist(lapply(other_bnd, gWithin, spgeom2=outer_bnd))

      if(!all(islands)){
        stop(paste("Polygon parts",
                   paste0(apply(inds[,inter, drop=FALSE], 2, paste0,
                                collapse=" and "),
                          collapse=", "), "intersect"))
      }

      islands <- all(islands)
    }
  }

  ## plot what the boundary is
  # highlighting the area to be modelled
  if(plot){
    # colourblind-safe colours from colorbrewer2 "qualitative" map
    red <- "#d95f02"
    # if the boundary is only 1 part, plotting is rather easier
    if(!islands){
      plot(bnd[[1]], type="l", main="Red indicates soap film surface", asp=1)
      lapply(bnd, polygon, col=red)
    }else{
      outer_bnd <- bnd[[outer_ind]]
      other_bnd <- bnd
      other_bnd[[outer_ind]] <- NULL
      plot(outer_bnd, type="n", main="Red indicates soap film surface", asp=1)
      # plot the outer loop
      polygon(outer_bnd, col=red)
      # plot the other polygons on top in white
      lapply(other_bnd, polygon, col="white")
    }
  }

  # function to check if points are inside the boundary
  point_check <- function(bnd, x, y, type){

    if(length(bnd)>1){
      # inSide doesn't deal with edge points very well
      # but does handle multiple rings better
      inout <- inSide(bnd, x, y)
    }else{
      # use sp::point.in.polygon
      # see ?point.in.polygon for returned codes, 1 is inside
      pip <- function(bnd, x, y){
        point.in.polygon(x, y, bnd$x, bnd$y)==1
      }
      # apply over the parts of the polygon
      inout <- pip(bnd[[1]], x, y)
    }
    if(!all(inout)){
      stop(paste(type, paste(which(!inout),collapse=", "),
                 "are outside the boundary."))
    }
  }

  ## check the knots
  if(!is.null(knots)){
    # check that the points have x and y elements
    stopifnot(c("x","y") %in% names(knots))
    point_check(bnd, knots$x, knots$y, "Knots")
    crunch_ind <- autocruncher(bnd, knots)
    if(!is.null(crunch_ind)){
      stop(paste0("Knots ", paste(crunch_ind, collapse=", "),
                 "are outside the boundary."))
    }
    if(plot){
      points(knots, col="#1b9e77", pch=19)
      if(!is.null(crunch_ind)) points(knots[crunch_ind, ], col="#1b9e77", pch=4)
    }
  }

  ## check the data
  if(!is.null(data)){
    # check that the points have x and y elements
    stopifnot(c("x","y") %in% names(data))
    point_check(bnd, data$x, data$y, "Data points")
    if(plot) points(data, col="#7570b3", pch=19)
  }

  return(TRUE)
}
