

# Check whether a soap film smoother boundary and knots make sense

# see Readme.md for details on how to use this


#' @import dplyr
#' @import graphics
#' @import mgcv
#' @import sf
#' @import utils
#' @export


soap_check <- function(bnd, knots = NULL,
                       data = NULL,
                       plot = TRUE,
                       tol = sqrt(.Machine$double.eps),
                       x_name = "x",
                       y_name = "y"){

  ## check that the boundary makes sense
  # check that boundary is a list
  stopifnot(is.list(bnd))

  # check that the boundary (or boundary part) have x and y elements
  lapply(bnd, function(x) stopifnot(c("x","y") %in% names(x)))
  # each boundary part must have at least 4 elements!
  lapply(bnd, function(x) stopifnot(length(x$x) > 3, length(x$y) > 3))

  # check that the boundary loops are actually loops
  check_ends <- function(x, tol){
    all.equal(c(x$x[1], x$y[1]),
              c(x$x[length(x$y)],
                x$y[length(x$y)]),
              tolerance = tol)
  }
  end_check <- unlist(lapply(bnd, check_ends, tol = tol))
  end_check_logical <- is.character(end_check)
  if(any(end_check_logical)){
    stop(paste("Boundary loop(s)", which(end_check_logical),
               "don't have identical start & end points", collapse = " "))
  }

  islands <- FALSE
  # check for intersections
  if(length(bnd) > 1){
    inds <- utils::combn(1:length(bnd), 2)

    # make the bnds into polys here
    make_bnd_poly <- function(bnd, crs = NULL){
      # make last point become the first point to complete the polygons
      bnd$x[length(bnd$x) + 1] <-  bnd$x[1]
      bnd$y[length(bnd$y) + 1] <-  bnd$y[1]

      # convert each bnd list object into sf object
      bnd <- as.data.frame(bnd) |>
        sf::st_as_sf(coords = c("x", "y")) |>
        dplyr::mutate(
          id = 1
        ) |>
        dplyr::group_by(id) |>
        dplyr::summarise(do_union = FALSE) |>
        sf::st_cast("POLYGON")
    }
    # loop over an dconvert each object in the list
    bnd_poly <- lapply(bnd, make_bnd_poly)

    # function to see if two polygons intersect
    intersects <- function(this.ind, bnd){
      poly1 <- bnd[[this.ind[1]]]
      poly2 <- bnd[[this.ind[2]]]

      sf::st_intersects(poly1, poly2, sparse = FALSE)
    }
    # apply over all the combinations
    inter <- apply(inds, 2, intersects, bnd = bnd_poly)

    if(any(inter)){
      # get the index for the prospective "outer" loop
      outer_ind <- which.max(unlist(lapply(bnd_poly, sf::st_area)))
      outer_bnd <- bnd_poly[[outer_ind]]

      other_bnd <- bnd_poly
      other_bnd[[outer_ind]] <- NULL

      # is everything else inside that?
      islands <- unlist(lapply(other_bnd, sf::st_within, y = outer_bnd,
                               sparse = FALSE))

      if(!all(islands)){
        warning(paste("Polygon parts",
                      paste0(apply(inds[,inter, drop=FALSE], 2, paste0,
                                   collapse=" and "),
                             collapse=", "), "intersect"))
      }

      islands <- all(islands)
    }
  }

  ## plot what the boundary is
  # highlighting the area to be modelled
}
# if(all(islands)){
# highlighting the area to be modeled
if(plot){
  # colourblind-safe colours from colorbrewer2 "qualitative" map
  red <- "#d95f02"
  # if the boundary is only 1 part, plotting is rather easier
  if(!islands){
    plot(bnd[[1]], type = "l", main = "Red indicates soap film surface",
         asp = 1,
         ylab = "y",
         xlab = "x")
    lapply(bnd, graphics::polygon, col = red)
  } else {
    outer_bnd <- bnd[[outer_ind]]
    other_bnd <- bnd
    other_bnd[[outer_ind]] <- NULL

    if(all(islands)){

      plot(outer_bnd, type = "n",
           main = "Red indicates soap film surface",
           ylab = "y",
           xlab = "x",
           asp = 1)
      # plot the outer loop
      graphics::polygon(outer_bnd, col = red)
      # plot the other polygons on top in white
      lapply(other_bnd, graphics::polygon, col = "white")
    } else {
      # if the boundary is only 1 part, plotting is rather easier
      plot(
        outer_bnd,
        type = "n",
        main = "Red indicates main boundary of film surface",
        ylab = "y",
        xlab = "x",
        asp = 1)

      # plot the outer loop
      graphics::polygon(outer_bnd, col = red)
      # sapply(1:length(other_bnd))
      col_n <- length(other_bnd)
      # if(is.null(col)) {
      cols <- terrain.colors(n = col_n)
      # } else
      # cols <- rainbow(n = col_n)
      # other_bnd <-
      # plot the other polygons on top in white
      lapply(1:length(other_bnd), function(x) graphics::polygon(other_bnd[[x]],
                                                                col = cols[[x]]))
    }
  }
}

  # function to check if points are inside the boundary
  point_check <- function(bnd, x, y, type){

    if(length(bnd)>1){
      # inSide doesn't deal with edge points very well
      # but does handle multiple rings better
      inout <- mgcv::inSide(bnd, x, y)
    }else{
      # use sp::point.in.polygon
      # see ?point.in.polygon for returned codes, 1 is inside
      pip <- function(bnd, x, y){
        sp::point.in.polygon(x, y, bnd$x, bnd$y)==1
      }
      # apply over the parts of the polygon
      inout <- pip(bnd[[1]], x, y)
    }
    if(!all(inout)){
      warning(paste(type, paste(which(!inout),collapse=", "),
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
