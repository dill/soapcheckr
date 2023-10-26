#' Check if boundary, knots, and data can be modeled using a soap-film smoother
#'
#' Setting up a soap film smoother is often a hard and frustrating process.
#' This function checks that the boundary, knots, and data that you feed to
#' a soap-film smoother are “correct”. The function will plot  the boundary,
#' knots and data that are trying to be modeled to ensure that they are are
#' appropriate
#'
#' @param bnd A list with sub-lists that will be the boundary the soap-film smoother
#' will smooth in. Coordinates need to be in meaningful units.
#' For example, coordinates could be in UTMs as UTMs rely on metres.
#' @param knots A dataframe with two columns that are the coordinates of the knots
#' that are to be supplied to the soap-film smoother. Coordinates need to be in
#' meaningful units. For example, coordinates could be in UTMs as UTMs rely on metres.
#' @param data A dataframe with two columns that are the coordinates of the data
#' that are to be supplied to the soap-film smoother. Coordinates need to be in
#' meaningful units. For example, coordinates could be in UTMs as UTMs rely on metres.
#' @param plot logical if plot of boundary, knots, and/or data should be plotted.
#' Default is `TRUE`
#' @param tol Tolerance value to check if boundaries are complete polygons.
#' Sometimes tolerance needs to be increased (e.g., tol = 1e-6)
#' @param x_name Column name of x coordinate for the knots and/or data object
#' @param y_name Column name of y coordinate for the knots and/or data object
#'
#' @return TRUE or FALSE depending on whether the boundary will be able to used
#' by a soap-film smoother. Addtionally if supplying knots and/or data. It will warn
#' the user which knots and/or data fall too close or outside the boundary.
#'
#' @examples
#' # library(mgcv)
#' # fsb <- list(fs.boundary())
#' # soap_check(fsb)
#'
#' soap_check(sissabagama_lake_ls)
#'
#' @import dplyr
#' @import graphics
#' @import grDevices
#' @import mgcv
#' @import sf
#' @import utils
#' @export
#'


soap_check <- function(bnd,
                       knots = NULL,
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
  # knots location check
  #  bnd supplied has to be list, data has to be dataframe object

  point_check <- function(bnd, data, type){
    if(length(bnd) > 1){
      # inSide doesn't deal with edge points very well
      # but does handle multiple rings better
      if(!all(names(data) %in% c("x", "y"))) {
        data$x <- data[[x_name]]
        data[[x_name]] <- NULL
        data$y <- data[[y_name]]
        data[[y_name]] <- NULL
      }
      x <- data$x
      y <- data$y
      inout <- mgcv::inSide(bnd, x = x, y = y)
    } else{
      # use sp::point.in.polygon
      # need to use sf points here
      # see ?point.in.polygon for returned codes, 1 is inside

      pip <- function(bnd, data){
        # convert bndry into sf object

        bnd_comb <- bnd |>
          as.data.frame() |>
          st_as_sf(coords = c("x", "y")) |>
          dplyr::mutate(
            id = 1
          ) |>
          dplyr::group_by(id) |>
          dplyr::summarise(do_union = FALSE) |>
          st_cast("POLYGON")

        if(!all(names(data) %in% c("x", "y"))) {
          data$x <- data[[x_name]]
          data[[x_name]] <- NULL
          data$y <- data[[y_name]]
          data[[y_name]] <- NULL
        }
        knt_sf <- data |>
          st_as_sf(coords = c("x", "y"))

        # to replace sp you need to put knots into
        # sp::point.in.polygon(x, y, bnd$x, bnd$y)==1
        st_intersects(bnd_comb, knt_sf, sparse = FALSE)[TRUE]
      }
      # apply over the parts of the polygon
      inout <- pip(bnd = bnd[1], data)
    }
    if(!all(inout)){
      warning(paste(type, paste(which(!inout), collapse = ", "),
                    " are outside the boundary."))
    }
  }

  # check the knots
  if(!is.null(knots)){
    # check that the points have x and y elements
    # stopifnot(c("x", "y") %in% names(knots))
    if(length(bnd) > 1) {
      suppressWarnings(
        point_check(bnd, knots, "Knots")
      )
    } else {
      point_check(bnd, knots, "Knots")
    }

    crunch_ind <- autocruncher(bnd, knots,
                               xname = x_name,
                               yname = y_name)

    if(!is.null(crunch_ind)){
      warning(paste0("Knots ", paste(crunch_ind, collapse = ", "),
                     "are outside the boundary."))
    }
    #DMILL colour "#1b9e77"
    if(plot){
      points(knots,
             col = "black", pch = 21)
      if(!is.null(crunch_ind)){

        points(knots[crunch_ind, ],
               col = "black", pch = 4)
      }
    }
  }
  #     # check the data
  if(!is.null(data)){
    stopifnot(c("x", "y") %in% names(data))

    # check that the points have x and y elements
    point_check(bnd, data, "Data points")
    if(plot)
      #DMILL colour #7570b3
      points(data, col = "black",
             pch = 19)
  }

  if(length(bnd) > 1 & all(islands)){
    return(TRUE)
  }
  if(length(bnd) > 1 & !all(islands)){
    return(FALSE)
  }
  if(length(bnd) %in% 1){
    return(TRUE)
  }
}

