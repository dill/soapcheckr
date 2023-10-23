#' The boundary list of Sissabagama Lake
#'
#' A list with sub-lists containing the latitude and longitude of each polygon that
#' make up Sissabagama Lake in northern Wisconsin.
#'
#' @format A list containing 5 sub-lists with each sub-list containing two vectors, x and y.
#'  \describe{
#'    \item{x}{The longitude of each polygon, name needs to be x}
#'    \item{y}{The latitude of each poloygon, name needs to be y}
#' }
#'
"sissabagama_lake_ls"


#' Knots for Sissabagama Lake
#'
#' The latitude and longitude of knots evenly spaced 200 m away from each other inside
#' Sissabagama Lake in northern Wisconsin.
#'
#' @format A dataframe containing 420 rows and 2 variables:
#'  \describe{
#'    \item{lon}{The longitude of each knot}
#'    \item{lat}{The latitude of each knot}
#' }
#'
"sissabagama_lake_knots"


#' Knots for Sissabagama Lake created through using a grid
#'
#' The latitude and longitude of knots evenly spaced 200 m away from each other
#' covering the boundary box of Sissabagama lake. Knots included are outside the
#' boundary of the lake.
#'
#' @format A dataframe containing 1,125 rows and 2 variables:
#'  \describe{
#'    \item{lon}{The longitude of each knot}
#'    \item{lat}{The latitude of each knot}
#' }
#'
"sissabagama_lake_grid_knots"


#' The boundary list of example overlapping polygons
#'
#' A list with sub-lists containing the x and y coordinates of 3 overlapping polygons.
#' This list of boundaries is used as an example to demonstrate how soap_check functions
#' when polygons overlap.
#'
#' @format A list containing 3 sub-lists with each sub-list containing two vectors, x and y.
#'  \describe{
#'    \item{x}{The x coordinates of each polygon, name needs to be x}
#'    \item{y}{The y coordinates of each poloygon, name needs to be y}
#' }
#'
"overlap_polygons"
