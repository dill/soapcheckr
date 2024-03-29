#' Sissabagama Lake located in northern Wisconsin
#'
#' Sissabagama Lake located in northern Wisconsin. The original shapefile was derived from
#' the Wisconsin DNR's watershed database. Latitude and longitude are in UTM zone 15N.
#'
#'
#' @format Simple feature collection with 1 features and 23 fields:
#'  \describe{
#'    \item{geometry}{POLYGON}
#' }
#'
"sissabagama_lake_sf"


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


#' Sampled bathymetry points for Sissabagama Lake
#'
#' The latitude and longitude of sampled depth derived from bathymetric map from the
#' [Wisconsin DNR](https://apps.dnr.wi.gov/lakes/lakepages/LakeDetail.aspx?wbic=2393500) for
#' Sissabagama Lake in northern Wisconsin.
#'
#' @format A dataframe containing 445 rows and 3 variables:
#'  \describe{
#'    \item{x}{The longitude of each sample}
#'    \item{y}{The latitude of each sample}
#'    \item{depth}{The latitude of each sample}
#' }
#'
"sissabagama_bath"
