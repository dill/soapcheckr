## ---- include = FALSE-----------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  fig.dim = c(7, 5), 
  comment = "#>"
)

## ----setup, message = FALSE, warning = FALSE------------------------------------------------------------------------------------------
{
  library(broom.mixed)
  library(fitdistrplus)
  library(dplyr)
  library(ggplot2)
  library(gratia)
  library(mgcv)
  library(purrr)
  library(soapcheckr)
  library(sf)
}

## ----ramsays horseshoe----------------------------------------------------------------------------------------------------------------
fsb <- list(fs.boundary())

## ----check ramsays horseshoe----------------------------------------------------------------------------------------------------------
soap_check(fsb)

## ----build fsb knots, message = FALSE-------------------------------------------------------------------------------------------------
# create knots 
knots <- expand.grid(x = seq(min(fsb[[1]]$x), 
                             max(fsb[[1]]$x), len = 15),
                     y = seq(min(fsb[[1]]$y) + 0.05,
                             max(fsb[[1]]$y), len = 10))
x <- knots$x
y <- knots$y

# identify the knots that are outside the boundary 
ind <- inSide(fsb, x = x, y = y)
# remove knots outside the boundary 
knots <- knots[ind, ]

## ----build fake data, message = FALSE-------------------------------------------------------------------------------------------------
set.seed(0)
n <- 600

# Our x and y data 
x <- runif(n) * 5 - 1
y <- runif(n) * 2 - 1

# create our response variable 
z <- fs.test(x, y, b = 1)

## remove outsiders

ind <- inSide(fsb, x = x, y = y) 

z <- z + rnorm(n) * 0.3 ## add noise

# create the data we want to model 
dat <- data.frame(z = z[ind],
                  x = x[ind],
                  y = y[ind])

## ----check our knots------------------------------------------------------------------------------------------------------------------
soap_check(fsb, knots = knots)

## ----remove knots that are offending--------------------------------------------------------------------------------------------------
crunch_index <- autocruncher(fsb, knots, k = 30)
crunch_index

# remove knots that are problematic
knots <- knots[-crunch_index, ] 

## -------------------------------------------------------------------------------------------------------------------------------------
soap_check(fsb, knots = knots)

## -------------------------------------------------------------------------------------------------------------------------------------
dat_2 <- dat[, 2:3]
soap_check(fsb, data = dat_2)

## ----run GAM--------------------------------------------------------------------------------------------------------------------------
m <- gam(z ~ s(x, y, k = 30 , bs = "so",
               xt = list(bnd = fsb)),
         knots = knots, 
         data = dat)


## ----check main effects of model, eval = FALSE----------------------------------------------------------------------------------------
#  anova(m)

## ----draw model effects---------------------------------------------------------------------------------------------------------------
draw(m)

## ----check model fit------------------------------------------------------------------------------------------------------------------
appraise(m)

## ----sissabagama lake crs, message = FALSE--------------------------------------------------------------------------------------------
sissabagama_lake_sf <- sissabagama_lake_sf %>% 
  st_transform(crs = 32615)

## ----sissabagama lake pnt sf, message = FALSE-----------------------------------------------------------------------------------------
bnd_pt_sf <- sissabagama_lake_sf %>%
  dplyr::select(geometry) %>%
  st_cast("MULTIPOINT") %>%
  mutate(
    id = 1:nrow(.)
  )

## ----sissabagama lake pnt, message = FALSE, warning = FALSE---------------------------------------------------------------------------
bnd_pt <- bnd_pt_sf %>%
  split(.$id) %>%
  purrr::map(~ st_cast(.x, "POINT") %>%
               mutate(
                 x = st_coordinates(.)[,"X"],
                 y = st_coordinates(.)[,"Y"]
               ) %>%
               st_drop_geometry() %>% 
               dplyr::select(-id)
  )

## ----create list of lists boundary, message = FALSE-----------------------------------------------------------------------------------
nr <- 1:length(bnd_pt)

sissabagama_bnd_ls <- lapply(nr, function(n) as.list.data.frame(bnd_pt[[n]]))

## ----check more complex boundary------------------------------------------------------------------------------------------------------
soap_check(sissabagama_bnd_ls)

## ----make knot grid, message = FALSE--------------------------------------------------------------------------------------------------
lake_grid <- sissabagama_lake_sf %>%
  st_make_grid(cellsize = 200, square = TRUE, what = "centers") %>%
  st_as_sf() 

st_geometry(lake_grid) <- "geometry"

## ----remove all knots that are outside boundary, message = FALSE, warning = FALSE-----------------------------------------------------
lake_intesects <- st_intersection(sissabagama_lake_sf, lake_grid)

## ----make knot dataframe, message = FALSE---------------------------------------------------------------------------------------------
lake_knots <- lake_intesects %>%
  mutate(
    lon = st_coordinates(.)[,"X"],
    lat = st_coordinates(.)[,"Y"]
  ) %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  dplyr::select(lon, lat)

## ----check knots using soap_check-----------------------------------------------------------------------------------------------------
soap_check(sissabagama_bnd_ls, knots = lake_knots, 
           x_name = "lon", y_name = "lat")

## -------------------------------------------------------------------------------------------------------------------------------------
nrow(lake_knots)
crunch_ind <- autocruncher(sissabagama_bnd_ls, lake_knots, k = 37, 
                           xname = "lon", yname = "lat")
crunch_ind

# remove knots that are problematic
lake_knots <- lake_knots[-crunch_ind, ] 

## ----recheck knots--------------------------------------------------------------------------------------------------------------------
soap_check(sissabagama_bnd_ls, knots = lake_knots, 
           x_name = "lon", y_name = "lat")


## ----remove depth for soap_check------------------------------------------------------------------------------------------------------
sissabagama_bath_pt <- sissabagama_bath %>% 
  dplyr::select(-depth)

soap_check(sissabagama_bnd_ls, data = sissabagama_bath_pt)

## ----check skewness and kurtosis------------------------------------------------------------------------------------------------------
depths <- sissabagama_bath$depth
descdist(depths)

## ----check distr----------------------------------------------------------------------------------------------------------------------
fit_gamma <- fitdist(depths, distr = "gamma", method = "mme")
plot(fit_gamma)

## ----add f to our boundary list-------------------------------------------------------------------------------------------------------
names(lake_knots) <- c("x", "y")

sissabagama_bnd_ls <- lapply(nr,
                             function(n)
                               sissabagama_bnd_ls[[n]] <- c(
                                 sissabagama_bnd_ls[[n]],
                                 list(f = rep(0, length(sissabagama_bnd_ls[[n]]$x))
                                 )
                               )
)

## ----model our depth, warning = FALSE-------------------------------------------------------------------------------------------------
m1 <- gam(depth ~ s(x, y,
                    bs = "so",
                    xt = list(bnd = sissabagama_bnd_ls)),
          family = Gamma(link = "identity"),
          knots = lake_knots,
          data = sissabagama_bath)

## ----check main effects---------------------------------------------------------------------------------------------------------------
anova(m1)
summary(m1)

## ----check partial effects------------------------------------------------------------------------------------------------------------
draw(m1)

## ----check model fit for m1-----------------------------------------------------------------------------------------------------------
appraise(m1)

## ----create 10m grid to efficiently predict over the boundary-------------------------------------------------------------------------
lake_pred <- sissabagama_lake_sf %>%
  st_make_grid(cellsize = 10, square = TRUE, what = "centers") %>% 
  st_as_sf() 
st_geometry(lake_pred) <- "geometry"

## ----remove all points that fall outside the boundary, warning = FALSE----------------------------------------------------------------
lake_pred <- st_intersection(lake_pred, sissabagama_lake_sf) %>% 
  dplyr::select(geometry)

## ----convert to a dataframe-----------------------------------------------------------------------------------------------------------
lake_pred_df <- lake_pred %>% 
  mutate(
    x = st_coordinates(.)[,"X"], 
    y = st_coordinates(.)[,"Y"], 
  ) %>% 
  st_drop_geometry()


## ----predict using augment from broom.mixed, message = FALSE--------------------------------------------------------------------------
pred <- augment(m1, newdata = lake_pred_df)
pred <- pred %>% 
  mutate(
    lower = (.fitted - 1.96 * .se.fit),
    higher = (.fitted + 1.96 * .se.fit)
  )
glimpse(pred)

## ----plot our predicted depths using ggplot-------------------------------------------------------------------------------------------
ggplot() +
  geom_raster(data = pred, aes(x = x, y = y, fill = .fitted)) +
  geom_sf(data = sissabagama_lake_sf, fill = NA, colour = "black") +
  scale_fill_viridis_c(name = "Depth (m)",
                       trans = "reverse",
                       breaks = rev(seq(0, 60, 15))
  ) + 
  theme_void(
    base_size = 15
  ) + 
  theme(
    legend.background = element_blank(),
    legend.position = c(0.98, 0.82),
  ) + 
  guides(fill = guide_colourbar(
    frame.linewidth = 0.3,
    ticks.colour = 'black', 
    frame.colour = 'black')) + 
  labs(x = "Longitude", 
       y = "Latitude")


