rm(list = ls())
library(terra)
library(tidyverse)
library(tmap)

# Import shapefiles
mojave <- vect("./data/shapefiles/RU/2011RecoveryUnitsDissolved.shp") |>
  project("WGS84")
conus <- vect("./data/shapefiles/CONUS/conus_ard_grid.shp") |>
  project("WGS84")

# For 2024, get landfills in a vector for CONUS
yoi <- 2024
landfills <- read.csv("./data/lmop_database.csv") |>
  dplyr::select(c(Longitude, Latitude, Landfill.ID, Year.Landfill.Opened,
                  Landfill.Closure.Year)) |>
  dplyr::rename(ID = Landfill.ID, Open = Year.Landfill.Opened, 
                Close = Landfill.Closure.Year) |>
  dplyr::filter(Open <= yoi & Close > yoi) |>
  vect(crs = "NAD83", geom = c("Longitude", "Latitude")) |>
  project(mojave) |>
  crop(conus)

# Rasterize landfills by distance
landfill_rast <- rast(ext(mojave), res = 0.1, crs = crs(mojave))
landfill_dist <- distance(landfill_rast, landfills, rasterize = TRUE) |>
  mask(mojave)