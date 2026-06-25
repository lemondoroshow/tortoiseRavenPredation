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
  mask(mojave) |>
  tidyterra::rename(distance = layer)

# Open tower data
co <- read.table("./data/towers/CO.dat", sep = "|", 
                 header = FALSE, fill = TRUE) |>
  dplyr::select(c("V5", "V11", "V16"))
ra <- read.table("./data/towers/RA.dat", sep = "|", 
                 header = FALSE, fill = TRUE) |>
  dplyr::select(c("V5", "V13", "V14"))
co$V5 <- as.character(co$V5)
towers <- left_join(co, ra, by = "V5") |>
  mutate(V13 = na_if(V13, "")) |>
  drop_na(V11, V16, V13)

# Wrangle data formats
towers$YearConstructed <- substr(towers$V13, 7, 10) |> 
  as.integer()
dismantled <- substr(towers$V14, 7, 10) |> 
  as.integer()
dismantled[is.na(dismantled)] <- 9999
towers$YearDismantled <- dismantled
ys <- towers$V11 / 3600
xs <- -1 * towers$V16 / 3600
towers$Latitude <- ys
towers$Longitude <- xs

# Convert to geom object
towers <- as.data.frame(towers) |>
  dplyr::filter(YearConstructed <= yoi & YearDismantled > yoi) |>
  vect(crs = "NAD83", geom = c("Longitude", "Latitude")) |>
  project(mojave) |>
  crop(conus)

# Rasterize towers by distance
tower_rast <- rast(ext(mojave), res = 0.1, crs = crs(mojave))
tower_dist <- distance(tower_rast, towers, rasterize = TRUE) |>
  mask(mojave) |>
  tidyterra::rename(distance = layer)
