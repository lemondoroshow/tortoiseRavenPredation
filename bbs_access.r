rm(list = ls())
library(tidyverse)
library(lubridate)
library(sp)
library(raster)
library(FedData)
library(sf)
library(stars)
library(terra)
library(tidyterra)

# Get common raven AOU
spec <- read.csv("./data/bbs/SpeciesList.csv") %>%
  filter(English_Common_Name == "Common Raven")
aou_cc <- spec$AOU

# Import weather data
weather <- read.csv("./data/bbs/Weather.csv") |>
  unite('date', sep = '-', Year, Month, Day, remove = FALSE) %>%
  mutate(date = as.Date(date, tz = "America/New_York")) %>%
  mutate(julian = as.numeric(format(date, '%j')))

# Import BBS data, and filter (we will use 2024 as an "example year")
bbs_obs <- read.csv("./data/bbs/States/Arizona.csv") %>%
  bind_rows(read.csv("./data/bbs/States/Califor.csv")) %>%
  bind_rows(read.csv("./data/bbs/States/Nevada.csv")) %>%
  bind_rows(read.csv("./data/bbs/States/Utah.csv")) %>%
  left_join(read.csv("./data/bbs/Routes.csv"), 
            by = c('Route', 'CountryNum', 'StateNum')) %>%
  filter(Year == 2024) %>%
  mutate(RouteNum = paste0(StateNum, formatC(Route, 2, flag = "0"))) %>%
  dplyr::select(RouteDataID, Latitude, Longitude, AOU, starts_with("Count")) %>%
  dplyr::select(-CountryNum) %>%
  complete(AOU, nesting(RouteDataID, Latitude, Longitude)) %>%
  replace(is.na(.), 0) %>%
  filter(AOU == aou_cc) %>%
  left_join(weather, by = c("RouteDataID"))

# Isolate BBS data in area of interest
mojave <- vect("./data/shapefiles/mojaveDesert/MojaveEcoregion_TNC_UTM83.shp") |>
  project("WGS84")
bbs_obs_mojave <- vect(bbs_obs, geom = c("Longitude", "Latitude"), 
                       crs = crs(mojave)) |>
  mask(mojave) |>
  as.data.frame(geom = "XY") |>
  rename(Longitude = x, Latitude = y) |>
  relocate(Latitude, .after = RouteDataID) |>
  relocate(Longitude, .after = Latitude)

# Create detection-nondetection matrix
y <- as.matrix(bbs_obs_mojave[, c('Count10', 'Count20', 'Count30', 'Count40', 'Count50')])
y <- ifelse(y > 0, 1, 0)

# Compile detection covariates
det_covs <- list(
  day = c(scale(bbs_obs_mojave$julian)),
  day.2 = c(scale(bbs_obs_mojave$julian) ^ 2),
  tod = c(scale(bbs_obs_mojave$StartTime)),
  obs = as.numeric(factor(bbs_obs_mojave$ObsN))
)

# Compile coordinates
coords <- bbs_obs_mojave[, c('Latitude', 'Longitude')]
coords_sp <- data.frame(coords, val = apply(y, 1, max)) %>%
  arrange(Longitude, Latitude)
coordinates(coords_sp) <- ~Longitude + Latitude
proj4string(coords_sp) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'