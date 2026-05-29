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

# Define year of interest
yoi = 2024

# Get common raven AOU
spec <- read.csv("./data/bbs/SpeciesList.csv") %>%
  filter(English_Common_Name == "Common Raven")
aou_cc <- spec$AOU

# Import weather data
weather <- read.csv("./data/bbs/Weather.csv") |>
  unite("date", sep = "-", Year, Month, Day, remove = FALSE) %>%
  mutate(date = as.Date(date, tz = "America/New_York")) %>%
  mutate(julian = as.numeric(format(date, "%j")))

# Import BBS data
bbs_obs <- read.csv("./data/bbs/States/Arizona.csv") %>%
  bind_rows(read.csv("./data/bbs/States/Califor.csv")) %>%
  bind_rows(read.csv("./data/bbs/States/Nevada.csv")) %>%
  bind_rows(read.csv("./data/bbs/States/Utah.csv")) %>%
  left_join(read.csv("./data/bbs/Routes.csv"), 
            by = c("Route", "CountryNum", "StateNum")) %>%
  filter(Year == yoi) %>%
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
y <- as.matrix(bbs_obs_mojave[, c("Count10", "Count20", "Count30", "Count40", "Count50")])
y <- ifelse(y > 0, 1, 0)

# Compile detection covariates
det_covs <- list(
  day = c(scale(bbs_obs_mojave$julian)),
  day.2 = c(scale(bbs_obs_mojave$julian) ^ 2),
  tod = c(scale(bbs_obs_mojave$StartTime)),
  obs = as.numeric(factor(bbs_obs_mojave$ObsN))
)

# Compile coordinates
coords <- bbs_obs_mojave[, c("Longitude", "Latitude")]

# Extract NDVIs
ndvis <- c()
for (row in 1:dim(bbs_obs_mojave)[1]) {
  
  # Import NDVI data
  date <- format(bbs_obs_mojave$date[row])
  ndvi <- rast(paste0("./data/ndvi/processed/", date, ".tif"))
  
  # Get NDVI at coords=
  ndvi_ext <- terra::extract(ndvi, coords[row,])[,2]
  ndvis <- c(ndvis, ndvi_ext)
}

# Add NDVIs to data frame
bbs_obs_mojave$NDVI <- ndvis

# Extract elevations, add to data frame
elev <- rast("./data/elevation/all_time.tif")
elev_ext <- terra::extract(elev, coords)
bbs_obs_mojave$Elevation <- unlist(elev_ext[,2])

# Extract road density, add to data frame
roads <- rast(paste0("./data/roads/", as.character(yoi), ".tif"))
roads_ext <- terra::extract(roads, coords)
bbs_obs_mojave$RoadDensity <- unlist(roads_ext[,2])

# Extract transmission line distance, add to data frame
lines <- rast("./data/transmissionLines/all_time.tif")
lines_ext <- terra::extract(lines, coords)
bbs_obs_mojave$TransmissionLineDistance <- unlist(lines_ext[,2])

# Extract % impervious, add to data frame
impervious <- rast(paste0("./data/impervious/processed/", as.character(yoi), ".tif"))
impervious_ext <- terra::extract(impervious, coords)
bbs_obs_mojave$PercentImpervious <- unlist(impervious_ext[,2])
