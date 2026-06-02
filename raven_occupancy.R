library(coda)
library(sf)
library(spOccupancy)

# Import data
load("./data/bbs_and_cov_data_bundle.R")

# Get coordinates and project to AEA for CONUS
coords_aea <- data.frame(bbs_data$coords) |>
  st_as_sf(coords = c("Longitude", "Latitude"), 
           crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") |>
  st_transform(crs = "EPSG:5070") |>
  st_coordinates()
bbs_data$coords <- coords_aea / 1000
