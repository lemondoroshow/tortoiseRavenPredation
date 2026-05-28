library(sf)
library(terra)
library(tidyterra)
library(tidyverse)
library(tmap)

#### Preparation ####

# Set year of interest
yoi = 2024

# Import Mojave shapefile
mojave <- vect("./data/shapefiles/mojaveDesert/MojaveEcoregion_TNC_UTM83.shp") |>
  project("WGS84")

# Create bounding box for plots
bbox <- st_bbox(mojave)
x_range <- bbox$xmax - bbox$xmin
y_range <- bbox$ymax - bbox$ymin
bbox[1] <- bbox[1] - 0.10 * x_range
bbox[2] <- bbox[2] - 0.10 * y_range
bbox[3] <- bbox[3] + 0.10 * x_range
bbox[4] <- bbox[4] + 0.10 * y_range
bbox <- st_as_sfc(bbox)
st_crs(bbox) <- crs(mojave)

# Get common raven AOU
spec <- read.csv("./data/bbs/SpeciesList.csv") |>
  filter(English_Common_Name == "Common Raven")
aou_cc <- spec$AOU

# Import weather data
weather <- read.csv("./data/bbs/Weather.csv") |>
  unite("date", sep = "-", Year, Month, Day, remove = FALSE) |>
  mutate(date = as.Date(date, tz = "America/New_York"))

# Import BBS data, and filter
bbs_obs <- read.csv("./data/bbs/States/Arizona.csv") |>
  bind_rows(read.csv("./data/bbs/States/Califor.csv")) |>
  bind_rows(read.csv("./data/bbs/States/Nevada.csv")) |>
  bind_rows(read.csv("./data/bbs/States/Utah.csv")) |>
  left_join(read.csv("./data/bbs/Routes.csv"), 
            by = c("Route", "CountryNum", "StateNum")) |>
  filter(AOU == aou_cc) |>
  left_join(weather, by = c("RouteDataID"))

# Isolate BBS data in area of interest; get dates of observations
mojave <- vect("./data/shapefiles/mojaveDesert/MojaveEcoregion_TNC_UTM83.shp") |>
  project("WGS84")
obs_dates <- vect(bbs_obs, geom = c("Longitude", "Latitude"), 
                  crs = crs(mojave)) |>
  mask(mojave) |>
  as.data.frame(geom = "XY") |>
  rename(Longitude = x, Latitude = y) |>
  dplyr::select(Latitude, Longitude, date, Year.x) |>
  dplyr::rename(Year = Year.x, Date = date)
write.csv(obs_dates, "./data/bbs_obs_mojave.csv")

# Transform dates
dates <- (filter(obs_dates, Year == yoi))$Date |>
  lapply(format, format = "%Y%m%d") |>
  unlist() |>
  sort()

# Clean up
rm(bbs_obs, spec, weather, x_range, y_range)

#### Roads ####

# Import all roads for year of interest
files <- list.files(paste0("./data/shapefiles/roads/", as.character(yoi)), 
                    pattern = "*.shp$", recursive = TRUE, full.names = TRUE)
roads <- vect()
for (file in files) {
  print(file)
  roads <- terra::union(roads, vect(file) |> project(mojave))
}

# Crop to an expanded Mojave box
roads_mjv <- project(roads, mojave) |>
  crop(ext(mojave))

# Create density raster
roads_rast <- rast(ext(roads_mjv), res = 0.1, crs = crs(roads_mjv))
roads_area <- cellSize(roads_rast)
roads_dens <- rasterizeGeom(roads_mjv, roads_rast, "length") / roads_area

# Crop and mask to Mojave; rename raster layer
roads_dens_mjv <- mask(roads_dens, mojave) |>
  crop(mojave) |>
  tidyterra::rename(density = length)

# Export raster
writeRaster(roads_dens_mjv, paste0("./data/roads/", as.character(yoi), ".tif"))

# Map density
roads_dens_map <- tm_shape(roads_dens_mjv, bbox = bbox) +
  tm_raster(col = "density", col.scale = tm_scale(values = "brewer.yl_or_br"),
            col.legend = tm_legend(position = tm_pos_out("right", "center"))) +
  tm_basemap("Esri.WorldImagery") +
  tm_title(paste0("Road density in the Mojave Desert, ", as.character(yoi)))
# tmap_save(roads_dens_map, "./plots/mojave_road_density.png")

#### Transmission lines ####

# Import transmission lines
lines <- vect(paste0("./data/shapefiles/transmissionLines/",
                     "Electric_Power_Transmission_Lines_A.shp")) |>
  project(mojave)

# Create distance raster
lines_rast <- rast(ext(mojave), res = 0.1, crs = crs(mojave))
lines_dist <- distance(lines_rast, lines, rasterize = TRUE) |>
  mask(mojave) |>
  tidyterra::rename(distance = layer)

# Export raster
writeRaster(lines_dist, "./data/transmissionLines/all_time.tif")

# Map distance
lines_dist_map <- tm_shape(lines_dist, bbox = bbox) +
  tm_raster(col = "distance", col.scale = tm_scale(values = "brewer.or_rd"),
            col.legend = tm_legend(position = tm_pos_out("right", "center"))) +
  tm_basemap("Esri.WorldImagery") +
  tm_title("Distance from nearest transmission line in the Mojave Desert")
# tmap_save(lines_dist_map, "./plots/mojave_transmission_line_distance.png")

#### NDVI ####

# Import all NDVI data
files <- list.files("./data/ndvi/raw/", pattern = as.character(yoi), full.names = TRUE)
for (file in files) {
  
  # Transform date
  date <- unlist(strsplit(file, "_"))[5]
  date <- as.Date(date, format = "%Y%m%d")
  
  # Import data; mask, crop
  ndvi <- rast(file) |>
    terra::subset("NDVI") |>
    project(mojave) |>
    mask(mojave) |>
    crop(mojave)
  
  # Export data
  writeRaster(ndvi, paste0("./data/ndvi/processed/", format(date), ".tif"),
              overwrite = TRUE)

}

#### % Impervious ####

# Import all % impervious files
files <- list.files("./data/impervious/raw/", pattern = as.character(yoi), full.names = TRUE)
for (file in files) {
  
  # Import data; mask, crop
  impervious <- rast(file) |>
    aggregate(10) |> # I'd love to keep the original resolution, but it's too big
    project(mojave) |>
    mask(mojave) |>
    crop(mojave)
    
  # Export raster
  writeRaster(impervious, paste0("./data/impervious/processed/",
                                 as.character(yoi), ".tif"))
}
