rm(list = ls())
library(elevatr)
library(sf)
library(terra)
library(tidyterra)
library(tidyverse)
library(tmap)

#### Preparation ####

# Set year of interest
yoi <- 2024

# Import Mojave shapefile
mojave <- vect("./data/shapefiles/RU/2011RecoveryUnitsDissolved.shp") |>
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

# Import BBS data
bbs_obs <- read.csv("./data/bbs/States/Arizona.csv") |>
  bind_rows(read.csv("./data/bbs/States/Califor.csv")) |>
  bind_rows(read.csv("./data/bbs/States/Nevada.csv")) |>
  bind_rows(read.csv("./data/bbs/States/Utah.csv")) |>
  left_join(read.csv("./data/bbs/Routes.csv"), 
            by = c("Route", "CountryNum", "StateNum")) |>
  mutate(RouteNum = paste0(StateNum, formatC(Route, 2, flag = "0"))) |>
  dplyr::select(RouteDataID, Latitude, Longitude, AOU, starts_with("Count")) |>
  dplyr::select(-CountryNum) |>
  complete(AOU, nesting(RouteDataID, Latitude, Longitude)) %>%
  replace(is.na(.), 0) |> # magrittr pipe on line 48 allows for this format
  filter(AOU == aou_cc) |>
  left_join(weather, by = c("RouteDataID"))

# Isolate BBS data in area of interest
obs_dates <- vect(bbs_obs, geom = c("Longitude", "Latitude"), 
                  crs = crs(mojave)) |>
  mask(mojave) |>
  as.data.frame(geom = "XY") |>
  rename(Longitude = x, Latitude = y) |>
  dplyr::select(Latitude, Longitude, date, Year) |>
  dplyr::rename(Date = date)
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
writeRaster(roads_dens_mjv, paste0("./data/roads/", as.character(yoi), ".tif"),
            overwrite = TRUE)

# Map density
roads_dens_map <- tm_shape(roads_dens_mjv, bbox = bbox) +
  tm_raster(col = "density", col.scale = tm_scale(values = "brewer.yl_or_br"),
            col.legend = tm_legend(position = tm_pos_out("right", "center"))) +
  tm_basemap("Esri.WorldImagery") +
  tm_title(paste0("Road density in the Mojave Desert, ", as.character(yoi)))
tmap_save(roads_dens_map, "./plots/mojave_road_density.png")

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
writeRaster(lines_dist, "./data/transmissionLines/all_time.tif", 
            overwrite = TRUE)

# Map distance
lines_dist_map <- tm_shape(lines_dist, bbox = bbox) +
  tm_raster(col = "distance", col.scale = tm_scale(values = "brewer.or_rd"),
            col.legend = tm_legend(position = tm_pos_out("right", "center"))) +
  tm_basemap("Esri.WorldImagery") +
  tm_title("Distance from nearest transmission line in the Mojave Desert")
tmap_save(lines_dist_map, "./plots/mojave_transmission_line_distance.png")

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
                                 as.character(yoi), ".tif"), overwrite = TRUE)
}

#### Elevation ####

# Create template raster for downloading data
rast_template <- rasterize(mojave, rast(crs = crs(mojave), 
                                        ext = ext(mojave), res = 0.01))

# Access data
elev <- rast(get_elev_raster(locations = rast_template, z = 10))

# Morph for Mojave
elev_mjv <- project(elev, mojave) |>
  mask(mojave) |>
  crop(mojave)

# Export
writeRaster(elev_mjv, "./data/elevation/all_time.tif", overwrite = TRUE)

#### NDVI -- all years ####

# Import Mojave shapefile
mojave <- vect("./data/shapefiles/RU/2011RecoveryUnitsDissolved.shp") |>
  project("WGS84")

# Iterate through years
years <- 2001:2024
years <- years[! years %in% c(2020)]
for (yoi in years) {
  
  # Import all NDVI data
  files <- list.files(paste0("./data/ndvi/raw/", as.character(yoi), "/"), 
                      full.names = TRUE)
  for (file in files) {
    
    # Transform date
    date <- unlist(strsplit(file, "_"))[5]
    date <- as.Date(date, format = "%Y%m%d")
    
    # Import data; mask, crop
    ndvi <- rast(file) |>
      terra::subset("NDVI") |>
      project(mojave) |>
      crop(bbox)
    
    # Export data
    writeRaster(ndvi, paste0("./data/ndvi/processed/", format(date), ".tif"),
                overwrite = TRUE)
    }
}

#### % Impervious -- all years ####

# Import Mojave shapefile
mojave <- vect("./data/shapefiles/RU/2011RecoveryUnitsDissolved.shp") |>
  project("WGS84")

# Iterate through years
years <- 2001:2024
years <- years[! years %in% c(2020)]
for (yoi in years) {
    
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
                                     as.character(yoi), ".tif"), overwrite = TRUE)
    }
}
#### Landfills -- all years ####

# Import shapefiles
mojave <- vect("./data/shapefiles/RU/2011RecoveryUnitsDissolved.shp") |>
  project("WGS84")
conus <- vect("./data/shapefiles/CONUS/conus_ard_grid.shp") |>
  project("WGS84")

# Iterate through years
years <- 2001:2024
years <- years[! years %in% c(2020)]
for (yoi in years) {
  
  # Read data and filter
  landfills <- read.csv("./data/lmop_database.csv") |>
    dplyr::select(c(Longitude, Latitude, Landfill.ID, Year.Landfill.Opened,
                    Landfill.Closure.Year)) |>
    dplyr::rename(ID = Landfill.ID, Open = Year.Landfill.Opened, 
                  Close = Landfill.Closure.Year) |>
    # Technically speaking, we can't be certain whether we are ignoring some 
    # landfills that are open with no expected landfill date. However, these 
    # landfills are ignored to avoid inconsistency with the same situation 
    # potentially occurring in past data
    dplyr::filter(Open <= yoi & Close > yoi) |>
    vect(crs = "NAD83", geom = c("Longitude", "Latitude")) |>
    project(mojave) |>
    crop(conus)
  
  # Rasterize landfills by distance
  landfill_rast <- rast(ext(mojave), res = 0.1, crs = crs(mojave))
  landfill_dist <- distance(landfill_rast, landfills, rasterize = TRUE) |>
    mask(mojave) |>
    tidyterra::rename(distance = layer)
  
  # Export raster
  writeRaster(landfill_dist, paste0("./data/landfills/", yoi, ".tif"),
              overwrite = TRUE)
  
}
#### Towers -- all years ####

# Import shapefiles
mojave <- vect("./data/shapefiles/RU/2011RecoveryUnitsDissolved.shp") |>
  project("WGS84")
conus <- vect("./data/shapefiles/CONUS/conus_ard_grid.shp") |>
  project("WGS84")

# Iterate through years
years <- 2001:2024
years <- years[! years %in% c(2020)]
for (yoi in years) {
  
  # Open tower data
  co <- read.table("./data/towers/raw/CO.dat", sep = "|", 
                   header = FALSE, fill = TRUE) |>
    dplyr::select(c("V5", "V11", "V16"))
  ra <- read.table("./data/towers/raw/RA.dat", sep = "|", 
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
  
  # Write raster
  writeRaster(tower_dist, paste0("./data/towers/processed/", yoi, ".tif"))
  
}