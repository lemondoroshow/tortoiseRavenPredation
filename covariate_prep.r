library(sf)
library(terra)
library(tidyterra)
library(tmap)

# Import Mojave shapefile
mojave <- vect("./data/shapefiles/mojaveDesert/MojaveEcoregion_TNC_UTM83.shp")

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

# Import all roads for 2024
files <- list.files("./data/shapefiles/roads/2024", pattern = "*.shp$",
                    recursive = TRUE, full.names = TRUE)
roads <- vect()
for (file in files) {
  print(file)
  roads <- union(roads, vect(file) |> project(mojave))
}

# Crop to an expanded Mojave box
roads_mjv <- project(roads, mojave) |>
  crop(ext(mojave))

# Create density raster
roads_rast <- rast(ext(roads_mjv), res = 10000, crs = crs(roads_mjv))
roads_area <- cellSize(roads_rast)
roads_dens <- rasterizeGeom(roads_mjv, roads_rast, "length") / roads_area

# Crop and mask to Mojave; rename raster layer
roads_dens_mjv <- mask(roads_dens, mojave) |>
  crop(mojave) |>
  tidyterra::rename(density = length)

# Map density
roads_dens_map <- tm_shape(roads_dens_mjv, bbox = bbox) +
  tm_raster(col = "density", col.scale = tm_scale(values = "brewer.yl_or_br"),
            col.legend = tm_legend(position = tm_pos_out("right", "center"))) +
  tm_basemap("Esri.WorldImagery") +
  tm_title("Road density in the Mojave Desert, 2024")
tmap_save(roads_dens_map, "./plots/mojave_road_density.png")

# Import transmission lines
lines <- vect(paste0("./data/shapefiles/transmissionLines/",
                     "Electric_Power_Transmission_Lines_A.shp")) |>
  project(mojave)

# Create distance raster
lines_rast <- rast(ext(mojave), res = 10000, crs = crs(mojave))
lines_dist <- distance(lines_rast, lines) |>
  mask(mojave) |>
  tidyterra::rename(distance = lyr.1)

# Map distance
lines_dist_map <- tm_shape(lines_dist, bbox = bbox) +
  tm_raster(col = "distance", col.scale = tm_scale(values = "brewer.or_rd"),
            col.legend = tm_legend(position = tm_pos_out("right", "center"))) +
  tm_basemap("Esri.WorldImagery") +
  tm_title("Distance from nearest transmission line in the Mojave Desert")
tmap_save(lines_dist_map, "./plots/mojave_transmission_line_distance.png")
