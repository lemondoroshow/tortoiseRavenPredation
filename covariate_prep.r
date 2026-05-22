library(terra)
library(tidyterra)
library(tmap)

# Import Mojave shapefile
mojave <- vect("./data/shapefiles/mojaveDesert/MojaveEcoregion_TNC_UTM83.shp")

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
roads_rast <- rast(ext(roads_mjv), res = 2000, crs = crs(roads_mjv))
roads_area <- cellSize(roads_rast)
roads_dens <- rasterizeGeom(roads_mjv, roads_rast, "length") / roads_area

# Crop and mask to Mojave; rename raster layer
roads_dens_mjv <- mask(roads_dens, mojave) |>
  crop(mojave) |>
  tidyterra::rename(density = length)

# Map density
roads_dens_map <- tm_shape(roads_dens_mjv) +
  tm_raster(col = "density", col.scale = tm_scale(values = "brewer.yl_or_br"),
            col.legend = tm_legend(position = tm_pos_out("right", "center"))) +
  tm_basemap("Esri.WorldGrayCanvas") +
  tm_title("Road density in the Mojave Desert")
tmap_save(roads_dens_map, "./plots/mojave_road_density.png")
