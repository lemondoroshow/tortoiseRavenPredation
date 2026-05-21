library(terra)
library(tmap)

# Import Mojave shapefile
mojave <- vect("./data/shapefiles/mojaveDesert/MojaveEcoregion_TNC_UTM83.shp")

# Import all roads for 2024
files <- list.files("./data/shapefiles/roads/2024", pattern = "*.shp$", recursive = TRUE, full.names = TRUE)
roads <- vect()
for (file in files) {
  print(file)
  roads <- union(roads, vect(file) |> project(mojave))
}

# Crop and mask to Mojave 
roads <- project(roads, mojave) |>
  mask(mojave) |>
  crop(mojave)
