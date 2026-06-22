rm(list = ls())
library(coda)
library(dplyr)
library(sf)
library(spOccupancy)
library(terra)
library(tmap)

set.seed(04021990)

#### One year ####

yoi <- 2024

# Import data
load(paste0("./data/bundles/", as.character(yoi), "_data_bundle.R"))

# Get coordinates and project to AEA for CONUS
coords_aea <- data.frame(bbs_data$coords) |>
  st_as_sf(coords = c("Longitude", "Latitude"), 
           crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") |>
  st_transform(crs = "EPSG:5070") |>
  st_coordinates()
bbs_data$coords <- coords_aea / 1000

# Set initial values
p_file_name <- "pfile-sp-1"
chain <- unlist(strsplit(p_file_name, "-"))[3]
p_file <- read.table(paste0("priors/", p_file_name),
                     sep = ' ', header = FALSE)
alpha.start <- p_file[p_file[, 1] == 'alpha', 2]
beta.start <- p_file[p_file[, 1] == 'beta', 2]
sigma.sq.start <- p_file[p_file[, 1] == 'sigma.sq', 2]
phi.start <- p_file[p_file[, 1] == 'phi', 2]
w.start <- p_file[p_file[, 1] == 'w', 2]

# Specify model
occ.formula <- ~ elevation + ndvi + impervious + roads + lines
det.formula <- ~ day + day.2 + tod + (1 | obs)
p.det <- length(bbs_data$det.covs)
p.occ <- ncol(bbs_data$occ.covs) + 1 # covs + intercept
dist.bbs <- dist(bbs_data$coords)
mean.dist <- mean(dist.bbs)
min.dist <- min(dist.bbs)
max.dist <- max(dist.bbs)
inits <- list(alpha = rep(alpha.start, p.det),
              beta = rep(beta.start, p.occ),
              sigma.sq = sigma.sq.start, 
              # phi = phi.start, # This line was giving me an error?
              w = rep(w.start, nrow(bbs_data$y)),
              z = apply(bbs_data$y, 1, max, na.rm = TRUE))
priors <- list(beta.normal = list(mean = rep(0, p.occ),
                                  var = rep(2.72, p.occ)),
               alpha.normal = list(mean = rep(0, p.det),
                                   var = rep(2.72, p.det)), 
               phi.unif = c(3 / max.dist, 3 / min.dist), 
               sigma.sq.ig = c(2, 5))
batch.length <- 25
n.batch <- 2000
n.burn <- 10000
n.thin <- 20
n.report <- 20
tuning <- list(phi = 1)

# Run model
out <- spPGOcc(occ.formula = occ.formula,
               det.formula = det.formula,
               data = bbs_data,
               inits = inits,
               batch.length = batch.length,
               n.batch = n.batch, 
               tuning = tuning,
               priors = priors,
               n.omp.threads = 1,
               verbose = TRUE,
               NNGP = TRUE, 
               n.neighbors = 10,
               cov.model = 'gaussian', 
               n.burn = n.burn,
               n.thin = n.thin,
               n.report = n.report,
               k.fold = 10, 
               k.fold.threads = 10)
run_str <- paste0(as.character(yoi), "_chain", chain, "_", 
                  format(Sys.time(), "%Y-%m-%dT%H:%M:%S%Z"))
save(out, file = paste("runs/bbs_spPGOcc_", run_str, ".R", sep = ''))

# Get points in Mojave to predict at
mojave <- vect("./data/shapefiles/RU/2011RecoveryUnitsDissolved.shp") |>
  project("EPSG:5070")
mojave_points <- rast(ext(mojave), res = 7500, crs = crs(mojave))

# Open and resample covariates
ndvi <- rast("./data/ndvi/processed/2024-05-30.tif") |> # Random date
  project("EPSG:5070") |>
  resample(mojave_points) |>
  scale(center = bbs_data$stats$ndvi[1], 
        scale = bbs_data$stats$ndvi[2])
elev <- rast("./data/elevation/all_time.tif") |>
  project("EPSG:5070") |>
  resample(mojave_points) |>
  scale(center = bbs_data$stats$elevation[1], 
        scale = bbs_data$stats$elevation[2])
roads <- rast("./data/roads/2024.tif") |>
  project("EPSG:5070") |>
  resample(mojave_points) |>
  scale(center = bbs_data$stats$roads[1], 
        scale = bbs_data$stats$roads[2])
lines <- rast("./data/transmissionLines/all_time.tif") |>
  project("EPSG:5070") |>
  resample(mojave_points) |>
  scale(center = bbs_data$stats$lines[1], 
        scale = bbs_data$stats$lines[2])
impervious <- rast("./data/impervious/processed/2024.tif") |>
  project("EPSG:5070") |>
  resample(mojave_points) |>
  scale(center = bbs_data$stats$impervious[1], 
        scale = bbs_data$stats$impervious[2])

# Put covariates together
covs <- as.data.frame(elev, xy = TRUE) |>
  left_join(
    as.data.frame(ndvi, xy = TRUE),
    by = c("x", "y")
  ) |> left_join(
    as.data.frame(impervious, xy = TRUE),
    by = c("x", "y")
  ) |> left_join(
    as.data.frame(roads, xy = TRUE),
    by = c("x", "y")
  ) |> left_join(
    as.data.frame(lines, xy = TRUE),
    by = c("x", "y")
  # ) |> left_join(
  #   as.data.frame(impervious * roads, xy = TRUE),
  #   by = c("x", "y")
  # ) |> left_join(
  #   as.data.frame(lines * roads, xy = TRUE),
  #   by = c("x", "y")
  ) |>
  na.omit()
colnames(covs) <- c("x", "y", "elevation", "ndvi", "impervious", "roads", "lines")
                    #"impervious_x_roads", "lines_x_roads")
n_locs <- dim(covs)[1]

# Prepare covariate matrix
covs$Intercept <- rep(1, n_locs)
X.0 <- dplyr::select(covs, -c("x", "y")) |>
  as.matrix()

# Prepare coordinate matrix
coords.0 <- dplyr::select(covs, c("x", "y")) |>
  as.matrix()

# Predict results
res <- predict(out, X.0, coords.0)
occ <- res$psi.0.samples |>
  colMeans()
occ_loc <- matrix(nrow = n_locs, ncol = 3)
occ_loc[,1:2] <- coords.0
occ_loc[,3] <- occ
colnames(occ_loc) <- c("x", "y", "psi")

# Create prediction raster
occ_rast <- as.data.frame(occ_loc) |>
  rast(crs = "EPSG:5070")
writeRaster(occ_rast, paste0("./results/", run_str, ".tif"), overwrite = TRUE)

#### All years ####

# Iterate through years
years <- 2001:2024
years <- years[! years %in% c(2020)]
for (yoi in years) {
  
  # Import data
  load(paste0("./data/bundles/", as.character(yoi), "_data_bundle.R"))
  
  # Get coordinates and project to AEA for CONUS
  coords_aea <- data.frame(bbs_data$coords) |>
    st_as_sf(coords = c("Longitude", "Latitude"), 
             crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") |>
    st_transform(crs = "EPSG:5070") |>
    st_coordinates()
  bbs_data$coords <- coords_aea / 1000
  
  # Set initial values
  p_file_name <- "pfile-sp-1"
  chain <- unlist(strsplit(p_file_name, "-"))[3]
  p_file <- read.table(paste0("priors/", p_file_name),
                       sep = ' ', header = FALSE)
  alpha.start <- p_file[p_file[, 1] == 'alpha', 2]
  beta.start <- p_file[p_file[, 1] == 'beta', 2]
  sigma.sq.start <- p_file[p_file[, 1] == 'sigma.sq', 2]
  phi.start <- p_file[p_file[, 1] == 'phi', 2]
  w.start <- p_file[p_file[, 1] == 'w', 2]
  
  # Specify model
  occ.formula <- ~ elevation + ndvi + impervious + roads + lines
  det.formula <- ~ day + day.2 + tod + (1 | obs)
  p.det <- length(bbs_data$det.covs)
  p.occ <- ncol(bbs_data$occ.covs) + 1 # covs + intercept
  dist.bbs <- dist(bbs_data$coords)
  mean.dist <- mean(dist.bbs)
  min.dist <- min(dist.bbs)
  max.dist <- max(dist.bbs)
  inits <- list(alpha = rep(alpha.start, p.det),
                beta = rep(beta.start, p.occ),
                sigma.sq = sigma.sq.start, 
                # phi = phi.start, # This line was giving me an error?
                w = rep(w.start, nrow(bbs_data$y)),
                z = apply(bbs_data$y, 1, max, na.rm = TRUE))
  priors <- list(beta.normal = list(mean = rep(0, p.occ),
                                    var = rep(2.72, p.occ)),
                 alpha.normal = list(mean = rep(0, p.det),
                                     var = rep(2.72, p.det)), 
                 phi.unif = c(3 / max.dist, 3 / min.dist), 
                 sigma.sq.ig = c(2, 5))
  batch.length <- 25
  n.batch <- 2000
  n.burn <- 10000
  n.thin <- 20
  n.report <- 20
  tuning <- list(phi = 1)
  
  # Run model
  out <- spPGOcc(occ.formula = occ.formula,
                 det.formula = det.formula,
                 data = bbs_data,
                 inits = inits,
                 batch.length = batch.length,
                 n.batch = n.batch, 
                 tuning = tuning,
                 priors = priors,
                 n.omp.threads = 1,
                 verbose = TRUE,
                 NNGP = TRUE, 
                 n.neighbors = 10,
                 cov.model = 'gaussian', 
                 n.burn = n.burn,
                 n.thin = n.thin,
                 n.report = n.report,
                 k.fold = 10, 
                 k.fold.threads = 10)
  run_str <- paste0(as.character(yoi), "_chain", chain, "_", 
                    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%Z"))
  save(out, file = paste("runs/bbs_spPGOcc_", run_str, ".R", sep = ''))
  
  # Get points in Mojave to predict at
  mojave <- vect("./data/shapefiles/RU/2011RecoveryUnitsDissolved.shp") |>
    project("EPSG:5070")
  mojave_points <- rast(ext(mojave), res = 7500, crs = crs(mojave))
  
  # Average NDVIs for year
  all_ndvis <- list.files("./data/ndvi/processed/", pattern = as.character(yoi),
                          full.names = TRUE)
  ndvi_rasts <- c()
  for (ndvi_path in all_ndvis) {
    ndvi_rasts <- c(ndvi_rasts, rast(ndvi_path))
  }
  mean_ndvi <- sprc(ndvi_rasts) |>
    na.omit() |>
    terra::mosaic()
  
  # Open and resample covariates
  ndvi <- project(mean_ndvi, "EPSG:5070") |>
    resample(mojave_points) |>
    scale(center = bbs_data$stats$ndvi[1], 
          scale = bbs_data$stats$ndvi[2])
  elev <- rast("./data/elevation/all_time.tif") |>
    project("EPSG:5070") |>
    resample(mojave_points) |>
    scale(center = bbs_data$stats$elevation[1], 
          scale = bbs_data$stats$elevation[2])
  roads <- rast("./data/roads/2024.tif") |> # Use 2024 roads for now
    project("EPSG:5070") |>
    resample(mojave_points) |>
    scale(center = bbs_data$stats$roads[1], 
          scale = bbs_data$stats$roads[2])
  lines <- rast("./data/transmissionLines/all_time.tif") |>
    project("EPSG:5070") |>
    resample(mojave_points) |>
    scale(center = bbs_data$stats$lines[1], 
          scale = bbs_data$stats$lines[2])
  impervious <- rast(paste0("./data/impervious/processed/", 
                            as.character(yoi), ".tif")) |>
    project("EPSG:5070") |>
    resample(mojave_points) |>
    scale(center = bbs_data$stats$impervious[1], 
          scale = bbs_data$stats$impervious[2])
  
  # Put covariates together
  covs <- as.data.frame(elev, xy = TRUE) |>
    left_join(
      as.data.frame(ndvi, xy = TRUE),
      by = c("x", "y")
    ) |> left_join(
      as.data.frame(impervious, xy = TRUE),
      by = c("x", "y")
    ) |> left_join(
      as.data.frame(roads, xy = TRUE),
      by = c("x", "y")
    ) |> left_join(
      as.data.frame(lines, xy = TRUE),
      by = c("x", "y")
      # ) |> left_join(
      #   as.data.frame(impervious * roads, xy = TRUE),
      #   by = c("x", "y")
      # ) |> left_join(
      #   as.data.frame(lines * roads, xy = TRUE),
      #   by = c("x", "y")
    ) |>
    na.omit()
  colnames(covs) <- c("x", "y", "elevation", "ndvi", "impervious", "roads", "lines")
  #"impervious_x_roads", "lines_x_roads")
  n_locs <- dim(covs)[1]
  
  # Prepare covariate matrix
  covs$Intercept <- rep(1, n_locs)
  X.0 <- dplyr::select(covs, -c("x", "y")) |>
    as.matrix()
  
  # Prepare coordinate matrix
  coords.0 <- dplyr::select(covs, c("x", "y")) |>
    as.matrix()
  
  # Predict results
  res <- predict(out, X.0, coords.0)
  occ <- res$psi.0.samples |>
    colMeans()
  occ_loc <- matrix(nrow = n_locs, ncol = 3)
  occ_loc[,1:2] <- coords.0
  occ_loc[,3] <- occ
  colnames(occ_loc) <- c("x", "y", "psi")
  
  # Create prediction raster
  occ_rast <- as.data.frame(occ_loc) |>
    rast(crs = "EPSG:5070")
  writeRaster(occ_rast, paste0("./results/", run_str, ".tif"), overwrite = TRUE)
  
}

#### Visualization ####

# Import Mojave shapefile
mojave <- vect("./data/shapefiles/RU/2011RecoveryUnitsDissolved.shp") |>
  project("WGS84")

# Create bounding box for plots
bbox <- st_bbox(mojave)
x_range <- bbox$xmax - bbox$xmin
y_range <- bbox$ymax - bbox$ymin
bbox[1] <- bbox[1] - 0.10 * x_range
bbox[2] <- bbox[2] - 0.13 * y_range
bbox[3] <- bbox[3] + 0.13 * x_range
bbox[4] <- bbox[4] + 0.10 * y_range
bbox <- st_as_sfc(bbox)
st_crs(bbox) <- crs(mojave)

# Load TCAs
tcas <- terra::vect("./data/shapefiles/TCA/USFWS_DesertTortoise_TCAs.shp") |>
  project(mojave)

# Iterate through files
files <- list.files("./results/")
for (f in files) {
  
  # Import data
  yoi <- substr(f, 1, 4)
  res <- rast(paste0("./results/", f)) |>
    project(mojave)
  
  # Map results
  res_map <- tm_shape(res, bbox = bbox) +
    tm_raster("psi", col.scale = tm_scale(
      values = "crest", n = 5
    ), col.legend = tm_legend(
      position = tm_pos_in("right", "bottom")
    )) +
    tm_basemap("Esri.WorldGrayCanvas") +
    tm_title(paste0("C. corax occupancy probability in ", yoi)) +
    tm_shape(tcas) +
    tm_borders(col = "black")
  
  # Save
  tmap_save(res_map, paste0("./plots/occupancy/", yoi, ".png"))
}
