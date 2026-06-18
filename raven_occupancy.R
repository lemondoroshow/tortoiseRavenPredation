rm(list = ls())
library(coda)
library(dplyr)
library(sf)
library(spOccupancy)
library(terra)
library(tmap)

# Import data
load("./data/bbs_and_cov_data_bundle.R")

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
occ.formula <- ~ elevation + ndvi + impervious + roads + lines + 
                  impervious * roads + impervious * lines
det.formula <- ~ day + day.2 + tod + (1 | obs)
p.det <- length(bbs_data$det.covs)
p.occ <- ncol(bbs_data$occ.covs) + 1 + 2 # covs + intercept + interactions
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
               n.report = n.report)
run_str <- paste0(chain, "_", Sys.time())
save(out, file = paste("runs/bbs_spPGOcc_", run_str, ".R", sep = ''))

# Get points in Mojave to predict at
mojave <- vect("./data/shapefiles/RU/2011RecoveryUnitsDissolved.shp") |>
  project("EPSG:5070")
mojave_points <- rast(ext(mojave), res = 7500, crs = crs(mojave))

# Open and resample covariates
ndvi <- rast("./data/ndvi/processed/2024-05-30.tif") |> # Random date
  project("EPSG:5070") |>
  resample(mojave_points) |>
  scale()
elev <- rast("./data/elevation/all_time.tif") |>
  project("EPSG:5070") |>
  resample(mojave_points) |>
  scale()
roads <- rast("./data/roads/2024.tif") |>
  project("EPSG:5070") |>
  resample(mojave_points) |>
  scale()
lines <- rast("./data/transmissionLines/all_time.tif") |>
  project("EPSG:5070") |>
  resample(mojave_points) |>
  scale()
impervious <- rast("./data/impervious/processed/2024.tif") |>
  project("EPSG:5070") |>
  resample(mojave_points) |>
  scale()

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
  ) |> left_join(
    as.data.frame(impervious * roads, xy = TRUE),
    by = c("x", "y")
  ) |> left_join(
    as.data.frame(lines * roads, xy = TRUE),
    by = c("x", "y")
  ) |>
  na.omit()
colnames(covs) <- c("x", "y", "elevation", "ndvi", "impervious", "roads", "lines",
                    "impervious_x_roads", "lines_x_roads")
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
