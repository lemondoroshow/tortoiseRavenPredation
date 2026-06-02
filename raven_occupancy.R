rm(list = ls())
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
p.occ <- ncol(bbs_data$occ.covs) + 1
dist.bbs <- dist(bbs_data$coords)
mean.dist <- mean(dist.bbs)
min.dist <- min(dist.bbs)
max.dist <- max(dist.bbs)
inits <- list(alpha = rep(alpha.start, p.det),
              beta = rep(beta.start, p.occ),
              sigma.sq = sigma.sq.start, 
              phi = phi.start, 
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
save(out, file = paste("runs/bbs_spPGOcc_", chain, "_", 
                       Sys.time(), ".R", sep = ''))
