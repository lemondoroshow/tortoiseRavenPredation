rm(list = ls())
library(car)
library(terra)
library(tidyverse)

#### Data prep -- largely taken from bbs_access.R ####

# Set year of interest
yoi <- 2024

# Get common raven AOU
spec <- read.csv("./data/bbs/SpeciesList.csv") |>
  filter(English_Common_Name == "Common Raven")
aou_cc <- spec$AOU

# Import weather data
weather <- read.csv("./data/bbs/Weather.csv") |>
  unite("date", sep = "-", Year, Month, Day, remove = FALSE) |>
  mutate(date = as.Date(date, tz = "America/New_York")) |>
  mutate(julian = as.numeric(format(date, "%j")))

# Import BBS data
bbs_obs <- read.csv("./data/bbs/States/Arizona.csv") |>
  bind_rows(read.csv("./data/bbs/States/Califor.csv")) |>
  bind_rows(read.csv("./data/bbs/States/Nevada.csv")) |>
  bind_rows(read.csv("./data/bbs/States/Utah.csv")) |>
  left_join(read.csv("./data/bbs/Routes.csv"),
            by = c("Route", "CountryNum", "StateNum")) |>
  filter(Year == yoi) |>
  mutate(RouteNum = paste0(StateNum, formatC(Route, 2, flag = "0"))) |>
  dplyr::select(RouteDataID, Latitude, Longitude, AOU, starts_with("Count")) |>
  dplyr::select(-CountryNum) |>
  complete(AOU, nesting(RouteDataID, Latitude, Longitude)) %>%
  replace(is.na(.), 0) |> # magrittr pipe on line 37 allows for this format
  filter(AOU == aou_cc) |>
  left_join(weather, by = c("RouteDataID"))

# Isolate BBS data in area of interest
mojave <- vect("./data/shapefiles/RU/2011RecoveryUnitsDissolved.shp") |>
  project("WGS84")
bbs_obs_mojave <- vect(bbs_obs, geom = c("Longitude", "Latitude"),
                       crs = crs(mojave)) |>
  mask(mojave) |>
  as.data.frame(geom = "XY") |>
  rename(Longitude = x, Latitude = y) |>
  relocate(Latitude, .after = RouteDataID) |>
  relocate(Longitude, .after = Latitude)
coords <- bbs_obs_mojave[, c("Longitude", "Latitude")]

# Create data frame with total counts
total_counts <- dplyr::select(bbs_obs_mojave, c(
  Count10, Count20, Count30, Count40, Count50
)) |>
  as.matrix() |>
  rowSums()
data <- dplyr::select(bbs_obs_mojave, c(Latitude, Longitude, date)) |>
  mutate(Counts = scale(total_counts))

# Extract NDVIs
ndvis <- c()
for (row in 1:dim(data)[1]) {

  # Import NDVI data
  date <- format(data$date[row])
  ndvi <- rast(paste0("./data/ndvi/processed/", date, ".tif"))

  # Get NDVI at coords
  ndvi_ext <- terra::extract(ndvi, coords[row,])[,2]
  ndvis <- c(ndvis, ndvi_ext)
}
data$NDVI <- ndvis |> scale()

# Extract elevations, add to data frame
elev <- rast("./data/elevation/all_time.tif")
elev_ext <- terra::extract(elev, coords)
data$Elevation <- unlist(elev_ext[,2]) |>
  scale()

# Extract road density, add to data frame
roads <- rast(paste0("./data/roads/", as.character(yoi), ".tif"))
roads_ext <- terra::extract(roads, coords)
data$RoadDensity <- unlist(roads_ext[,2]) |>
  scale()

# Extract transmission line distance, add to data frame
lines <- rast("./data/transmissionLines/all_time.tif")
lines_ext <- terra::extract(lines, coords)
data$TransmissionLineDistance <- unlist(lines_ext[,2]) |>
  scale()

# Extract % impervious, add to data frame
impervious <- rast(paste0("./data/impervious/processed/",
                          as.character(yoi), ".tif"))
impervious_ext <- terra::extract(impervious, coords)
data$PercentImpervious <- unlist(impervious_ext[,2]) |>
  scale()

# Extract tower distance, add to data frame
towers <- rast(paste0("./data/towers/processed/",
                      as.character(yoi), ".tif"))
towers_ext <- terra::extract(towers, coords)
data$TowerDistance <- unlist(towers_ext[,2])

# Extract tower distance, add to data frame
landfills <- rast(paste0("./data/landfills/",
                         as.character(yoi), ".tif"))
landfills_ext <- terra::extract(landfills, coords)
data$LandfillDistance <- unlist(landfills_ext[,2])

#### Conduct analysis ####

# Define basic model
model <- lm(Counts ~ NDVI + Elevation + RoadDensity + TransmissionLineDistance +
              PercentImpervious + TowerDistance + LandfillDistance, data = data)
summary(model)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.43676 -0.11903 -0.05622  0.14374  0.52890 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               6.381e-01  5.214e-01   1.224 0.260617    
# NDVI                      9.040e-03  1.732e-01   0.052 0.959829    
# Elevation                 7.294e-02  1.584e-01   0.460 0.659175    
# RoadDensity               8.827e-01  1.415e-01   6.238 0.000429 ***
#   TransmissionLineDistance  4.714e-01  3.513e-01   1.342 0.221608    
# PercentImpervious        -1.522e-01  1.459e-01  -1.043 0.331553    
# TowerDistance            -3.504e-05  2.176e-05  -1.611 0.151282    
# LandfillDistance         -6.969e-07  5.933e-06  -0.117 0.909790    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.391 on 7 degrees of freedom
# Multiple R-squared:  0.9236,	Adjusted R-squared:  0.8471 
# F-statistic: 12.08 on 7 and 7 DF,  p-value: 0.001974
vif(model)
# NDVI                Elevation              RoadDensity TransmissionLineDistance
# 2.746589                 2.298054                 1.833766                11.303909
# PercentImpervious            TowerDistance         LandfillDistance
# 1.949094                12.087701                 2.852308

# Remove towers
model <- lm(Counts ~ NDVI + Elevation + RoadDensity + TransmissionLineDistance +
              PercentImpervious + LandfillDistance, data = data)
summary(model)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.5864 -0.1961  0.0006  0.1375  0.6240 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              -8.441e-02  2.911e-01  -0.290 0.779214    
# NDVI                      1.320e-01  1.702e-01   0.775 0.460396    
# Elevation                 3.651e-02  1.717e-01   0.213 0.836934    
# RoadDensity               9.237e-01  1.524e-01   6.059 0.000303 ***
#   TransmissionLineDistance -4.004e-02  1.647e-01  -0.243 0.814065    
# PercentImpervious        -1.355e-01  1.594e-01  -0.850 0.419848    
# LandfillDistance          1.956e-06  6.241e-06   0.313 0.761952    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4282 on 8 degrees of freedom
# Multiple R-squared:  0.8952,	Adjusted R-squared:  0.8166 
# F-statistic: 11.39 on 6 and 8 DF,  p-value: 0.001518
vif(model)
# NDVI                Elevation              RoadDensity TransmissionLineDistance 
# 2.212871                 2.251215                 1.774513                 2.072192 
# PercentImpervious         LandfillDistance 
# 1.939284                 2.632395 

# Remove lines
model <- lm(Counts ~ NDVI + Elevation + RoadDensity + TowerDistance +
              PercentImpervious + LandfillDistance, data = data)
summary(model)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.60709 -0.14672  0.02069  0.11088  0.63429 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        5.847e-02  3.062e-01   0.191 0.853295    
# NDVI               7.689e-02  1.737e-01   0.443 0.669778    
# Elevation          6.211e-02  1.659e-01   0.374 0.717884    
# RoadDensity        8.901e-01  1.483e-01   6.002 0.000323 ***
#   TowerDistance     -8.665e-06  9.770e-06  -0.887 0.401008    
# PercentImpervious -1.795e-01  1.515e-01  -1.185 0.270174    
# LandfillDistance   2.129e-06  5.816e-06   0.366 0.723780    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4101 on 8 degrees of freedom
# Multiple R-squared:  0.9039,	Adjusted R-squared:  0.8318 
# F-statistic: 12.54 on 6 and 8 DF,  p-value: 0.00109
vif(model)
# NDVI         Elevation       RoadDensity     TowerDistance PercentImpervious 
# 2.512384          2.292093          1.830997          2.215873          1.911185 
# LandfillDistance 
# 2.492637 