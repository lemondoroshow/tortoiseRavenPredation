library(tidyverse)
library(tmap)
library(sf)

# Load Mojave routes
mojave_routes <- read.csv("./data/bbs_routes_mojave.csv")

# Create x- and y-lists
route_xs <- mojave_routes$x
route_ys <- mojave_routes$y
names(route_xs) <- mojave_routes$RTENO
names(route_ys) <- mojave_routes$RTENO

# Get common raven AOU
spec <- read.csv("./data/bbs/SpeciesList.csv") |>
  filter(English_Common_Name == "Common Raven")
aou_cc <- spec$AOU

# Import BBS data and filter for common ravens in the Mojave
bbs_obs_mojave <- read.csv("./data/bbs/States/Arizona.csv") |>
  rbind(read.csv("./data/bbs/States/Califor.csv")) |>
  rbind(read.csv("./data/bbs/States/Nevada.csv")) |>
  rbind(read.csv("./data/bbs/States/Utah.csv")) |>
  mutate(RouteNum = paste0(StateNum, formatC(Route, 2, flag = "0"))) |>
  filter(RouteNum %in% mojave_routes$RTENO) |>
  filter(AOU == aou_cc) |>
  select(Year, RouteNum) |>
  arrange(Year) |>
  mutate(x = route_xs[RouteNum], y = route_ys[RouteNum])

# Calculate frequencies per route, total span
obs_freq <- data.frame(route = unique(bbs_obs_mojave$RouteNum))
counts <- c()
xs <- c()
ys <- c()
names <- c()
for (i in 1:dim(obs_freq)[1]) {
  counts <- c(
    counts,
    (bbs_obs_mojave |>
      filter(RouteNum == obs_freq$route[i]) |>
      dim())[1]
  )

  # NOTE -- Barstow and Inyokern are repeated in name, but not in coordinates
  # This WILL need to be fixed / figured out later, but for this high-level analysis
  # I'm not as worried
  xs <- c(
    xs,
    subset(mojave_routes, RTENO == obs_freq$route[i])$x[1]
  )
  ys <- c(
    ys, 
    subset(mojave_routes, RTENO == obs_freq$route[i])$y[1]
  )
  names <- c(
    names,
    subset(mojave_routes, RTENO == obs_freq$route[i])$RTENAME[1]
  )
}

# Add data to data frame
obs_freq$count <- counts
obs_freq$x <- xs
obs_freq$y <- ys
obs_freq$name <- names

# Create spatial data
obs_freq_sf <- st_as_sf(obs_freq, coords = c("x", "y"))
st_crs(obs_freq_sf) <- "WGS84"

# Create bounding box
bbox <- st_bbox(obs_freq_sf)
x_range <- bbox$xmax - bbox$xmin
y_range <- bbox$ymax - bbox$ymin
bbox[1] <- bbox[1] - 0.25 * x_range
bbox[3] <- bbox[3] + 0.25 * x_range
bbox[2] <- bbox[2] - 0.25 * y_range
bbox[4] <- bbox[4] + 0.25 * y_range
bbox <- st_as_sfc(bbox)
st_crs(bbox) <- "WGS84"

# Map data
#tm_shape(World, bbox = st_bbox(obs_freq_sf) * 1.1) +
corvax_map<- tm_shape(obs_freq_sf, bbox = bbox) +
  tm_symbols(size = "count", size.scale = tm_scale_continuous(values.range=c(0.1, 3)), fill = "purple") +
  tm_basemap("OpenStreetMap") +
  tm_title("C. corax counts for the Mojave, 1968 - Present")
tmap_save(corvax_map, "./plots/corvax_count_map.png")