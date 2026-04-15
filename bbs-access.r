library(tidyverse)
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
