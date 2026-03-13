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
counts = c()
for (i in 1:dim(obs_freq)[1]) {
  counts = c(
    counts,
    (bbs_obs_mojave |>
      filter(RouteNum == obs_freq$route[i]) |>
      unique() |>
      dim())[1]
  )
}
obs_freq$count = counts # This only results in 
# > sum(obs_freq$count) = 885
# TODO: Figure out why it isn't length of bbs_obs_mojave? (888)
