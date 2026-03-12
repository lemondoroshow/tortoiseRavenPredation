library(bbsAssistant)
library(tidyverse)
library(sf)

# Load BBS data
bbs <- grab_bbs_data()

# Load route data
mojave_routes <- read.csv("./data/bbsRoutesMojave.csv")
bbs_rte <- bbs$routes |>
  mutate(RTENO_MJV = paste0(StateNum, Route)) |>
  filter(RTENO_MJV %in% mojave_routes$RTENO)

# Filter BBS data for routes in Mojave
bbs_obs_mojave <- bbs$observations |>
  filter(Route %in% bbs_rte$Route)
