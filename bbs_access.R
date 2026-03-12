library(bbsAssistant)
library(tidyverse)
library(sf)

# Load BBS data
bbs <- grab_bbs_data()

# Load Mojave routes
mojave_routes <- read.csv("./data/bbsRoutesMojave.csv")

# Get common raven AOU
spec <- bbs$species_list |>
  filter(Scientific_Name == "Corvus corax")
aou_cc <- spec$AOU

# Filter BBS data for common ravens in the Mojave
bbs_obs_mojave <- bbs$observations |>
  mutate(RouteNum = paste0(StateNum, Route)) |>
  filter(RouteNum %in% mojave_routes$RTENO) |>
  filter(AOU == aou_cc) |>
  select(Year, RouteNum) |>
  arrange(Year)


