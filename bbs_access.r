library(tidyverse)

# Load Mojave routes
mojave_routes <- read.csv("./data/bbs_routes_mojave.csv")

# Get common raven AOU
spec <- read.csv("./data/bbs/SpeciesList.csv") %>%
  filter(English_Common_Name == "Common Raven")
aou_cc <- spec$AOU

# Import weather data
weather <- read.csv("./data/bbs/Weather.csv") |>
  unite('date', sep = '-', Year, Month, Day, remove = FALSE) %>%
  mutate(date = as.Date(date, tz = "America/New_York")) %>%
  mutate(julian = as.numeric(format(date, '%j')))

# Import BBS data, and filter
bbs_obs_mojave <- read.csv("./data/bbs/States/Arizona.csv") %>%
  bind_rows(read.csv("./data/bbs/States/Califor.csv")) %>%
  bind_rows(read.csv("./data/bbs/States/Nevada.csv")) %>%
  bind_rows(read.csv("./data/bbs/States/Utah.csv")) %>%
  left_join(read.csv("./data/bbs/Routes.csv"), 
            by = c('Route', 'CountryNum', 'StateNum')) %>%
  filter(Year == 2024) %>%
  mutate(RouteNum = paste0(StateNum, formatC(Route, 2, flag = "0"))) %>%
  filter(RouteNum %in% mojave_routes$RTENO) %>%
  dplyr::select(RouteDataID, Latitude, Longitude, AOU, starts_with("Count")) %>%
  dplyr::select(-CountryNum) %>%
  complete(AOU, nesting(RouteDataID, Latitude, Longitude)) %>%
  replace(is.na(.), 0) %>%
  filter(AOU == aou_cc) %>%
  left_join(weather, by = c("RouteDataID"))

# %>%
#   mutate(RouteNum = paste0(StateNum, formatC(Route, 2, flag = "0"))) %>%
#   filter(RouteNum %in% mojave_routes$RTENO) %>%
#   filter(AOU == aou_cc) %>%
#   filter(Year >= 2001) %>%
#   select(Year, RouteNum, SpeciesTotal) %>%
#   arrange(Year) %>%
#   mutate(x = route_xs[RouteNum], y = route_ys[RouteNum])