rm(list = ls())
library(tidyverse)
library(terra)
library(tidyterra)

# Import BBS data
bbs_obs <- read.csv("./data/bbs/States/Arizona.csv") |>
  bind_rows(read.csv("./data/bbs/States/Califor.csv")) |>
  bind_rows(read.csv("./data/bbs/States/Nevada.csv")) |>
  bind_rows(read.csv("./data/bbs/States/Utah.csv")) |>
  left_join(read.csv("./data/bbs/Routes.csv"),
            by = c("Route", "CountryNum", "StateNum")) |>
  mutate(RouteNum = paste0(StateNum, formatC(Route, 2, flag = "0"))) |>
  dplyr::filter(Year >= 2001) |>
  dplyr::select(RouteDataID, Latitude, Longitude, AOU, starts_with("Count")) |>
  dplyr::select(-CountryNum) |>
  complete(AOU, nesting(RouteDataID, Latitude, Longitude)) %>%
  replace(is.na(.), 0)

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

# Get total counts
aou_counts_mojave = mutate(bbs_obs_mojave,
  Total.Counts = dplyr::select(
    bbs_obs_mojave, 
    c(Count10, Count20, Count30, Count40, Count50)
  ) |>
    as.matrix() |>
    rowSums()
) |>
  dplyr::select(c("AOU", "Total.Counts")) |>
  group_by(AOU) |>
  summarise(Total.Count = sum(Total.Counts)) |>
  arrange(.by = desc(Total.Count))

# Add species names to data frame
species <- read.csv("./data/bbs/SpeciesList.csv")
species_from_aou <- function(aou, species_data, name = "common") {
  df <- dplyr::filter(species, AOU == aou)
  if (name == "common") {
    return(df$English_Common_Name)
  } else if (name == "binomial") {
    return(paste0(df$Genus, " ", df$Species))
  } else if (name == "order") {
    return(df$Order)
  } else if (name == "family") {
    return(df$Family)
  } else {
    return(NA)
  }
}
species_counts_mojave <- mutate(
  aou_counts_mojave,
  Common.Name = lapply(aou_counts_mojave$AOU, species_from_aou,
                       species, "common") |>
    unlist(),
  Binomial.Name = lapply(aou_counts_mojave$AOU, species_from_aou,
                         species, "binomial") |>
    unlist(),
  Order = lapply(aou_counts_mojave$AOU, species_from_aou,
                 species, "order") |>
    unlist(),
  Family = lapply(aou_counts_mojave$AOU, species_from_aou,
                  species, "family") |>
    unlist(),
)
write.csv(species_counts_mojave, "./data/species_counts_mojave.csv")
