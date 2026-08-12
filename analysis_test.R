rm(list = ls())
library(dplyr)
library(fastDummies)
library(ggplot2)
library(nlme)
library(stringr)
library(terra)
library(tidyterra)
library(tmap)

#### Curvilinear GLS for tortoises ####

# Import data
ga_data <- read.csv("./data/tortoise_densities.csv")
years <- 1:24

# Set up variables
tcas_list <- list(
  "Western_Mojave" = c("FK", "SC", "OR"),
  "Colorado_Desert" = c("PV", "FE", "CM", "PT", "JT", "CK", "AG"),
  "Eastern_Mojave" = c("EV", "IV"),
  "Northeastern_Mojave" = c("CS", "MM", "GB", "BD")
)
rus_list <- names(tcas_list)
m_ga <- c()
b_ga <- c()

# Set colors for plots
fit_colors <- list(
  "Western_Mojave" = c("FK" = "#003f5c", "SC" = "#bc5090", "OR" = "#ffa600",
                       "tortoise trend" = "#006CD1"),
  "Colorado_Desert" = c("PV" = "#003f5c", "FE" = "#374c80", "CM" = "#7a5195", 
                        "PT" = "#bc5090", "JT" = "#ef5675", "CK" = "#ff764a", 
                        "AG" = "#ffa600", "tortoise trend" = "#006CD1"),
  "Eastern_Mojave" = c("EV" = "#003f5c", "IV" = "#bc5090", 
                       "tortoise trend" = "#006CD1"),
  "Northeastern_Mojave" = c("CS" = "#ef5675", "MM" = "#374c80", "GB" = "#7a5195", 
                            "BD" = "#bc5090", "tortoise trend" = "#006CD1") 
)

# Set coordinates for labels
label_loc <- list(
  "Western_Mojave" = list(x = 2020, y = 10),
  "Colorado_Desert" = list(x = 2020, y = 14),
  "Eastern_Mojave" = list(x = 2016, y = 4.5),
  "Northeastern_Mojave" = list(x = 2017.5, y = 6.75)
)

# Iterate through RUs
res <- data.frame(
  ru = rus_list
)
for (ru in rus_list) {
  
  # Get TCAs of interest
  tcas_tmp <- tcas_list[ru] |>
    unlist() |>
    unname()
  
  # Prepare GLS dataframe
  df_tmp <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(df_tmp) <- c("y", "time", "time.2", "stratum")
  for (tca in tcas_tmp) {
    df_tmp <- rbind(
      df_tmp,
      data.frame(
        y = ga_data[tca] |>
          unlist() |>
          unname(),
        time = years,
        time.2 = years ^ 2,
        stratum = tca
      )
    )
  }
  
  # Add dummy variables and create formula
  df_tmp <- group_by(df_tmp, time) |>
    arrange(.by_group = TRUE) |>
    dummy_cols(select_columns = "stratum", remove_first_dummy = TRUE) |>
    dplyr::select(-c("stratum"))
  fmla <- reformulate(colnames(df_tmp)[-1], response = "y") 
  
  # Fit model
  gls_fit <- gls(fmla, data = df_tmp, na.action = na.omit)
  coef <- gls_fit$coefficients
  
  # Extract fit for pop-avg
  b_ga <- coef[1]
  y <- coef[1] + coef[2] * years + coef[3] * years ^ 2
  for (i in 4:length(coef)) {
    y <- y + coef[i] * 1 / (length(coef) - 1)
    b_ga <- b_ga + coef[i] * 1 / (length(coef) - 1)
  }
  
  # Create visualization data frames
  df <- data.frame(
    year = years + 2000,
    tortoise_fit = y
  )
  
  # Plot data and fits
  res_err <- summary(gls_fit)$sigma
  x = label_loc[[ru]]$x
  y = label_loc[[ru]]$y
  colors <- fit_colors[ru] |>
    unname() |>
    unlist()
  points_list <- lapply(tcas_tmp, function(tca) {
    geom_point(data = ga_data, aes(y = .data[[tca]], color = tca))
  })
  plot <- ggplot(data = na.omit(df), aes(x = year)) + 
    points_list +
    geom_line(aes(y = tortoise_fit, color = "tortoise trend")) +
    annotate("label", x = x, y = y,
             label = paste0(
               "RSE = ", round(res_err, 2), "\n",
               "intercept (2000) = ", round(b_ga, 2)
             )
    ) +
    scale_color_manual(values = colors, name = "variable") +
    ylab("tortoises / km²") +
    ggtitle(paste0(
      "Tortoise densities and GLS fit for ",
      str_replace(ru, "_", " ")
    )) +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid = element_line(color = "grey")
    )
  ggsave(paste0("./plots/curveGLS/", ru, "_fit.png"), width = 10, height = 6,
         unit = "in")
  
}

#### Cell-level OLS for ravens ####

# Define OLS slope retrieval function
# Must account for NAs in raster
getSlope <- function(y, x) {
  if (!all(is.na(y)) & !all(is.na(x))) {
    ols_fit <- lm(y ~ x)
    return(unname(ols_fit$coefficients[2]))
  } else {
    return(NA)
  }
}

# Import all results
res_folder <- "./results/el+nd+im+ro+to+la+imXro+imXla+toXro"
files <- list.files(res_folder)
rasts <- c()
for (file in files) {
  rasts <- c(rasts, rast(paste0(res_folder, "/", file)))
}

# Apply OLS to each cell
years <- c(1:24)[-20]
res <- sds(rasts) |>
  terra::app(getSlope, years)

# Map results
tcas <- vect("data/shapefiles/TCA/USFWS_DesertTortoise_TCAs.shp") |>
  project(res)
rus <- vect("data/shapefiles/RU/2011RecoveryUnits.shp") |>
  project(res)
map <- tm_shape(res) + 
  tm_raster(
    "psi",
    col.legend = tm_legend("change in psi / year")
  ) +
  tm_basemap("Esri.WorldImagery") + 
  tm_shape(tcas) +
  tm_lines() +
  tm_shape(rus) +
  tm_lines() +
  tm_title("OLS slopes (psi / year) for raven occupancies, 2001-2024")
tmap_save(map, "./plots/cell_level_ols_cc.png")
