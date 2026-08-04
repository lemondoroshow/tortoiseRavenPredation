rm(list = ls())
library(dplyr)
library(fastDummies)
library(ggplot2)
library(nlme)
library(stringr)
library(terra)
library(tmap)

#### Ravens -- preparation ####

# Load shapefiles
tcas <- vect("./data/shapefiles/TCA/USFWS_DesertTortoise_TCAs.shp") |>
  project("WGS84")
rus <- vect("./data/shapefiles/RU/2011RecoveryUnits.shp") |>
  project("WGS84")

# Create data frame
raven_occ_tca <- data.frame(matrix(2001:2024, nrow = 24, ncol = 1))
colnames(raven_occ_tca) <- "year"
raven_occ_ru <- data.frame(raven_occ_tca)

# Iterate through TCAs
model_name <- "el+nd+im+ro+to+la+imXro+imXla+toXro"
for (i in 1:dim(tcas)[1]) {
  
  # Subset TCA vector
  tca <- subset(tcas, i)
  name <- tca$stratum
  psis <- c()
  
  # Iterate through years
  for (yoi in 2001:2024) {
    if (yoi != 2020) {
      
      # Load year's results raster; project, mask, and crop
      f <- list.files(paste0("./results/", model_name), 
                      pattern = as.character(yoi), full.names = TRUE)[1]
      res <- rast(f) |>
        project(tcas) |>
        mask(tca) |>
        crop(tca)
      
      # Calculate average probability of occurence
      psi_avg <- global(res, fun = "mean", na.rm = TRUE)
      psis <- c(psis, psi_avg)
    } else {
      psis <- c(psis, NA)
    }
  }
  
  # Add TCA to data frame
  raven_occ_tca <- mutate(raven_occ_tca, !!name := unlist(psis))
}

# Iterate through RUs
for (i in 1:dim(rus)[1]) {
  
  # Subset TCA vector
  ru <- subset(rus, i)
  name <- gsub(" ", "", ru$Unit_Name)
  psis <- c()
  
  # Iterate through years
  for (yoi in 2001:2024) {
    if (yoi != 2020) {
      
      # Load year's results raster; project, mask, and crop
      f <- list.files(paste0("./results/", model_name), 
                      pattern = as.character(yoi), full.names = TRUE)[1]
      res <- rast(f) |>
        project(rus) |>
        mask(ru) |>
        crop(ru)
      
      # Calculate average probability of occurence
      psi_avg <- global(res, fun = "mean", na.rm = TRUE)
      psis <- c(psis, psi_avg)
    } else {
      psis <- c(psis, NA)
    }
  }
  
  # Add TCA to data frame
  raven_occ_ru <- mutate(raven_occ_ru, !!name := unlist(psis))
}

# Export data frames
write.csv(raven_occ_tca, "./data/raven_psi_tca.csv", row.names = FALSE)
write.csv(raven_occ_ru, "./data/raven_psi_ru.csv", row.names = FALSE)

#### Ravens -- visualization ####

# Import data
tca_data_ga <- read.csv("./data/tortoise_densities.csv")
tca_data_cc <- read.csv("./data/raven_psi_tca.csv")
ru_data_cc <- read.csv("./data/raven_psi_ru.csv")

# Iterate through TCAs
years <- 2001:2024
for (tca in colnames(tca_data_ga)[-1]) {
  
  # Isolate data for TCA
  tortoise <- dplyr::select(tca_data_ga, !!tca) |>
    dplyr::rename(tortoise = !!tca)
  raven <- dplyr::select(tca_data_cc, !!tca) |>
    dplyr::rename(raven = !!tca)
  data <- data.frame(
    year = years, 
    tortoise = scale(tortoise), 
    raven = scale(raven)
  )
  
  # Plot data
  colors <- c("G. agassizii\ndensity" = "blue", "C. corax\nocc. prob." = "red")
  plot <- ggplot(na.omit(data), aes(x = year)) +
    geom_point(aes(y = tortoise, color = "G. agassizii\ndensity")) +
    geom_line(aes(y = tortoise), color = "blue") +
    geom_point(aes(y = raven, color = "C. corax\nocc. prob.")) +
    geom_line(aes(y = raven), color = "red") +
    scale_color_manual(values = colors, name = "variable") +
    ylab("scaled presence variable") +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid.minor = element_line(color = "grey"),
      panel.grid.major = element_line(color = "grey")
    ) +
    ggtitle(paste0("Tortoise densities and raven occurrence probabilities in ", 
                   tca))
  
  # Save plot
  ggsave(paste0("./plots/tcaTortoiseRaven/", tca, ".png"), plot)
}

# Iterate through RUs
for (ru in colnames(ru_data_cc)[-1]) {
  
  # Isolate data for RU
  data <- dplyr::select(ru_data_cc, c(year, !!ru)) |>
    rename(raven = !!ru)
  
  # Plot data
  plot <- ggplot(na.omit(data), aes(x = year, y = raven)) +
    geom_point() + 
    geom_line() +
    ylab("mean raven occ. prob.") + 
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid = element_line(color = "grey")
    ) +
    ggtitle(paste0("Raven occ. prob. in ", ru))
  
  # Save plot
  ggsave(paste0("./plots/ruRaven/", ru, ".png"))
  
}

#### Tortoises -- raw visualization ####

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

# Set colors for plots
fit_colors <- list(
  "Western_Mojave" = c("FK" = "#003f5c", "SC" = "#bc5090", "OR" = "#ffa600"),
  "Colorado_Desert" = c("PV" = "#003f5c", "FE" = "#374c80", "CM" = "#7a5195", 
                        "PT" = "#bc5090", "JT" = "ef5675#", "CK" = "#ff764a", 
                        "AG" = "#ffa600"),
  "Eastern_Mojave" = c("EV" = "#003f5c", "IV" = "#bc5090"),
  "Northeastern_Mojave" = c("CS" = "#003f5c", "MM" = "#374c80", "GB" = "#7a5195", 
                            "BD" = "#bc5090") 
)

# Iterate through RUs
for (ru in rus_list) {
  
  # Get TCAs of interest
  tcas_tmp <- tcas_list[ru] |>
    unlist() |>
    unname()
  
  # Plot data and fits
  colors <- fit_colors[ru] |>
    unname() |>
    unlist()
  points_list <- lapply(tcas_tmp, function(tca) {
    geom_point(data = ga_data, aes(y = .data[[tca]], color = tca))
  })
  lines_list <- lapply(tcas_tmp, function(tca) {
    geom_line(data = dplyr::select(ga_data, c("year", !!tca)) |>
                na.omit(), aes(y = .data[[tca]], color = tca), na.rm = T)
  })
  plot <- ggplot(data = data.frame(years = years + 2000), aes(x = year)) + 
    points_list +
    lines_list + 
    scale_color_manual(values = colors, name = "TCA") +
    ylab("tortoises / km²") +
    ggtitle(paste0(
      "Annual tortoise densities per USFWS for ",
      str_replace(ru, "_", " ")
    )) +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid = element_line(color = "grey")
    )
  ggsave(paste0("./finalPlots/", ru, "_tortoise_densities.png"))
  
}

#### Tortoises and ravens -- RU level ####

# Import data
ga_data <- read.csv("./data/tortoise_densities.csv")
cc_data <- read.csv("./data/raven_psi_ru.csv")
years <- 1:24

# Set up variables
tcas_list <- list(
  "WesternMojave" = c("FK", "SC", "OR"),
  "ColoradoDesert" = c("PV", "FE", "CM", "PT", "JT", "CK", "AG"),
  "EasternMojave" = c("EV", "IV"),
  "NortheasternMojave" = c("CS", "MM", "GB", "BD")
)
rus_list <- names(tcas_list)
m_ga <- c()
b_ga <- c()
m_cc <- c()
b_cc <- c()

# Set colors for plots
fit_colors <- list(
  "WesternMojave" = c("FK" = "#003f5c", "SC" = "#bc5090", "OR" = "#ffa600",
                      "tortoise trend" = "#006CD1", "raven trend" = "#994F00", 
                      "ravens" = "red"),
  "ColoradoDesert" = c("PV" = "#003f5c", "FE" = "#374c80", "CM" = "#7a5195", 
                       "PT" = "#bc5090", "JT" = "#ef5675", "CK" = "#ff764a", 
                       "AG" = "#ffa600", "tortoise trend" = "#006CD1", 
                       "raven trend" = "#994F00", "ravens" = "red"),
  "EasternMojave" = c("EV" = "#003f5c", "IV" = "#bc5090", 
                      "tortoise trend" = "#006CD1", "raven trend" = "#994F00", 
                      "ravens" = "red"),
  "NortheasternMojave" = c("CS" = "#003f5c", "MM" = "#374c80", "GB" = "#7a5195", 
                           "BD" = "#bc5090", "tortoise trend" = "#006CD1",
                           "raven trend" = "#994F00", "ravens" = "red") 
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
          unname() |>
          scale(),
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
  # Remove [-2] to include curvilinear time relationship
  fmla <- reformulate(colnames(df_tmp)[-1][-2], response = "y") 
  
  # Fit model
  gls_fit <- gls(fmla, data = df_tmp, na.action = na.omit)
  coef <- gls_fit$coefficients
  b_ga <- c(b_ga, unname(coef[1]))
  
  # Extract fit for pop-avg
  y <- coef[1] + coef[2] * years
  for (i in 3:length(coef)) {
    y <- y + coef[i] * 1 / (length(coef) - 1)
  }
  
  # Add slope to list from fit
  m_ga <- c(
    m_ga,
    (y[length(y)] - y[1]) / (length(years) - 1)
  )
  
  # Fit raven model
  y_cc <- cc_data[ru] |>
    unlist() |>
    unname() |>
    scale()
  ols_cc <- lm(y_cc ~ years)
  m_cc <- c(m_cc, unname(ols_cc$coefficients[2]))
  b_cc <- c(b_cc, unname(ols_cc$coefficients[1]))
  
  # Create visualization data frames
  df <- data.frame(
    year = years + 2000,
    ravens = y_cc,
    raven_fit = ols_cc$coefficients[2] * years + ols_cc$coefficients[1],
    tortoise_fit = y
  )

  # Plot data and fits
  colors <- fit_colors[ru] |>
    unname() |>
    unlist()
  points_list <- lapply(tcas_tmp, function(tca) {
    geom_point(data = ga_data, aes(y = scale(.data[[tca]]), color = tca))
  })
  plot <- ggplot(data = na.omit(df), aes(x = year)) + 
    points_list +
    geom_point(aes(y = ravens, color = "ravens"), shape = 15, size = 3) +
    geom_line(aes(y = raven_fit, color = "raven trend")) +
    geom_line(aes(y = tortoise_fit, color = "tortoise trend")) +
    scale_color_manual(values = colors, name = "variable") +
    ylab("scaled presence variable") +
    ggtitle(paste0(
      "Average tortoise density and raven occupancy probability for ",
      ru
    )) +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid = element_line(color = "grey")
    )
  ggsave(paste0("./plots/", ru, "_fit.png"))
  
}

# Add data to results
res <- dplyr::mutate(
  res,
  SlopeTortoise = m_ga,
  InterceptTortoise = b_ga,
  SlopeRaven = m_cc,
  InterceptRaven = b_cc,
  TxR = m_ga * m_cc
)
write.csv(res, "./results/fit_tortoise_raven.csv", row.names = FALSE)

#### Tortoises and ravens -- TCA level ####

# Import data
tcas <- vect("./data/shapefiles/TCA/USFWS_DesertTortoise_TCAs.shp") |>
  project("WGS84")
ga_data <- read.csv("./data/tortoise_densities.csv")
cc_data <- read.csv("./data/raven_psi_tca.csv")
years <- 1:24

# Iterate through TCAs
res <- c()
for (tca in colnames(ga_data[-1])) {
  
  # Fit tortoise model
  y_ga <- ga_data[tca] |>
    unlist() |>
    unname()
  ols_ga <- lm(y_ga ~ years)
  m_ga <- unname(ols_ga$coefficients[2])
  
  # Fit raven model
  y_cc <- cc_data[tca] |>
    unlist() |>
    unname()
  ols_cc <- lm(y_cc ~ years)
  m_cc <- unname(ols_cc$coefficients[2])
  
  # Add metrics to results
  res <- c(res, m_ga * m_cc) |>
    unlist()
}

# Format, scale data
df <- data.frame(matrix(nrow = length(colnames(ga_data)[-1])))
df$stratum <- colnames(ga_data)[-1]
df$res <- scale(res)
df <- dplyr::select(df, c(stratum, res))
tcas_new <- left_join(tcas, df, by = c("stratum"))
tcas_new$res <- as.numeric(tcas_new$res)

# Visualize
map <- tm_shape(tcas_new) +
  tm_polygons(
    fill = "res", 
    fill.scale = tm_scale(values = "-tableau.classic_blue"),
    fill.legend = tm_legend(title = "TxR slopes")) +
  tm_basemap("Esri.WorldImagery") +
  tm_title(
    "Scaled tortoise times raven slopes for OLS fits of annual presence",
    size = 1
  )
tmap_save(map, "./plots/txr_tca.png")

#### Tortoises -- RU level, GLS ####

# Import tortoise data
df <- read.csv("./data/tortoise_densities.csv")
years <- 1:24

# Set up variables
tcas_list <- list(
  "WesternMojave" = c("FK", "SC", "OR"),
  "ColoradoDesert" = c("PV", "FE", "CM", "PT", "JT", "CK", "AG"),
  "EasternMojave" = c("EV", "IV"),
  "NortheasternMojave" = c("CS", "MM", "GB", "BD")
)
rus_list <- names(tcas_list)

# Set colors for plots
fit_colors <- list(
  "WesternMojave" = c("FK" = "#003f5c", "SC" = "#bc5090", "OR" = "#ffa600",
                      "GLS" = "#006CD1"),
  "ColoradoDesert" = c("PV" = "#003f5c", "FE" = "#374c80", "CM" = "#7a5195", 
                       "PT" = "#bc5090", "JT" = "#ef5675", "CK" = "#ff764a", 
                       "AG" = "#ffa600", "GLS" = "#006CD1"),
  "EasternMojave" = c("EV" = "#003f5c", "IV" = "#bc5090", "GLS" = "#006CD1"),
  "NortheasternMojave" = c("CS" = "#003f5c", "MM" = "#374c80", "GB" = "#7a5195", 
                           "BD" = "#bc5090", "GLS" = "#006CD1") 
)

# Iterate through RUs
res <- data.frame(year = years)
for (ru in rus_list) {
  
  # Prepare GLS dataframe
  tcas_tmp <- tcas_list[ru] |>
    unlist() |>
    unname()
  df_tmp <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(df_tmp) <- c("y", "time", "time.2", "stratum")
  for (tca in tcas_tmp) {
    df_tmp <- rbind(
      df_tmp,
      data.frame(
        y = df[tca] |>
          unlist() |>
          unname() #|>
          #scale(), # Uncomment these lines when calculating results
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
  # Remove [-2] to include curvilinear time relationship
  fmla <- reformulate(colnames(df_tmp)[-1][-2], response = "y") 
  
  # Fit model
  gls_fit <- gls(fmla, data = df_tmp, na.action = na.omit)
  coef <- gls_fit$coefficients
  
  # Extract fit for pop-avg
  y <- coef[1] + coef[2] * years
  for (i in 3:length(coef)) {
    y <- y + coef[i] * 1 / (length(coef) - 1)
  }
  gls_res <- data.frame(
    year = 2001:2024,
    y = y
  )
  res[ru] <- gls_res$y
  
  # Plot results
  colors <- fit_colors[ru] |>
    unname() |>
    unlist()
  plot <- ggplot(df, aes(x = year))
  points_list <- lapply(tcas_tmp, function(tca) {
    geom_point(data = df, aes(y = .data[[tca]], color = tca))
  })
  plot <- plot +
    points_list +
    geom_line(aes(y = y, color = "GLS"), data = gls_res) +
    scale_color_manual(values = colors, name  = "TCA") +
    theme(
      panel.background = element_rect(fill = "white"), 
      panel.grid = element_line(color = "grey")
    ) +
    xlab("Year") +
    ylab("Density (tortoises / km-sqared)") + 
    ggtitle(paste0(ru, " tortoise densities, GLS fit"))
  ggsave(paste0("./plots/", ru, "_GLS_tortoise.png"))
}

# Export fits # Do not export without uncommenting lines 397 and 398
# May need to be fixed as of 07/07/2026
# write.csv(res, "./results/gls_tortoise_fits.csv")
