rm(list = ls())
library(dplyr)
library(fastDummies)
library(ggplot2)
library(nlme)
library(stringr)
library(terra)
library(tmap)

#### 0. Ravens -- preparation ####

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

#### 1. Tortoises -- final ####

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
slopes_ga <- list()
length(slopes_ga) <- 4
names(slopes_ga) <- names(tcas_list)
int_ga <- list()
length(int_ga) <- 4
names(int_ga) <- names(tcas_list)
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
  
  # Add slope to list from fit
  m_ga <- (y[length(y)] - y[1]) / (length(years) - 1)
  b_ga <- y[1] - m_ga
  slopes_ga[[ru]] <- m_ga
  int_ga[[ru]] <- b_ga
  
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
                "intercept (2000) = ", round(b_ga, 2), "\n",
                "slope = ", round(m_ga, 2)
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
  ggsave(paste0("./finalPlots/ga_", ru, "_fit.png"), width = 10, height = 6,
         unit = "in")
  
}

#### 2. Ravens -- final ####

# Import data
cc_data <- read.csv("./data/raven_psi_ru.csv")
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
fit_colors <- c("raven trend" = "#994F00", "ravens" = "red")

# Set coordinates for labels
label_loc <- list(
  "Western_Mojave" = list(x = 2012, y = 0.5),
  "Colorado_Desert" = list(x = 2017, y = 0.275),
  "Eastern_Mojave" = list(x = 2015, y = 0.375),
  "Northeastern_Mojave" = list(x = 2017.5, y = 0.32)
)

# Iterate through RUs
res <- data.frame(
  ru = rus_list
)
slopes_cc <- list()
length(slopes_cc) <- 4
names(slopes_cc) <- names(tcas_list)
int_cc <- list()
length(int_cc) <- 4
names(int_cc) <- names(tcas_list)
for (ru in rus_list) {
  
  # Fit raven model
  y_cc <- cc_data[ru] |>
    unlist() |>
    unname()
  ols_cc <- lm(y_cc ~ years)
  slopes_cc[[ru]] <- unname(ols_cc$coefficients[2])
  int_cc[[ru]] <- unname(ols_cc$coefficients[1])
  
  # Create visualization data frames
  df <- data.frame(
    year = years + 2000,
    ravens = y_cc,
    raven_fit = ols_cc$coefficients[2] * years + ols_cc$coefficients[1]
  )
  
  # Plot fit
  res_err <- summary(ols_cc)$sigma
  x = label_loc[[ru]]$x
  y = label_loc[[ru]]$y
  plot <- ggplot(data = na.omit(df), aes(x = year)) + 
    geom_point(aes(y = ravens, color = "ravens"), shape = 15, size = 3) +
    geom_line(aes(y = raven_fit, color = "raven trend")) +
    annotate("label", x = x, y = y,
             label = paste0(
               "RSE = ", round(res_err, 2), "\n",
               "intercept (2000) = ", round(ols_cc$coefficients[1], 2), "\n",
               "slope = ", round(ols_cc$coefficients[2], 4)
             )
    ) +
    scale_color_manual(values = fit_colors, name = "variable") +
    ylab("average probability of raven occupancy") +
    ggtitle(paste0(
      "Raven occupancy probabilities and OLS fit for ",
      str_replace(ru, "_", " ")
    )) +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid = element_line(color = "grey")
    )
  ggsave(paste0("./finalPlots/cc_", ru, "_fit.png"), width = 10, height = 6,
         unit = "in")
  
}
