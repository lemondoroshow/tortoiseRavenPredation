rm(list = ls())
library(dplyr)
library(ggplot2)
library(nlme)
library(terra)

#### Tortoises ####

# Import tortoise data
df <- read.csv("./data/tortoise_densities.csv")
years <- df$year - 2000

# Log-transform densities
# df[2:dim(df)[2]] <- lapply(df[2:dim(df)[2]], log)

# Define RUs and TCAs, for my own reference
western_mojave <- c("FK", "SC", "OR")
colorado_desert <- c("PV", "FE", "CM", "PT", "JT", "CK", "AG")
eastern_mojave <- c("EV", "IV")
northeastern_mojave <- c("CS", "MM", "GB", "BD")

## Western Mojave

# Average densities
y_wm <- rowMeans(data.frame(df$FK, df$SC, df$OR), na.rm = TRUE)
ols_wm <- lm(y_wm ~ df$year)

# Isolate OLS fit
b_wm <- ols_wm$coefficients[1]
m_wm <- ols_wm$coefficients[2]
ols_wm_res <- data.frame(year = df$year, y = df$year * m_wm + b_wm)

# Create GLS dataframe
df_wm <- data.frame(
  y = df$FK,
  time = df$year - 2000,
  stratum = 1
) |> rbind(data.frame(
  y = df$SC,
  time = df$year - 2000,
  stratum = 2
)) |> rbind(data.frame(
  y = df$OR,
  time = df$year - 2000,
  stratum = 3
)) |> 
  group_by(time) |>
  arrange(.by_group = TRUE)

# Fit model
gls_wm <- gls(y ~ time + stratum + time * stratum,
              data = df_wm, na.action = na.omit)

# Extract OLS fit
gls_wm_res <- data.frame(
  year = df$year,
  y = gls_wm$coefficients[1] +
      gls_wm$coefficients[2] * years +
      gls_wm$coefficients[3] * (1 + 2 + 3) / 3 + # Population-averaged stratum effect
      gls_wm$coefficients[4] * (1 + 2 + 3) / 3 * years # Population-averaged stratum x time effect
)

# Plot fits
wm_fit_colors <- c("FK" = "#003f5c", "SC" = "#bc5090", "OR" = "#ffa600", ".OLS Fit" = "#994F00",
                   ".GLS Fit" = "#006CD1")
wm_fit_plot <- ggplot(df, aes(x = year)) +
  geom_point(aes(y = FK, color = "FK"), shape = 16) +
  geom_point(aes(y = SC, color = "SC"), shape = 17) +
  geom_point(aes(y = OR, color = "OR"), shape = 18) +
  geom_line(aes(y = y, color = ".OLS Fit"), data = ols_wm_res) +
  geom_line(aes(y = y, color = ".GLS Fit"), data = gls_wm_res) +
  scale_color_manual(values = wm_fit_colors, name  = "TCA") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey")) +
  xlab("Year") +
  ylab("Density (tortoises / km-sqared)") + 
  ggtitle("Western Mojave tortoise densities")
ggsave("./plots/tortoiseDensities/western_mojave_fit_nolog.png", wm_fit_plot)

## Colorado Desert

# Average densities
y_cd <- rowMeans(data.frame(df$PV, df$FE, df$CM, df$PT, df$JT, df$CK, df$AG), na.rm = TRUE)
ols_cd <- lm(y_cd ~ df$year)

# Isolate OLS fit
b_cd <- ols_cd$coefficients[1]
m_cd <- ols_cd$coefficients[2]
ols_cd_res <- data.frame(year = df$year, y = df$year * m_cd + b_cd)

# Create GLS dataframe
df_cd <- data.frame(
  y = df$PV,
  time = df$year - 2000,
  stratum = 1
) |> rbind(data.frame(
  y = df$FE,
  time = df$year - 2000,
  stratum = 2
)) |> rbind(data.frame(
  y = df$CM,
  time = df$year - 2000,
  stratum = 3
)) |> rbind(data.frame(
  y = df$PT,
  time = df$year - 2000,
  stratum = 4
)) |> rbind(data.frame(
  y = df$JT,
  time = df$year - 2000,
  stratum = 5
)) |> rbind(data.frame(
  y = df$CK,
  time = df$year - 2000,
  stratum = 6
)) |> rbind(data.frame(
  y = df$AG,
  time = df$year - 2000,
  stratum = 7
)) |> 
  group_by(time) |>
  arrange(.by_group = TRUE)

# Fit model
gls_cd <- gls(y ~ time + stratum + time * stratum,
              data = df_cd, na.action = na.omit)

# Extract OLS fit
gls_cd_res <- data.frame(
  year = df$year,
  y = gls_cd$coefficients[1] +
      gls_cd$coefficients[2] * years +
      gls_cd$coefficients[3] * (1 + 2 + 3 + 4 + 5 + 6 + 7) / 7 + # Population-averaged stratum effect
      gls_cd$coefficients[4] * (1 + 2 + 3 + 4 + 5 + 6 + 7) / 7 * years # Population-averaged stratum x time effect
)

# Plot fits
cd_fit_colors <- c("PV" = "#003f5c", "FE" = "#374c80", "CM" = "#7a5195", "PT" = "#bc5090", 
                   "JT" = "#ef5675", "CK" = "#ff764a", "AG" = "#ffa600", ".OLS Fit" = "#994F00", 
                   ".GLS Fit" = "#006CD1")
cd_fit_plot <- ggplot(df, aes(x = year)) +
  geom_point(aes(y = PV, color = "PV"), shape = 16) +
  geom_point(aes(y = FE, color = "FE"), shape = 17) +
  geom_point(aes(y = CM, color = "CM"), shape = 18) +
  geom_point(aes(y = PT, color = "PT"), shape = 4) +
  geom_point(aes(y = JT, color = "JT"), shape = 3) +
  geom_point(aes(y = CK, color = "CK"), shape = 15) +
  geom_point(aes(y = AG, color = "AG"), shape = 8) +
  geom_line(aes(y = y, color = ".OLS Fit"), data = ols_cd_res) +
  geom_line(aes(y = y, color = ".GLS Fit"), data = gls_cd_res) +
  scale_color_manual(values = cd_fit_colors, name  = "TCA") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey")) +
  xlab("Year") +
  ylab("Density (tortoises / km-sqared)") + 
  ggtitle("Colorado Desert tortoise densities")
ggsave("./plots/tortoiseDensities/colorado_desert_fit_nolog.png", cd_fit_plot)

## Eastern Mojave

# Average densities
y_em <- rowMeans(data.frame(df$EV, df$IV), na.rm = TRUE)
ols_em <- lm(y_em ~ df$year)

# Isolate OLS fit
b_em <- ols_em$coefficients[1]
m_em <- ols_em$coefficients[2]
ols_em_res <- data.frame(year = df$year, y = df$year * m_em + b_em)

# Create GLS dataframe
df_em <- data.frame(
  y = df$EV,
  time = df$year - 2000,
  stratum = 1
) |> rbind(data.frame(
  y = df$IV,
  time = df$year - 2000,
  stratum = 2
)) |> 
  group_by(time) |>
  arrange(.by_group = TRUE)

# Fit model
gls_em <- gls(y ~ time + stratum + time * stratum,
              data = df_em, na.action = na.omit)

# Extract OLS fit
gls_em_res <- data.frame(
  year = df$year,
  y = gls_em$coefficients[1] +
      gls_em$coefficients[2] * years +
      gls_em$coefficients[3] * (1 + 2) / 2 + # Population-averaged stratum effect
      gls_em$coefficients[4] * (1 + 2) / 2 * years # Population-averaged stratum x time effect
)

# Plot fits
em_fit_colors <- c("EV" = "#003f5c", "IV" = "#bc5090", ".OLS Fit" = "#994F00", ".GLS Fit" = "#006CD1")
em_fit_plot <- ggplot(df, aes(x = year)) +
  geom_point(aes(y = EV, color = "EV"), shape = 16) +
  geom_point(aes(y = IV, color = "IV"), shape = 17) +
  geom_line(aes(y = y, color = ".OLS Fit"), data = ols_em_res) +
  geom_line(aes(y = y, color = ".GLS Fit"), data = gls_em_res) +
  scale_color_manual(values = em_fit_colors, name  = "TCA") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey")) +
  xlab("Year") +
  ylab("Density (tortoises / km-sqared)") + 
  ggtitle("Eastern Mojave tortoise densities")
ggsave("./plots/tortoiseDensities/eastern_mojave_fit_nolog.png", em_fit_plot)

## Northeastern Mojave

# Average densities
y_nm <- rowMeans(data.frame(df$CS, df$MM, df$GB, df$BD), na.rm = TRUE)
ols_nm <- lm(y_nm ~ df$year)

# Isolate OLS fit
b_nm <- ols_nm$coefficients[1]
m_nm <- ols_nm$coefficients[2]
ols_nm_res <- data.frame(year = df$year, y = df$year * m_nm + b_nm)

# Create GLS dataframe
df_nm <- data.frame(
  y = df$CS,
  time = df$year - 2000,
  stratum = 1
) |> rbind(data.frame(
  y = df$MM,
  time = df$year - 2000,
  stratum = 2
)) |> rbind(data.frame(
  y = df$GB,
  time = df$year - 2000,
  stratum = 3
)) |> rbind(data.frame(
  y = df$BD,
  time = df$year - 2000,
  stratum = 4
)) |> 
  group_by(time) |>
  arrange(.by_group = TRUE)

# Fit model
gls_nm <- gls(y ~ time + stratum + time * stratum,
              data = df_nm, na.action = na.omit)

# Extract OLS fit
gls_nm_res <- data.frame(
  year = df$year,
  y = gls_nm$coefficients[1] +
      gls_nm$coefficients[2] * years +
      gls_nm$coefficients[3] * (1 + 2 + 3 + 4) / 4 + # Population-averaged stratum effect
      gls_nm$coefficients[4] * (1 + 2 + 3 + 4) / 4 * years # Population-averaged stratum x time effect
)

# Plot fits
nm_fit_colors <- c("CS" = "#003f5c", "MM" = "#374c80", "GB" = "#7a5195", "BD" = "#bc5090", 
                   ".OLS Fit" = "#994F00", ".GLS Fit" = "#006CD1")
nm_fit_plot <- ggplot(df, aes(x = year)) +
  geom_point(aes(y = CS, color = "CS"), shape = 16) +
  geom_point(aes(y = MM, color = "MM"), shape = 17) +
  geom_point(aes(y = GB, color = "GB"), shape = 18) +
  geom_point(aes(y = BD, color = "BD"), shape = 4) +
  geom_line(aes(y = y, color = ".OLS Fit"), data = ols_nm_res) +
  geom_line(aes(y = y, color = ".GLS Fit"), data = gls_nm_res) +
  scale_color_manual(values = nm_fit_colors, name  = "TCA") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey")) +
  xlab("Year") +
  ylab("Density (tortoises / km-sqared)") + 
  ggtitle("Northeastern Mojave tortoise densities")
ggsave("./plots/tortoiseDensities/northeastern_mojave_fit_nolog.png", nm_fit_plot)

# Clean up
rm(list=setdiff(ls(), "df"))

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
