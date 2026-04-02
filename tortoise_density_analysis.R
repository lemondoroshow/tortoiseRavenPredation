library(dplyr)
library(ggplot2)
library(nlme)

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
  geom_point(aes(y = FK, color = "FK")) +
  geom_point(aes(y = SC, color = "SC")) +
  geom_point(aes(y = OR, color = "OR")) +
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
  geom_point(aes(y = PV, color = "PV")) +
  geom_point(aes(y = FE, color = "FE")) +
  geom_point(aes(y = CM, color = "CM")) +
  geom_point(aes(y = PT, color = "PT")) +
  geom_point(aes(y = JT, color = "JT")) +
  geom_point(aes(y = CK, color = "CK")) +
  geom_point(aes(y = AG, color = "AG")) +
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
  geom_point(aes(y = EV, color = "EV")) +
  geom_point(aes(y = IV, color = "IV")) +
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
  geom_point(aes(y = CS, color = "CS")) +
  geom_point(aes(y = MM, color = "MM")) +
  geom_point(aes(y = GB, color = "GB")) +
  geom_point(aes(y = BD, color = "BD")) +
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
