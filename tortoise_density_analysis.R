library(dplyr)
library(ggplot2)

# Import tortoise data
df <- read.csv("./data/tortoise_densities.csv")

# Log-transform densities
df[2:dim(df)[2]] <- lapply(df[2:dim(df)[2]], log)

# Define RUs and TCAs, for my own reference
western_mojave <- c("FK", "SC", "OR")
colorado_desert <- c("PV", "FE", "CM", "PT", "JT", "CK", "AG")
eastern_mojave <- c("EV", "IV")
northeastern_mojave <- c("CS", "MM", "GB", "BD")

# Plot WM log-densities
wm_colors <- c("FK" = "#003f5c", "SC" = "#bc5090", "OR" = "#ffa600")
wm_plot <- ggplot(df, aes(x = year)) +
  geom_point(aes(y = FK, color = "FK")) +
  geom_point(aes(y = SC, color = "SC")) +
  geom_point(aes(y = OR, color = "OR")) +
  scale_color_manual(values = wm_colors, name  = "TCA") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey")) +
  xlab("Year") +
  ylab("Log-density (tortoises / km-sqared)") + 
  ggtitle("Western Mojave tortoise log-densities")
ggsave("./plots/tortoiseDensities/western_mojave.png", wm_plot)

# Average WM densities
y_wm <- (df$FK + df$SC + df$OR) * 1 / 3
fit_wm <- lm(y_wm ~ df$year)

# Isolate fit
b_wm <- fit_wm$coefficients[1]
m_wm <- fit_wm$coefficients[2]
line_wm <- data.frame(year = df$year, y = df$year * m_wm + b_wm)

# Plot fit
wm_fit_colors <- c("FK" = "#003f5c", "SC" = "#bc5090", "OR" = "#ffa600", "Fit" = "black")
wm_fit_plot <- ggplot(df, aes(x = year)) +
  geom_point(aes(y = FK, color = "FK")) +
  geom_point(aes(y = SC, color = "SC")) +
  geom_point(aes(y = OR, color = "OR")) +
  geom_line(aes(y = y, color = "Fit"), data = line_wm) +
  scale_color_manual(values = wm_fit_colors, name  = "TCA") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey")) +
  xlab("Year") +
  ylab("Log-density (tortoises / km-sqared)") + 
  ggtitle("Western Mojave tortoise log-densities")
ggsave("./plots/tortoiseDensities/western_mojave_fit.png", wm_fit_plot)

# Plot CD log-densities
cd_colors <- c("PV" = "#003f5c", "FE" = "#374c80", "CM" = "#7a5195", "PT" = "#bc5090", 
              "JT" = "#ef5675", "CK" = "#ff764a", "AG" = "#ffa600")
cd_plot <- ggplot(df, aes(x = year)) +
  geom_point(aes(y = PV, color = "PV")) +
  geom_point(aes(y = FE, color = "FE")) +
  geom_point(aes(y = CM, color = "CM")) +
  geom_point(aes(y = PT, color = "PT")) +
  geom_point(aes(y = JT, color = "JT")) +
  geom_point(aes(y = CK, color = "CK")) +
  geom_point(aes(y = AG, color = "AG")) +
  scale_color_manual(values = cd_colors, name = "TCA") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey")) + 
  xlab("Year") +
  ylab("Log-density (tortoises / km-squared") + 
  ggtitle("Colorado Desert tortoise log-densities")
ggsave("./plots/tortoiseDensities/colorado_desert.png", cd_plot)

# Average CD densities
y_cd <- (df$PV + df$FE + df$CM + df$PT + df$JT + df$CK + df$AG) * 1 / 7
fit_cd <- lm(y_cd ~ df$year)

# Isolate fit
b_cd <- fit_cd$coefficients[1]
m_cd <- fit_cd$coefficients[2]
line_cd <- data.frame(year = df$year, y = m_cd * df$year + b_cd)

# Plot CD log-densities
cd_fit_colors <- c("PV" = "#003f5c", "FE" = "#374c80", "CM" = "#7a5195", "PT" = "#bc5090", 
              "JT" = "#ef5675", "CK" = "#ff764a", "AG" = "#ffa600", "Fit" = "black")
cd_fit_plot <- ggplot(df, aes(x = year)) +
  geom_point(aes(y = PV, color = "PV")) +
  geom_point(aes(y = FE, color = "FE")) +
  geom_point(aes(y = CM, color = "CM")) +
  geom_point(aes(y = PT, color = "PT")) +
  geom_point(aes(y = JT, color = "JT")) +
  geom_point(aes(y = CK, color = "CK")) +
  geom_point(aes(y = AG, color = "AG")) +
  geom_line(aes(y = y, color = "Fit"), data = line_cd) +
  scale_color_manual(values = cd_colors, name = "TCA") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey")) + 
  xlab("Year") +
  ylab("Log-density (tortoises / km-squared") + 
  ggtitle("Colorado Desert tortoise log-densities")
ggsave("./plots/tortoiseDensities/colorado_desert_fit.png", cd_fit_plot)

# Plot EM log-densities
em_colors <- c("EV" = "#003f5c", "IV" = "#bc5090")
em_plot <- ggplot(df, aes(x = year)) +
  geom_point(aes(y = EV, color = "EV")) +
  geom_point(aes(y = IV, color = "IV")) +
  scale_color_manual(values = em_colors, name  = "TCA") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey")) +
  xlab("Year") +
  ylab("Log-density (tortoises / km-sqared)") + 
  ggtitle("Eastern Mojave tortoise log-densities")
ggsave("./plots/tortoiseDensities/eastern_mojave.png", em_plot)

# Average EM densities
y_em <- (df$EV + df$IV) * 1 / 2
fit_em <- lm(y_em ~ df$year)

# Isolate fit
b_em <- fit_em$coefficients[1]
m_em <- fit_em$coefficients[2]
line_em <- data.frame(year = df$year, y = df$year * m_em + b_em)

# Plot fit
em_fit_colors <- c("EV" = "#003f5c", "IV" = "#bc5090", "Fit" = "black")
em_fit_plot <- ggplot(df, aes(x = year)) +
  geom_point(aes(y = EV, color = "EV")) +
  geom_point(aes(y = IV, color = "IV")) +
  geom_line(aes(y = y, color = "Fit"), data = line_em) +
  scale_color_manual(values = em_fit_colors, name  = "TCA") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey")) +
  xlab("Year") +
  ylab("Log-density (tortoises / km-sqared)") + 
  ggtitle("Eastern Mojave tortoise log-densities")
ggsave("./plots/tortoiseDensities/eastern_mojave_fit.png", em_fit_plot)

# Plot NM log-densities
nm_colors <- c("CS" = "#003f5c", "MM" = "#374c80", "GB" = "#7a5195", "BD" = "#bc5090")
nm_plot <- ggplot(df, aes(x = year)) +
  geom_point(aes(y = CS, color = "CS")) +
  geom_point(aes(y = MM, color = "MM")) +
  geom_point(aes(y = GB, color = "GB")) +
  geom_point(aes(y = BD, color = "BD")) +
  scale_color_manual(values = nm_colors, name = "TCA") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey")) + 
  xlab("Year") +
  ylab("Log-density (tortoises / km-squared") + 
  ggtitle("Northeastern Mojave tortoise log-densities")
ggsave("./plots/tortoiseDensities/northeastern_mojave.png", nm_plot)

# Average NM densities
y_nm <- (df$CS + df$MM + df$GB + df$BD) * 1 / 4
fit_nm <- lm(y_nm ~ df$year)

# Isolate fit
b_nm <- fit_nm$coefficients[1]
m_nm <- fit_nm$coefficients[2]
line_nm <- data.frame(year = df$year, y = m_nm * df$year + b_nm)

# Plot NM log-densities
nm_fit_colors <- c("CS" = "#003f5c", "MM" = "#374c80", "GB" = "#7a5195", "BD" = "#bc5090",
                   "Fit" = "black")
nm_fit_plot <- ggplot(df, aes(x = year)) +
  geom_point(aes(y = CS, color = "CS")) +
  geom_point(aes(y = MM, color = "MM")) +
  geom_point(aes(y = GB, color = "GB")) +
  geom_point(aes(y = BD, color = "BD")) +
  geom_line(aes(y = y, color = "Fit"), data = line_nm) +
  scale_color_manual(values = nm_colors, name = "TCA") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey")) + 
  xlab("Year") +
  ylab("Log-density (tortoises / km-squared") + 
  ggtitle("Northeastern Mojave tortoise log-densities")
ggsave("./plots/tortoiseDensities/northeastern_mojave_fit.png", nm_fit_plot)