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
colors = c("FK" = "#003f5c", "SC" = "#bc5090", "OR" = "#ffa600")
wm_plot <- ggplot(df, aes(x = year)) +
  geom_point(aes(y = FK, color = "FK")) +
  geom_point(aes(y = SC, color = "SC")) +
  geom_point(aes(y = OR, color = "OR")) +
  scale_color_manual(values = colors, name  = "TCA") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey"),
        legend.title = element_text(face = "bold")) +
  xlab("Year") +
  ylab("Log-density (tortoises / km-sqared)") + 
  ggtitle("Western Mojave tortoise log-densities")
ggsave("./plots/tortoiseDensities/western_mojave.png", wm_plot)