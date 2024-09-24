# Analyzing and visualizing the seasonal weather data
#################################################################
# Library dependencies
library(ggplot2)
library(dplyr)
library(dunn.test)
library(ggthemes)
library(scales)
library(patchwork)

# Plotting theme
theme_pub <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = margin(0, 0, 0, 0),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))

}

# Standard error calculation function of choice
standard_error <- function(x) {
  sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
}

###### SPRING
spring_weather_data <- read.csv('Weather Data/Spring historical weather data.csv')

# Calculate average and standard error for temperature and humidity
spring_summary_stats <- spring_weather_data %>%
  summarise(
    average_temperature = mean(Temp, na.rm = TRUE),
    se_temperature = standard_error(Temp),
    average_humidity = mean(RH., na.rm = TRUE),
    se_humidity = standard_error(RH.)
  )

# Print the results
print(spring_summary_stats)

### TEMPERATURE
# Checking normality
shapiro.test(spring_weather_data$Temp[spring_weather_data$Year == 2020])
shapiro.test(spring_weather_data$Temp[spring_weather_data$Year == 2021])
shapiro.test(spring_weather_data$Temp[spring_weather_data$Year == 2022])
shapiro.test(spring_weather_data$Temp[spring_weather_data$Year == 2023])

# Kruskal-Wallis test
kruskal.test(Temp ~ Year, data = spring_weather_data)

# Dunn's test
dunn.test(spring_weather_data$Temp, spring_weather_data$Year, method = "bonferroni")

### HUMIDITY
# Checking normality
shapiro.test(spring_weather_data$RH.[spring_weather_data$Year == 2020])
shapiro.test(spring_weather_data$RH.[spring_weather_data$Year == 2021])
shapiro.test(spring_weather_data$RH.[spring_weather_data$Year == 2022])
shapiro.test(spring_weather_data$RH.[spring_weather_data$Year == 2023])

# Kruskal-Wallis test
kruskal.test(RH. ~ Year, data = spring_weather_data)

# Dunn's test
dunn.test(spring_weather_data$RH., spring_weather_data$Year, method = "bonferroni")


###### SUMMER
summer_weather_data <- read.csv('Weather Data/Summer historical weather data.csv')

summer_summary_stats <- summer_weather_data %>%
  summarise(
    average_temperature = mean(Temp, na.rm = TRUE),
    se_temperature = standard_error(Temp),
    average_humidity = mean(RH., na.rm = TRUE),
    se_humidity = standard_error(RH.)
  )

print(summer_summary_stats)

### TEMPERATURE
shapiro.test(summer_weather_data$Temp[summer_weather_data$Year == 2019])
shapiro.test(summer_weather_data$Temp[summer_weather_data$Year == 2020])
shapiro.test(summer_weather_data$Temp[summer_weather_data$Year == 2021])
shapiro.test(summer_weather_data$Temp[summer_weather_data$Year == 2022])

kruskal_test_result <- kruskal.test(Temp ~ Year, data = summer_weather_data)
kruskal_test_result

dunn.test(summer_weather_data$Temp, summer_weather_data$Year, method = "bonferroni")

### HUMIDITY
shapiro.test(summer_weather_data$RH.[summer_weather_data$Year == 2019])
shapiro.test(summer_weather_data$RH.[summer_weather_data$Year == 2020])
shapiro.test(summer_weather_data$RH.[summer_weather_data$Year == 2021])
shapiro.test(summer_weather_data$RH.[summer_weather_data$Year == 2022])

kruskal.test(RH. ~ Year, data = summer_weather_data)

dunn.test(summer_weather_data$RH., summer_weather_data$Year, method = "bonferroni")


###### FALL
fall_weather_data <- read.csv('Weather Data/Fall historical weather data.csv')

fall_summary_stats <- fall_weather_data %>%
  summarise(
    average_temperature = mean(Temp, na.rm = TRUE),
    se_temperature = standard_error(Temp),
    average_humidity = mean(RH., na.rm = TRUE),
    se_humidity = standard_error(RH.)
  )

print(fall_summary_stats)

### TEMPERATURE
shapiro.test(fall_weather_data$Temp[fall_weather_data$Year == 2019])
shapiro.test(fall_weather_data$Temp[fall_weather_data$Year == 2020])
shapiro.test(fall_weather_data$Temp[fall_weather_data$Year == 2021])
shapiro.test(fall_weather_data$Temp[fall_weather_data$Year == 2022])

kruskal.test(Temp ~ Year, data = fall_weather_data)

dunn.test(fall_weather_data$Temp, fall_weather_data$Year, method = "bonferroni")

### HUMIDITY
shapiro.test(fall_weather_data$RH.[fall_weather_data$Year == 2019])
shapiro.test(fall_weather_data$RH.[fall_weather_data$Year == 2020])
shapiro.test(fall_weather_data$RH.[fall_weather_data$Year == 2021])
shapiro.test(fall_weather_data$RH.[fall_weather_data$Year == 2022])

kruskal.test(RH. ~ Year, data = fall_weather_data)

dunn.test(fall_weather_data$RH., fall_weather_data$Year, method = "bonferroni")


###### COMPARING ALL SEASONS (STUDY YEARS)
all_weather_data <- read.csv('Weather Data/Studies weather data compiled.csv')

### TEMPERATURE
shapiro.test(all_weather_data$Temperature[all_weather_data$Study == "Summer 2022"])
shapiro.test(all_weather_data$Temperature[all_weather_data$Study == "Fall 2022"])
shapiro.test(all_weather_data$Temperature[all_weather_data$Study == "Spring 2023"])

kruskal.test(Temperature ~ Study, data = all_weather_data)

dunn.test(all_weather_data$Temperature, all_weather_data$Study, method = "bonferroni")

### HUMIDITY
shapiro.test(all_weather_data$RH.[all_weather_data$Study == "Summer 2022"])
shapiro.test(all_weather_data$RH.[all_weather_data$Study == "Fall 2022"])
shapiro.test(all_weather_data$RH.[all_weather_data$Study == "Spring 2023"])

kruskal.test(RH. ~ Study, data = all_weather_data)

dunn.test(all_weather_data$RH., all_weather_data$Study, method = "bonferroni")



######### VISUALIZATIONS

###### WEATHER DATA PLOTS

### SUMMER
# Temperature plot
summer_weather_data$datetime <- as.POSIXct(paste("2022", summer_weather_data$Date, summer_weather_data$Time), format="%Y %m-%d %H:%M")

summer_weather_data <- summer_weather_data %>%
  mutate(Temp = as.numeric(Temp),
         Precipitation = as.numeric(Precipitation),
         RH = as.numeric(`RH.`))

summerhistoricalweatherplot <- ggplot(summer_weather_data, aes(x = datetime, y = Temp)) +
  geom_line(data = filter(summer_weather_data, Year != 2022), aes(color = interaction(Year), linetype = interaction(Year)), alpha = 0.9) +
  geom_line(data = filter(summer_weather_data, Year == 2022), aes(color = interaction(Year), linetype = interaction(Year)), linewidth = 1.2) +
  scale_color_manual(values = c("darkgrey", "darkgrey", "darkgrey", "#6ca3ff"), name = "Year") +
  scale_linetype_manual(values = c(4, 2, 3, 1), name = "Year") +
  labs(
    x = "Date",
    y = "Temperature (°C)"
  ) +
  scale_x_datetime(
    expand = expansion(mult = c(0, 0)),
    date_breaks = "1 day",
    date_minor_breaks = "6 hours"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)),
                     limits = c(-5, 35)) +
  theme_pub() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    legend.key.size = unit(.5, "cm"),
    axis.text.y = element_text(size = 18)
  )

summerhistoricalweatherplot

# Humidity + precipitation (2022)
summer_data_2022 <- filter(summer_weather_data, Year == 2022)
max_precipitation <- max(summer_data_2022$Precipitation, na.rm = TRUE)
max_rh <- max(summer_data_2022$RH, na.rm = TRUE)

scaling_factor <- 110 / max_rh

summerweathervariableplot <- ggplot(summer_data_2022, aes(x = datetime)) +
  geom_ribbon(aes(ymin = 0, ymax = RH, fill = "Relative Humidity"), alpha = 0.3) +
  geom_line(aes(y = RH, color = "Relative Humidity"), size = .7, linetype = 1, show.legend = FALSE) +
  geom_bar(aes(y = Precipitation, fill = "Precipitation"), stat = "identity", alpha = 0.9) +
  scale_y_continuous(
    name = "Precipitation (mm)",
    expand = expansion(mult = c(0, 0)),
    limits = c(0, 110),
    sec.axis = sec_axis(
      transform = ~ .,
      name = "Relative Humidity (%)",
      breaks = seq(0, 110, by = 25)
    )
  ) +
  scale_x_datetime(
    expand = expansion(mult = c(0, 0)),
    date_breaks = "1 day",
    date_minor_breaks = "6 hours"
  ) +
  scale_fill_manual(values = c("Relative Humidity" = "#f8776e", "Precipitation" = "#6ca3ff")) +
  scale_color_manual(values = c("Relative Humidity" = "#f8776e")) +
  labs(
    x = "Date",
    fill = "Weather Variable",
  ) +
  theme_pub() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    legend.key.size = unit(.5, "cm"),
    axis.text.y = element_text(size = 18)
  )

summerweathervariableplot

# Comibination
summerweatherplot_combo <- (summerhistoricalweatherplot | summerweathervariableplot)

print(summerweatherplot_combo)

ggsave("summerweatherplot_combo.png", width = 30, height = 10, units = "cm", dpi = 300, path = "", summerweatherplot_combo)

###### We can turn this into functions for ease of use with other seasonal datasets
### Make sure you adjust the plot colors in the function to suit your needs

historical_weather_plot <- function(data) {
  ggplot(summer_weather_data, aes(x = datetime, y = Temp)) +
    geom_line(data = filter(summer_weather_data, Year != 2022), aes(color = interaction(Year), linetype = interaction(Year)), alpha = 0.9) +
    geom_line(data = filter(summer_weather_data, Year == 2022), aes(color = interaction(Year), linetype = interaction(Year)), linewidth = 1.2) +
    scale_color_manual(values = c("darkgrey", "darkgrey", "darkgrey", "#6ca3ff"), name = "Year") +
    scale_linetype_manual(values = c(4, 2, 3, 1), name = "Year") +
    labs(
      x = "Date",
      y = "Temperature (°C)"
    ) +
    scale_x_datetime(
      expand = expansion(mult = c(0, 0)),
      date_breaks = "1 day",
      date_minor_breaks = "6 hours"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)),
                       limits = c(-5, 35)) +
    theme_pub() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 18),
      axis.title.x = element_text(margin = margin(t = 5)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      legend.key.size = unit(.5, "cm"),
      axis.text.y = element_text(size = 18)
    )
}

weather_humidity_precipitation_plot <- function(data) {
  ggplot(summer_data_2022, aes(x = datetime)) +
    geom_ribbon(aes(ymin = 0, ymax = RH, fill = "Relative Humidity"), alpha = 0.3) +
    geom_line(aes(y = RH, color = "Relative Humidity"), size = .7, linetype = 1, show.legend = FALSE) +
    geom_bar(aes(y = Precipitation, fill = "Precipitation"), stat = "identity", alpha = 0.9) +
    scale_y_continuous(
      name = "Precipitation (mm)",
      expand = expansion(mult = c(0, 0)),
      limits = c(0, 110),
      sec.axis = sec_axis(
        transform = ~ .,
        name = "Relative Humidity (%)",
        breaks = seq(0, 110, by = 25)
      )
    ) +
    scale_x_datetime(
      expand = expansion(mult = c(0, 0)),
      date_breaks = "1 day",
      date_minor_breaks = "6 hours"
    ) +
    scale_fill_manual(values = c("Relative Humidity" = "#f8776e", "Precipitation" = "#6ca3ff")) +
    scale_color_manual(values = c("Relative Humidity" = "#f8776e")) +
    labs(
      x = "Date",
      fill = "Weather Variable",
    ) +
    theme_pub() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 18),
      axis.title.x = element_text(margin = margin(t = 5)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      legend.key.size = unit(.5, "cm"),
      axis.text.y = element_text(size = 18)
    )
}

# Example usage
summer_weather_data$datetime <- as.POSIXct(paste("2022", summer_weather_data$Date, summer_weather_data$Time), format="%Y %m-%d %H:%M")

summer_weather_data <- summer_weather_data %>%
  mutate(Temp = as.numeric(Temp),
         Precipitation = as.numeric(Precipitation),
         RH = as.numeric(`RH.`))

historical_weather_plot(summer_weather_data)

summer_data_2022 <- filter(summer_weather_data, Year == 2022)
max_precipitation <- max(summer_data_2022$Precipitation, na.rm = TRUE)
max_rh <- max(summer_data_2022$RH, na.rm = TRUE)

scaling_factor <- 110 / max_rh

weather_humidity_precipitation_plot(summer_data_2022)



### Temperature boxplots comparing all seasonal weather data

COMPILED_weatherdata <- read.csv('Weather Data/Studies weather data compiled.csv')

weather_boxplots <- ggplot(COMPILED_weatherdata, aes(x = Study, y = Temperature, color = Study, fill = Study)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.7) +
  stat_summary(fun = median, geom = "crossbar", size =.5, alpha = 0.7) +
  labs(x = "", y = "Temperature (°C)") +
  theme_pub() +
  theme(axis.text.x = element_text(size = 20, face = "bold")) +
  theme(axis.text.y = element_text(size = 18)) +
  theme(plot.background = element_rect(fill = "white")) +
  theme(panel.background = element_rect(fill = "white")) +
  theme(legend.background = element_rect(fill = "white")) +
  theme(text = element_text(family = "Helvetica")) +
  theme(axis.title = element_text(size = 20)) +
  guides(fill = "none", color = "none")

weather_boxplots

ggsave("weather_boxplots.png", width = 25, height = 25, units = "cm", dpi = 300, path = "", weather_boxplots)

### Temperature rainfall plot

weather_rainfall_plot <- ggplot(COMPILED_weatherdata, aes(x = Study, y = Temperature, fill = Study, color = Study)) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.7, size = 1.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7, color = "black") +  # Set boxplot outline to black
  ggdist::stat_halfeye(adjust = 0.5, width = 0.6, justification = -0.2, .width = 0, point_colour = NA) +
  labs(x = "", y = "Temperature (°C)") +  # Swap the labels for flipped axes
  theme_pub() +
  theme(
    axis.text.y = element_text(size = 20, face = "bold"),  # Adjusted to y-axis
    axis.text.x = element_text(size = 18),  # Adjusted to x-axis
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white"),
    text = element_text(family = "Helvetica"),
    axis.title = element_text(size = 20)
  ) +
  guides(fill = "none", color = "none") +  # Remove legends
  coord_flip()  # Flip the coordinates

plot(weather_rainfall_plot)

ggsave("weather_rainfall_plot.png", width = 25, height = 25, units = "cm", dpi = 300, path = "", weather_rainfall_plot)
