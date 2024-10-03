# ADH visualization
library(ggpubr)
library(ggplot2)
library(scales)
library(ggnewscale)

### ADH calculations
#Summer 2022
summeradh <- read.csv("Summer Temp Data for ADH (pre).csv")
summeradh$adjusted_temperature <- pmax(summeradh$Temperature - summeradh$Baseline, 0)
summeradh$adh <- c(0, cumsum(summeradh$adjusted_temperature[-1]))
summary(summeradh)

write.csv(summeradh, "summeradh.csv", row.names = FALSE)


#Fall 2022

falladh <- read.csv("Fall Temp Data for ADH (pre).csv")
falladh$adjusted_temperature <- pmax(falladh$Temperature - falladh$Baseline, 0)
falladh$adh <- c(0, cumsum(falladh$adjusted_temperature[-1]))
summary(falladh)

write.csv(falladh, "falladh.csv", row.names = FALSE)


#Spring 2023

springadh <- read.csv("Spring Temp Data for ADH (pre).csv")
springadh$adjusted_temperature <- pmax(springadh$Temperature - springadh$Baseline, 0)
springadh$adh <- c(0, cumsum(springadh$adjusted_temperature[-1]))
summary(springadh)

write.csv(springadh, "springadh.csv", row.names = FALSE)



summertaph <- read.csv("Summer 2022 Taphonomy.csv")

stage_order <- c("Dry", "Advanced", "Active", "Bloat", "Fresh")

summertaph$Stage <- factor(summertaph$Stage, levels = stage_order)

my_palette <- viridis_pal()(length(unique(summertaph$Stage)))

ggplot(summertaph, aes(x = ADHSUB, y = factor(Carcass, levels = rev(c("C1", "C2", "C3", "C4", "C5", "C6"))), fill = Stage)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_palette) +
  theme_pub() +
  labs(title = "Decomposition Stages for Each Pig Carcass",
       x = "Timepoint",
       y = "Carcass ID",
       fill = "Stage") +
  theme(
    text = element_text(family = "Econ Sans Cnd")
  )

# ...

# ALL

alltaph <- read.csv("Taphonomy ADH all seasons.csv")

alltaph$Stage <- factor(alltaph$Stage, levels = stage_order)

alltaph$SEASON <- factor(alltaph$SEASON,
                                     levels = c("Fall 2022", "Summer 2022", "Spring 2023"))

carcass_study <- data.frame(
  Carcass = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C16", "C17", "C18", "C19", "C20", "C21"),
  Study = c(rep("Summer", 6), rep("Fall", 6), rep("Spring", 6))
)

# Define color palettes
stage_palette <- my_palette
season_palette <- c("Summer" = "#6ca3ff", "Fall" = "#b11509", "Spring" = "#a4d13a")
season_stage_palette <- c("Fresh" = "#FDE725FF", "Bloat" = "#5DC863FF", "Active" = "#21908CFF", "Advanced" = "#3B528BFF", "Dry" = "#440154FF","Summer" = "#6ca3ff", "Fall" = "#b11509", "Spring" = "#a4d13a")


# Season backgrounds data frame
season_backgrounds <- data.frame(
  ymin = c(0.5, 6.5, 12.5),
  ymax = c(6.5, 12.5, 18.5),
  Season = c("Summer", "Fall", "Spring")
)

all_taphonomy_per_carcass <- ggplot(alltaph, aes(x = ADHSUB, y = factor(Carcass, levels = rev(c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C16", "C17", "C18", "C19", "C20", "C21"))))) +
  geom_rect(data = season_backgrounds,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = Season),
            alpha = 0.2,
            inherit.aes = FALSE) +
  scale_fill_manual(values = season_palette, name = "Season") +
  new_scale_fill() +
  geom_bar(aes(fill = Stage), stat = "identity") +
  scale_fill_manual(values = stage_palette, name = "Stage") +
  labs(x = "ADH",
       y = "Carcass ID") +
  theme_pub() +
  theme(
    text = element_text(family = "Econ Sans Cnd"),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(
    fill = guide_legend(reverse = TRUE, keywidth = 3, keyheight = 1),
    fill_new = guide_legend(reverse = TRUE, keywidth = 3, keyheight = 1)
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))

all_taphonomy_per_carcass

# Range plots

alltaphonomyrange <- read.csv("Taphonomy ADH all seasons proper ranges.csv")

alltaphonomyrange <- alltaphonomyrange %>%
  filter_all(all_vars(!grepl("N/A", .)))

alltaphonomyrange$ADHCUM <- as.numeric(alltaphonomyrange$ADHCUM)


range_data <- alltaphonomyrange %>%
  group_by(SEASON, Stage) %>%
  summarize(mean_ADHSUB = mean(ADHSUB),
            mean_ADHCUM = mean(ADHCUM),
            se_ADHSUB = sd(ADHSUB) / sqrt(n()),
            se_ADHCUM = sd(ADHCUM) / sqrt(n()))

range_data$Stage <- factor(range_data$Stage, levels = c("Fresh", "Bloat", "Active", "Advanced", "Dry"))


create_range_plot <- function(data, season) {
  ggplot(data, aes(x = Stage, ymin = mean_ADHSUB, ymax = mean_ADHCUM, group = 1, color = Stage)) +
    geom_linerange(position = position_dodge(width = 0.9), alpha = 1, size = 3) +
    geom_errorbar(aes(ymin = mean_ADHSUB - ifelse(se_ADHSUB != 0, se_ADHSUB, 0),
                      ymax = mean_ADHSUB + ifelse(se_ADHSUB != 0, se_ADHSUB, 0)),
                  position = position_dodge(width = 0.9),
                  width = 0, linetype = 1, length = unit(0.05, "npc")) +  # Adjust length here
    geom_errorbar(aes(ymin = mean_ADHCUM - ifelse(se_ADHCUM != 0, se_ADHCUM, 0),
                      ymax = mean_ADHCUM + ifelse(se_ADHCUM != 0, se_ADHCUM, 0)),
                  position = position_dodge(width = 0.9),
                  width = 0, linetype = 1, length = unit(0.05, "npc")) +  # Adjust length here
    geom_point(aes(y = mean_ADHSUB, color = Stage), position = position_dodge(width = 0.9), size = 5) +
    geom_point(aes(y = mean_ADHCUM, color = Stage), position = position_dodge(width = 0.9), size = 5) +
    labs(title = season, x = "Decomposition Stage", y = "Accumulated Degree Hours (ADH)") +
    theme_pub() +
    theme(
          legend.position = "none",
          text = element_text(family = "Econ Sans Cnd", size = 14),
          axis.title.y = element_text(margin = margin(r = 10)),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14)) +
    scale_y_continuous(limits = c(0, 2000)) +
    scale_color_manual(values = c("Fresh" = "#fde725", "Bloat" = "#5dc862", "Active" = "#20908c",
                                  "Advanced" = "#3b528b", "Dry" = "#440d54"))
}


spring_range_plot <- create_range_plot(subset(range_data, SEASON == "SPRING"), "Spring 2023")
summer_range_plot <- create_range_plot(subset(range_data, SEASON == "SUMMER"), "Summer 2022")
fall_range_plot <- create_range_plot(subset(range_data, SEASON == "FALL"), "Fall 2022")

# Arrange range plots
combined_range_plots <- plot_grid(spring_range_plot, summer_range_plot, fall_range_plot, nrow = 1)

# Print the combined range plots
print(combined_range_plots)

ggsave("decomprange.png", width = 12, height = 9, units = "in", dpi = 300, path = "/Volumes/GRIGSBYGRAD/Plots for publication 2024/Decomp range", combined_range_plots)


alltaphonomyrange$Stage <- factor(alltaphonomyrange$Stage, levels = c("Fresh", "Bloat", "Active", "Advanced", "Dry"))

create_rangeish_plot <- function(data, season) {
  ggplot(data, aes(x = Stage)) +
    geom_point(aes(y = ADHSUB, shape = "ADHSUB"), color = "black", alpha = 0.7, size = 2) +
    geom_point(aes(y = ADHCUM, shape = "ADHCUM"), color = "black", alpha = 0.7, size = 2) +
    labs(title = season, x = "Decomposition Stage", y = "Accumulated Degree Hours (ADH)") +
    scale_shape_manual(name = "", values = c(ADHSUB = 3, ADHCUM = 1),  # 3 for X and 1 for circle
                       labels = c("Upper limit", "Lower limit")) +
    theme_pub() +
    theme(legend.position = c(0.95, 0.05),  # Position the legend in bottom right corner
          legend.justification = c(1, 0),   # Anchor point of the legend box
          legend.box.just = "right",        # Align legend text to the right
          axis.title.y = element_text(margin = margin(r = 10)),
          plot.title = element_text(size = 30),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_blank(),  # Remove plot background
          panel.background = element_blank(),
          text = element_text(family = "Econ Sans Cnd", size = 14),
    ) +
    scale_y_continuous(limits = c(0, 2000)) +
    guides(shape = guide_legend(ncol = 1))  # Set legend to one column
}

spring_rangeish_plot <- create_rangeish_plot(subset(alltaphonomyrange, SEASON == "SPRING"), "Spring")
summer_rangeish_plot <- create_rangeish_plot(subset(alltaphonomyrange, SEASON == "SUMMER"), "Summer")
fall_rangeish_plot <- create_rangeish_plot(subset(alltaphonomyrange, SEASON == "FALL"), "Fall")

# Arrange range plots
combined_rangeish_plots <- plot_grid(spring_rangeish_plot, summer_rangeish_plot, fall_rangeish_plot, nrow = 1)

# Print the combined range plots
print(combined_rangeish_plots)

ggsave("decompdots.png", width = 12, height = 9, units = "in", dpi = 300, path = "/Volumes/GRIGSBYGRAD/Plots for publication 2024/Decomp range", combined_rangeish_plots)




### Combining plots
combined_adh_data_plot <- (combined_range_plots / all_taphonomy_per_carcass) +
  plot_layout(widths = c(1, 1))  # Adjust the relative widths as needed

combined_adh_data_plot



### ADH over time

summeradh <- read.csv("summeradh.csv")
falladh <- read.csv("falladh.csv")
springadh <- read.csv("springadh.csv")


# Function to create ADH plot
create_adh_plot <- function(data, season, color) {
  ggplot(data, aes(x = Hour, y = adh)) +
    geom_line(color = color, size = 1) +
    labs(x = "Hours", y = "ADH") +
    theme_pub() +
    theme(
      text = element_text(family = "Econ Sans Cnd"),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 14),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 2000))
}

fall_adh_plot <- create_adh_plot(falladh, "Fall 2022", "#b11509")
fall_adh_plot

spring_adh_plot <- create_adh_plot(springadh, "Spring 2023", "#a4d13a")
spring_adh_plot

summer_adh_plot <- create_adh_plot(summeradh, "Summer 2022", "#6ca3ff")
summer_adh_plot



combined_adh_vs_time_plot <- ( spring_adh_plot / summer_adh_plot / fall_adh_plot )
combined_adh_vs_time_plot


