# ANCOM-BC was performed in QIIME2 and the data was extracted for manipulation and visualization
# For help developing an ANCOM analysis refer to QIIME2 documentation

# Libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggthemes)
library(scales)
library(patchwork)
library(extrafont)
library(stringr)
library(wesanderson)

font_import()
loadfonts()

# Custom plotting theme
theme_pub <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(family = "Econ Sans Cnd"),
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

scale_fill_pub <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

scale_colour_pub <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}




### ANCOM visualization

# Load data
lfc_data <- read.csv("lfc_slice.csv")
q_val_data <- read.csv("q_val_slice.csv")
se_data <- read.csv("se_slice.csv")

# Remove rows where ID ends with '__'
lfc_data <- lfc_data %>% filter(!grepl("__$", id))
se_data <- se_data %>% filter(!grepl("__$", id))
q_val_data <- q_val_data %>% filter(!grepl("__$", id))

# Pivot data
melted_lfc_data <- pivot_longer(lfc_data,
                                cols = starts_with("Stage"),
                                names_to = "Stage",
                                values_to = "LogFoldChange")

melted_se_data <- pivot_longer(se_data,
                               cols = starts_with("Stage"),
                               names_to = "Stage",
                               values_to = "SE")

melted_q_val_data <- pivot_longer(q_val_data,
                                  cols = starts_with("Stage"),
                                  names_to = "Stage",
                                  values_to = "Q")

# Extract taxon IDs
melted_lfc_data$id <- sub(".*;g__", "", melted_lfc_data$id)
melted_se_data$id <- sub(".*;g__", "", melted_se_data$id)
melted_q_val_data$id <- sub(".*;g__", "", melted_q_val_data$id)

# Fix stage names
melted_lfc_data$Stage <- gsub("Stage", "", melted_lfc_data$Stage)
melted_se_data$Stage <- gsub("Stage", "", melted_se_data$Stage)
melted_q_val_data$Stage <- gsub("Stage", "", melted_q_val_data$Stage)

# Merge data
merged_ancom_data <- merge(melted_lfc_data, melted_se_data, by = c("id", "Stage"))
merged_ancom_data <- merge(merged_ancom_data, melted_q_val_data, by = c("id", "Stage"))
head(merged_ancom_data)

# Filter data
merged_ancom_data_filtered <- merged_ancom_data %>%
  group_by(Stage) %>%
  arrange(LogFoldChange) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 10 | rank > n() - 10) %>%
  ungroup()

# Plotting
ancom_plots <- list()

for (stage in unique(merged_ancom_data_filtered$Stage)) {
  ancom_plot_data <- subset(merged_ancom_data_filtered, Stage == stage)

  # Add asterisk to the ID if it's significant
  ancom_plot_data$significant <- ifelse(ancom_plot_data$Q < 0.05, "*", "")

  # Concatenate asterisk with ID
  ancom_plot_data$id_with_significance <- paste0(ancom_plot_data$significant, ancom_plot_data$id)

  # Reorder id factor levels based on the sign of LogFoldChange
  ancom_plot_data$id_with_significance <- reorder(ancom_plot_data$id_with_significance, +ancom_plot_data$LogFoldChange)

  fdp <- ggplot(ancom_plot_data, aes(x = id_with_significance, y = LogFoldChange, fill = LogFoldChange)) +
    geom_bar(stat = "identity", position = "identity", width = 0.8) +
    geom_errorbar(aes(ymin = LogFoldChange - SE, ymax = LogFoldChange + SE), width = 0.2, position = position_dodge(0.5)) +
    labs(title = "",
         x = element_blank(),
         y = "Log Fold Change (LFC)") +
    scale_fill_gradient2(low = "#fd533d", mid = "white", high = "#5898c9", midpoint = 0, limits = c(-9, 9),  breaks = seq(-9, 9, by = 3)) +  # Color scale
    theme_pub() +
    theme(legend.position = "bottom",
          panel.grid.major = element_line(color = "#ececec", linetype = "solid"),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.key.size = unit(0.9, "lines"),  # Adjust legend key size
          legend.text = element_text(size = 10),
          legend.direction = "horizontal",  # Horizontal gradient legend
          legend.box = "horizontal",
          legend.margin = margin(-7, 0, 4, 0),
          text = element_text(family = "Helvetica"),
          plot.title = element_text(size = 30, hjust = 0),
          plot.margin = margin(0, 0, -0.4, 0),
          axis.text.y = element_text(vjust = .4, hjust = 1, size = 14),
          axis.text.x = element_text(size = 14)) +
    coord_flip() +
    scale_y_continuous(breaks = seq(-10, 10, by = 1), limits = c(-10, 10))

  ancom_plots[[stage]] <- fdp
}

ancom_plots[[1]]
ancom_plots[[2]]
ancom_plots[[3]]
ancom_plots[[4]]


combined_ancom_plots <- wrap_plots(ancom_plots[[4]], ancom_plots[[3]], ancom_plots[[1]], ancom_plots[[2]])

combined_ancom_plots

ggsave("combined_ancom_plots.png", width = 20, height = 20, units = "cm", dpi = 300, combined_ancom_plots)


### Optional

combined_ancom_plots_title <- ggdraw() +
  draw_label("_______ 2022", size = 32,
             color = "black", fontface = "bold", fontfamily = "Helvetica", hjust = 2.38, vjust = 1)

# Arrange the title and the combined plot
combined_ancom_plots_titled <- plot_grid(combined_ancom_plots_title, combined_ancom_plots, ncol = 1, align = "v", rel_heights = c(0.1, 1))

combined_ancom_plots_titled

ggsave("combined_ancom_plots_titled.png", width = 13, height = 9.5, units = "in", dpi = 300, combined_ancom_plots_titled)

### Export data as csv
write.csv(ancom_plot_data, "ancom_plot_data.csv")


### Alternative visualization, condenses information into one plot

# Merge data
merged_ancom_data <- merge(melted_lfc_data, melted_se_data, by = c("id", "Stage"))
merged_ancom_data <- merge(merged_ancom_data, melted_q_val_data, by = c("id", "Stage"))

# Filter data to include top 10 positive and top 10 negative LFC taxa
merged_ancom_data_filtered <- merged_ancom_data %>%
  group_by(id) %>%
  mutate(max_lfc = max(LogFoldChange),
         min_lfc = min(LogFoldChange)) %>%
  ungroup() %>%
  arrange(desc(max_lfc)) %>%
  slice_head(n = 10 * 4) %>%  # Top 10 positive LFC taxa * 4 stages
  bind_rows(
    merged_ancom_data %>%
      group_by(id) %>%
      mutate(max_lfc = max(LogFoldChange),
             min_lfc = min(LogFoldChange)) %>%
      ungroup() %>%
      arrange(min_lfc) %>%
      slice_head(n = 10 * 4)  # Top 10 negative LFC taxa * 4 stages
  )

# Add significance indicator
merged_ancom_data_filtered$significant <- ifelse(merged_ancom_data_filtered$Q < 0.05, "*", "")

merged_ancom_data_filtered <- merged_ancom_data_filtered %>%
  mutate(id = str_replace(id, "Clostridium_sensu_stricto_1", "Clostridium"))

# Reorder id factor levels based on the overall magnitude of LogFoldChange
id_order <- merged_ancom_data_filtered %>%
  group_by(id) %>%
  summarize(max_abs_lfc = max(abs(LogFoldChange))) %>%
  arrange(desc(max_abs_lfc)) %>%
  pull(id)

merged_ancom_data_filtered$id <- factor(merged_ancom_data_filtered$id,
                                        levels = rev(unique(id_order)))


stage_palette <- rev(c("#5dc863", "#21908c", "#3b528b", "#440154"))
stage_order <- rev(c("Bloat", "Active", "Advanced", "Dry"))
merged_ancom_data_filtered$Stage <- factor(merged_ancom_data_filtered$Stage, levels = stage_order)

ancom_plot <- ggplot(merged_ancom_data_filtered, aes(x = id, y = LogFoldChange, fill = Stage)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1), width = 1) +
  geom_errorbar(aes(ymin = LogFoldChange - SE, ymax = LogFoldChange + SE),
                position = position_dodge(width = 1), width = 0.25) +
  geom_text(aes(label = significant, y = ifelse(LogFoldChange >= 0, LogFoldChange + SE, LogFoldChange - SE)),
            position = position_dodge(width = 1),
            vjust = ifelse(merged_ancom_data_filtered$LogFoldChange >= 0, .7, .7),
            hjust = ifelse(merged_ancom_data_filtered$LogFoldChange >= 0, -0.5, 1.5),
            size = 4) +
  scale_fill_manual(values = stage_palette) +
  coord_flip() +
  labs(x = "Bacterial Taxa",
       y = "Log Fold Change (LFC)",
       fill = "Decomposition Stage") +
  theme_pub() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(face = "bold", size = 14),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5),
        text = element_text(family = "Econ Sans Cnd"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")
  ) +
  scale_y_continuous(breaks = seq(-7, 7, by = 2), limits = c(-7, 7)) +
  guides(fill = guide_legend(reverse = TRUE))


# Display the plot
print(ancom_plot)

# Save the plot
ggsave("combined_ancom_plot.png", ancom_plot, width = 20, height = 20, units = "cm", dpi = 300)

# Export data as csv
write.csv(merged_ancom_data_filtered, "ancom_plot_data.csv", row.names = FALSE)



### Making one for seasonal comparisons
# Load data
lfc_data <- read.csv("lfc_slice.csv")
q_val_data <- read.csv("q_val_slice.csv")
se_data <- read.csv("se_slice.csv")

# Remove rows where ID ends with '__'
lfc_data <- lfc_data %>% filter(!grepl("__$", id))
se_data <- se_data %>% filter(!grepl("__$", id))
q_val_data <- q_val_data %>% filter(!grepl("__$", id))

# Pivot data
melted_lfc_data <- pivot_longer(lfc_data,
                                cols = starts_with("Season"),
                                names_to = "Season",
                                values_to = "LogFoldChange")

melted_se_data <- pivot_longer(se_data,
                               cols = starts_with("Season"),
                               names_to = "Season",
                               values_to = "SE")

melted_q_val_data <- pivot_longer(q_val_data,
                                  cols = starts_with("Season"),
                                  names_to = "Season",
                                  values_to = "Q")

# Extract taxon IDs
melted_lfc_data$id <- sub(".*;g__", "", melted_lfc_data$id)
melted_se_data$id <- sub(".*;g__", "", melted_se_data$id)
melted_q_val_data$id <- sub(".*;g__", "", melted_q_val_data$id)

# Fix stage names
melted_lfc_data$Season <- gsub("Season", "", melted_lfc_data$Season)
melted_se_data$Season <- gsub("Season", "", melted_se_data$Season)
melted_q_val_data$Season <- gsub("Season", "", melted_q_val_data$Season)

# Merge data
merged_ancom_data <- merge(melted_lfc_data, melted_se_data, by = c("id", "Season"))
merged_ancom_data <- merge(merged_ancom_data, melted_q_val_data, by = c("id", "Season"))

# Filter data to include top 10 positive and top 10 negative LFC taxa
merged_ancom_data_filtered <- merged_ancom_data %>%
  group_by(id) %>%
  summarize(max_lfc = max(LogFoldChange),
            min_lfc = min(LogFoldChange)) %>%
  ungroup() %>%
  arrange(desc(max_lfc)) %>%
  slice_head(n = 10) %>%  # Top 10 positive LFC taxa
  bind_rows(
    merged_ancom_data %>%
      group_by(id) %>%
      summarize(max_lfc = max(LogFoldChange),
                min_lfc = min(LogFoldChange)) %>%
      ungroup() %>%
      arrange(min_lfc) %>%
      slice_head(n = 10)  # Top 10 negative LFC taxa
  ) %>%
  inner_join(merged_ancom_data, by = "id")

# Add significance indicator
merged_ancom_data_filtered$significant <- ifelse(merged_ancom_data_filtered$Q < 0.05, "*", "")

merged_ancom_data_filtered <- merged_ancom_data_filtered %>%
  mutate(id = str_replace(id, "Clostridium_sensu_stricto_1", "Clostridium"))

# Reorder id factor levels based on the overall magnitude of LogFoldChange
id_order <- merged_ancom_data_filtered %>%
  group_by(id) %>%
  summarize(max_abs_lfc = max(abs(LogFoldChange))) %>%
  arrange(desc(max_abs_lfc)) %>%
  pull(id)

merged_ancom_data_filtered$id <- factor(merged_ancom_data_filtered$id,
                                        levels = rev(unique(id_order)))

season_palette <- rev(c("#a4d13a", "#b11509"))

# merged_ancom_data_filtered$Season <- factor(merged_ancom_data_filtered$Season, levels = stage_order)

season_ancom_plot <- ggplot(merged_ancom_data_filtered, aes(x = id, y = LogFoldChange, fill = Season)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1), width = 1) +
  geom_errorbar(aes(ymin = LogFoldChange - SE, ymax = LogFoldChange + SE),
                position = position_dodge(width = 1), width = 0.25) +
  geom_text(aes(label = significant, y = ifelse(LogFoldChange >= 0, LogFoldChange + SE, LogFoldChange - SE)),
            position = position_dodge(width = 1),
            vjust = ifelse(merged_ancom_data_filtered$LogFoldChange >= 0, .7, .7),
            hjust = ifelse(merged_ancom_data_filtered$LogFoldChange >= 0, -0.5, 1.5),
            size = 4) +
  scale_fill_manual(values = season_palette) +
  coord_flip() +
  labs(x = "Bacterial Taxa",
       y = "Log Fold Change (LFC)",
       fill = "Decomposition Season") +
  theme_pub() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(face = "bold", size = 14),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5),
        text = element_text(family = "Econ Sans Cnd"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")
  ) +
  scale_y_continuous(breaks = seq(-7, 7, by = 2), limits = c(-7, 7)) +
  guides(fill = guide_legend(reverse = TRUE))


# Display the plot
print(season_ancom_plot)

# Save the plot
ggsave("combined_season_ancom_plot.png", season_ancom_plot, width = 22, height = 20, units = "cm", dpi = 300)

# Export data as csv
write.csv(merged_ancom_data_filtered, "season_ancom_plot_data.csv", row.names = FALSE)
