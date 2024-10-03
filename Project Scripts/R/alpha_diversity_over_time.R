### Alpha diversity over time

library(mia)

## Summer
physeq_summer_body <- qza_to_phyloseq(
  features="SUB-dada-filtered-nmnc-table.qza",
  tree="SUB-dada-rooted-tree.qza",
  "SUB-taxonomy.qza",
  metadata = "SUMMERBODYMETAWITHSAMPLEIDCOLUMN.tsv"
)

summer_body_TSE_justshannon <- makeTreeSummarizedExperimentFromPhyloseq(physeq_summer_body)

summer_body_TSE_justshannon <- mia::estimateDiversity(summer_body_TSE_justshannon,
                                                      assay.type = "counts",
                                                      index = "shannon",
                                                      name = "shannon")

summer_body_alpha_plot_data_justshannon <- data.frame(
  Carcass.id = colData(summer_body_TSE_justshannon)$Carcass.id,
  Stage = colData(summer_body_TSE_justshannon)$Stage,
  Shannon = colData(summer_body_TSE_justshannon)$shannon,
  ADH = colData(summer_body_TSE_justshannon)$ADH
)

summer_body_alpha_plot_data_justshannon$ADH <- as.numeric(as.character(summer_body_alpha_plot_data_justshannon$ADH))


shannon_summary <- summer_body_alpha_plot_data_justshannon %>%
  group_by(ADH) %>%
  summarise(mean_shannon = mean(Shannon),
            se_shannon = sd(Shannon) / sqrt(n()))

summeralphaovertime <- ggplot(data = shannon_summary, aes(x = ADH, y = mean_shannon)) +
  annotate("rect", xmin=0, xmax=150, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#fde724")  +
  annotate("rect", xmin=150, xmax=510, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#5cc862")  +
  annotate("rect", xmin=510, xmax=1050, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#20908c")  +
  annotate("rect", xmin=1050, xmax=1380, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#3b528b")  +
  annotate("rect", xmin=1380, xmax=2000, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#440154")  +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(data = shannon_summary, aes(y = mean_shannon, ymin = mean_shannon - se_shannon, ymax = mean_shannon + se_shannon), width = 0.1, size = 1) +
  labs(x = "ADH", y = "Mean Alpha Diversity (Shannon)", title = "Summer 2022", color = "Stage") +
  theme_pub() +
  theme(
    axis.title = element_text(size=14, face = "bold"),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "Econ Sans Cnd"),
    plot.title = element_text(hjust = 0, size = 18, face = "bold")  # Add this line
  ) +
  scale_x_continuous(breaks = c(0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000), expand = c(0,0)) +
  ylim(0, 4)

summeralphaovertime

### Linear regressions

summer_body_alpha_plot_data_justshannon$ADH <- as.numeric(as.character(summer_body_alpha_plot_data_justshannon$ADH))

lm_model_summer_shannon <- lm(Shannon ~ ADH, data = summer_body_alpha_plot_data_justshannon)
summary(lm_model_summer_shannon)

regressionsummeralphabodyovertime <- ggplot(data = summer_body_alpha_plot_data_justshannon, aes(x = ADH, y = Shannon)) +
  annotate("rect", xmin=0, xmax=150, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#fde724")  +
  annotate("rect", xmin=150, xmax=510, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#5cc862")  +
  annotate("rect", xmin=510, xmax=1050, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#20908c")  +
  annotate("rect", xmin=1050, xmax=1380, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#3b528b")  +
  annotate("rect", xmin=1380, xmax=2000, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#440154")  +
  geom_point() +  # Scatter plot of the data
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +  # Fitted regression line
  labs(x = "ADH", y = "Alpha Diversity (Shannon)") +
  theme_pub() +
  theme(text = element_text(family = "Econ Sans Cnd")) +
  theme(axis.title = element_text(size=14, face = "bold")) +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = c(0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000), expand = c(0,0)) +
  ylim(0, 4)

regressionsummeralphabodyovertime

summer_alpha_over_time <- (summeralphaovertime / regressionsummeralphabodyovertime)

summer_alpha_over_time


### Fall
physeq_fall_body <- qza_to_phyloseq(
  features="FAB-dada-filtered-nmnc-table.qza",
  tree="FAB-dada-rooted-tree.qza",
  "FAB-taxonomy.qza",
  metadata = "FABMETA-STAGES.tsv"
)

# Convert the phyloseq object to Tree Summarized Experiment (TSE) object for future use
fall_body_TSE_justshannon <- makeTreeSummarizedExperimentFromPhyloseq(physeq_fall_body)

# Generating our estimate for shannon diversity
fall_body_TSE_justshannon <- mia::estimateDiversity(fall_body_TSE_justshannon,
                                                    assay.type = "counts",
                                                    index = "shannon",
                                                    name = "shannon")

fall_body_alpha_plot_data_justshannon <- data.frame(
  Carcass.id = colData(fall_body_TSE_justshannon)$Carcass.id,
  Stage = colData(fall_body_TSE_justshannon)$Stage,
  Shannon = colData(fall_body_TSE_justshannon)$shannon,
  ADH = colData(fall_body_TSE_justshannon)$ADH
)

fall_body_alpha_plot_data_justshannon$ADH <- as.numeric(as.character(fall_body_alpha_plot_data_justshannon$ADH))


shannon_summary_fall <- fall_body_alpha_plot_data_justshannon %>%
  group_by(ADH) %>%
  summarise(mean_shannon = mean(Shannon),
            se_shannon = sd(Shannon) / sqrt(n()))

fallalphaovertime <- ggplot(data = shannon_summary_fall, aes(x = ADH, y = mean_shannon)) +
  annotate("rect", xmin=0, xmax=70, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#fde724")  +
  annotate("rect", xmin=70, xmax=500, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#5cc862")  +
  annotate("rect", xmin=500, xmax=1100, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#20908c")  +
  annotate("rect", xmin=1100, xmax=1550, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#3b528b")  +
  annotate("rect", xmin=1550, xmax=1750, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#440154")  +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(data = shannon_summary_fall, aes(y = mean_shannon, ymin = mean_shannon - se_shannon, ymax = mean_shannon + se_shannon), width = 0.1, size = 1) +
  labs(x = "ADH", y = "Mean Alpha Diversity (Shannon)", title = "Fall 2022", color = "Stage") +
  theme_pub() +
  theme(
    axis.title = element_text(size=14, face = "bold"),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "Econ Sans Cnd"),
    plot.title = element_text(hjust = 0, size = 18, face = "bold")  # Add this line
  ) +
  scale_x_continuous(breaks = c(0, 250, 500, 750, 1000, 1250, 1500, 1750), expand = c(0,0)) +
  ylim(0, 4)

fallalphaovertime


fall_body_alpha_plot_data_justshannon$ADH <- as.numeric(as.character(fall_body_alpha_plot_data_justshannon$ADH))

lm_model_fall_shannon <- lm(Shannon ~ ADH, data = fall_body_alpha_plot_data_justshannon)
summary(lm_model_fall_shannon)

regressionfallalphabodyovertime <- ggplot(data = fall_body_alpha_plot_data_justshannon, aes(x = ADH, y = Shannon)) +
  annotate("rect", xmin=0, xmax=70, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#fde724")  +
  annotate("rect", xmin=70, xmax=500, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#5cc862")  +
  annotate("rect", xmin=500, xmax=1100, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#20908c")  +
  annotate("rect", xmin=1100, xmax=1550, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#3b528b")  +
  annotate("rect", xmin=1550, xmax=1750, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#440154")  +
  geom_point() +  # Scatter plot of the data
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +  # Fitted regression line
  labs(x = "ADH", y = "Alpha Diversity (Shannon)") +
  theme_pub() +
  theme(text = element_text(family = "Econ Sans Cnd")) +
  theme(axis.title = element_text(size=14, face = "bold")) +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = c(0, 250, 500, 750, 1000, 1250, 1500, 1750), expand = c(0,0)) +
  ylim(0, 4)

regressionfallalphabodyovertime

fall_alpha_over_time <- (fallalphaovertime / regressionfallalphabodyovertime)

fall_alpha_over_time


### Spring
physeq_spring_body <-qza_to_phyloseq(
  features="SPB-dada-filtered-nmnc-table.qza",
  tree="SPB-dada-rooted-tree.qza",
  "SPB-taxonomy.qza",
  metadata = "SPRINGBODYMETA.tsv"
)

spring_body_TSE <- makeTreeSummarizedExperimentFromPhyloseq(physeq_spring_body)

spring_body_TSE_justshannon <- mia::estimateDiversity(spring_body_TSE,
                                                      assay.type = "counts",
                                                      index = "shannon",
                                                      name = "shannon")


spring_body_alpha_plot_data_justshannon <- data.frame(
  Shannon = colData(spring_body_TSE_justshannon)$shannon,
  ADH = colData(spring_body_TSE_justshannon)$ADH
)

spring_body_alpha_plot_data_justshannon$ADH <- as.numeric(as.character(spring_body_alpha_plot_data_justshannon$ADH))


shannon_summary_spring <- spring_body_alpha_plot_data_justshannon %>%
  group_by(ADH) %>%
  summarise(mean_shannon = mean(Shannon),
            se_shannon = sd(Shannon) / sqrt(n()))

springalphaovertime <- ggplot(data = shannon_summary_spring, aes(x = ADH, y = mean_shannon)) +
  annotate("rect", xmin=0, xmax=600, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#fde724")  +
  annotate("rect", xmin=600, xmax=720, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#5cc862")  +
  annotate("rect", xmin=720, xmax=1050, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#20908c")  +
  annotate("rect", xmin=1050, xmax=1190, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#3b528b")  +
  annotate("rect", xmin=1190, xmax=1250, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#440154")  +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(data = shannon_summary_spring, aes(y = mean_shannon, ymin = mean_shannon - se_shannon, ymax = mean_shannon + se_shannon), width = 0.1, size = 1) +
  labs(x = "ADH", y = "Mean Alpha Diversity (Shannon)", title = "Spring 2023", color = "Stage") +
  theme_pub() +
  theme(
    axis.title = element_text(size=14, face = "bold"),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "Econ Sans Cnd"),
    plot.title = element_text(hjust = 0, size = 18, face = "bold")  # Add this line
  ) +
  scale_x_continuous(breaks = c(0, 250, 500, 750, 1000, 1250), expand = c(0,0)) +
  ylim(0, 4)

springalphaovertime


spring_body_alpha_plot_data_justshannon$ADH <- as.numeric(as.character(spring_body_alpha_plot_data_justshannon$ADH))

lm_model_spring_shannon <- lm(Shannon ~ ADH, data = spring_body_alpha_plot_data_justshannon)
summary(lm_model_spring_shannon)

regressionspringalphabodyovertime <- ggplot(data = spring_body_alpha_plot_data_justshannon, aes(x = ADH, y = Shannon)) +
  annotate("rect", xmin=0, xmax=600, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#fde724")  +
  annotate("rect", xmin=600, xmax=720, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#5cc862")  +
  annotate("rect", xmin=720, xmax=1050, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#20908c")  +
  annotate("rect", xmin=1050, xmax=1190, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#3b528b")  +
  annotate("rect", xmin=1190, xmax=1250, ymin=-Inf, ymax=Inf, alpha=0.1, fill="#440154")  +
  geom_point() +  # Scatter plot of the data
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +  # Fitted regression line
  labs(x = "ADH", y = "Alpha Diversity (Shannon)") +
  theme_pub() +
  theme(text = element_text(family = "Econ Sans Cnd")) +
  theme(axis.title = element_text(size=14, face = "bold")) +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = c(0, 250, 500, 750, 1000, 1250), expand = c(0,0)) +
  ylim(0, 4)

regressionspringalphabodyovertime

spring_alpha_over_time <- (springalphaovertime / regressionspringalphabodyovertime)

spring_alpha_over_time


all_alpha_over_time <- (spring_alpha_over_time | summer_alpha_over_time | fall_alpha_over_time)

all_alpha_over_time

ggsave("all_alpha_over_time.png", width = 44, height = 24, units = "cm", dpi = 300, all_alpha_over_time)

