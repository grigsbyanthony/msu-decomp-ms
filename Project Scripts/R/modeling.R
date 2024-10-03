### Modeling attempts are documented here

# Libraries
library(ggplot2)
library(dplyr)
library(lubridate)
library(gridExtra)
library(phyloseq)
library(qiime2r)
library(randomForest)
library(caret)
library(ggdendro)
library(gridExtra)
library(dendextend)
library(ggtree)
library(ape)
library(patchwork)
library(ggExtra)

### Carcass modeling

### Making classifier with seasonality
all_carcasses <- qza_to_phyloseq(
  features="ALLBODY-dada-filtered-nmnc-table.qza",
  tree="ALLBODY-dada-rooted-tree.qza",
  "taxonomy.qza",
  metadata = "ALLSEASONSNOTEMPS.tsv"
)

otu_table_all_carcasses <- as.data.frame(otu_table(all_carcasses))
all_carcasses_sample_data <- as.data.frame(sample_data(all_carcasses))

otu_table_all_carcasses_t <- t(otu_table_all_carcasses)

# Combine
all_carcasses_data <- cbind(all_carcasses_sample_data, otu_table_all_carcasses_t)

# Extract and format response data
response <- all_carcasses_data$Season
response <- as.factor(response)

features <- all_carcasses_data[, !names(all_carcasses_data) %in% c("ADH", "SAMPLE", "Carcass.id", "date", "sampling.event",
                                                                             "Season", "Sampletype", "ADH.numeric", "Stage", "Larvae",
                                                                             "Temperature", "Humidity")]

# Model parameters
set.seed(123)
trainIndex <- createDataPartition(response, p = 0.8, list = FALSE)
trainData <- features[trainIndex, ]
trainResponse <- response[trainIndex]
testData <- features[-trainIndex, ]
testResponse <- response[-trainIndex]

carcass_seasonality_model <- randomForest(x = trainData, y = trainResponse, ntree = 500, importance = TRUE)

print(carcass_seasonality_model)

# Store predictions and generate a preliminary matrix
predictions <- predict(carcass_seasonality_model, newdata = testData)

confusionMatrix <- table(Predicted = predictions, Actual = testResponse)

accuracy <- sum(diag(confusionMatrix)) / sum(confusionMatrix)

print(confusionMatrix)

cat("Accuracy:", accuracy, "\n")

# Look at variable importance
importance <- importance(carcass_seasonality_model)

var_importance <- data.frame(Variable = rownames(importance), Importance = importance[, 1])

var_importance <- var_importance[order(var_importance$Importance, decreasing = TRUE), ]

# Filtering to the 10 most important features
top_10_importance <- head(var_importance, 10)
top_10_importance

# Map the table and merge
top_10_importance_otus <- data.frame(
  OTU = c("b3d2cc3b9f0184830bec2ef73bcc6166", "4850a238b7fbbd5ffa204c5a5a396f2c",
          "4d72007c70f4abc0bf9e8fe8b826d07d", "8d0bca9f657ef580a9e8b864f89a99bf",
          "d15bc449222795a9ff230013aa633686", "04524511556ae43c41266c364d115459",
          "4005f039311a1265dd8f02e81c3bbd8f", "cd64f9b653d0444b04264fd6f09b715e",
          "945184b6386c192c0066e0a98a154780", "9580aa8a52ec18a4e92e701cdb595faa"),
  Genus = c("Lactococcus garvieae", "Hafnia-Obesumbacterium sp.", "Enterobacteriaceae sp 1.", "Psychrobacter sp.", "Acinetobacter sp.",
            "Ignatzschineria sp. 1", "Enterococcus sp.", "Streptococcus sp.", "Enterobacteriaceae sp 2.", "Clostridium perfringens")
)

top_10_importance_with_genus <- merge(top_10_importance, top_10_importance_otus, by.x = "Variable", by.y = "OTU")
print(top_10_importance_with_genus)

# Plot importance
feature_importance_seasonality <- ggplot(top_10_importance_with_genus, aes(x = reorder(Genus, Importance), y = Importance, fill = Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  ylab("Importance") +
  theme_pub() +
  scale_x_discrete(expand = c(0, 0)) +  # Remove space around x-axis
  scale_y_continuous(expand = c(0, 0), limits = c(0, 25)) +  # Remove space around y-axis
  scale_fill_gradient(low = "#5ec962", high = "darkgreen") +  # Gradient from light to dark green
  theme(text = element_text(family = "Econ Sans Cnd"),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_blank(),
        legend.position = "none")   # Remove legend

feature_importance_seasonality

# with tree
tree <- phy_tree(all_carcasses)
top_10_otus <- top_10_importance_with_genus$Variable
pruned_tree <- ape::keep.tip(tree, top_10_otus)
tree_order <- pruned_tree$tip.label

top_10_importance_with_genus <- top_10_importance_with_genus %>%
  mutate(Variable = factor(Variable, levels = tree_order))

tree_plot <- ggtree(pruned_tree) +
  theme_tree2(family = "Econ Sans Cnd") +
  theme(legend.position = "none") +
  geom_tiplab(align = TRUE, size = 0)  # Adjust size as needed

tree_plot

feature_importance_seasonality <- ggplot(top_10_importance_with_genus, aes(x = Variable, y = Importance, fill = Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  ylab("Importance") +
  theme_pub() +
  scale_x_discrete(expand = c(0, 0), labels = top_10_importance_with_genus$Genus) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 25)) +
  scale_fill_gradient(low = "#5ec962", high = "darkgreen") +
  theme(
    text = element_text(family = "Econ Sans Cnd"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold", hjust = 0.5),  # Center vertically
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_blank(),
    legend.position = "none",
  ) +
  geom_text(aes(label = Genus, y = 0), hjust = 1, nudge_y = -1, size = 4)  # Add genus labels

feature_importance_seasonality

combined_plot_seasonality <- tree_plot + feature_importance_seasonality +
  plot_layout(widths = c(1, 1)) &
  theme(plot.margin = margin(5.5, 5.5, 5.5, 0))

combined_plot_seasonality

# Making a nicer confusion matrix + class.error plot w/ ggplot
confusion_matrix <- carcass_seasonality_model$confusion
class_error <- confusion_matrix[, "class.error"]
conf_matrix_df <- as.data.frame(as.table(confusion_matrix[, -ncol(confusion_matrix)]))
names(conf_matrix_df) <- c("Reference", "Prediction", "Freq")

# Add class error to the data frame
conf_matrix_df$ClassError <- class_error[conf_matrix_df$Reference]

seasonality_model_confusion_matrix <- ggplot(data = conf_matrix_df, aes(x = Prediction, y = Reference, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%d", Freq)), vjust = 1, family = "Econ Sans Cnd", fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#5ec962") +
  xlab("Predicted Season") +
  ylab("Actual Season") +
  theme_pub() +
  theme(aspect.ratio = 1,
        text = element_text(family = "Econ Sans Cnd"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 20, face = "bold"),
        legend.position = "none") +
  coord_cartesian(clip = "off") + # Allow drawing outside the plot area
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))

# Display the plot
print(seasonality_model_confusion_matrix)


# Making a nicer confusion matrix w/ ggplot
conf_matrix <- confusionMatrix(predictions, testResponse)
conf_matrix_df <- as.data.frame(as.table(conf_matrix$table))

seasonality_confusion_matrix <- ggplot(data = conf_matrix_df, aes(x = Prediction, y = Reference, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.0f", Freq)), vjust = 1, family = "Econ Sans Cnd", fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#5ec962") +
  xlab("Predicted Season") +
  ylab("Actual Season") +
  theme_pub() +
  theme(aspect.ratio = 1,
        text = element_text(family = "Econ Sans Cnd"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 20, face = "bold"),
        legend.position = "none") +  # Remove legend
  scale_x_discrete(expand = c(0, 0)) +  # Remove space around x-axis
  scale_y_discrete(expand = c(0, 0))    # Remove space around y-axis

seasonality_confusion_matrix

# Combine the matrix and importance into one figure
carcass_seasonality_model_summary <- (seasonality_model_confusion_matrix | seasonality_confusion_matrix) / combined_plot_seasonality

print(carcass_seasonality_model_summary)

ggsave("carcass_seasonality_model_summary.png", width = 23, height = 23, units = "cm", dpi = 300, carcass_seasonality_model_summary)


### Making a regressor with ADH
response <- all_carcasses_data$ADH.numeric

# Prepare data for random forest
features <- all_carcasses_data[, !names(all_carcasses_data) %in% c("ADH", "SAMPLE", "Carcass.id", "date", "sampling.event", "Season", "Sampletype", "ADH.numeric", "Stage", "Larvae", "Temperature", "Humidity")]

# Model parameters
set.seed(123)
trainIndex <- createDataPartition(response, p = 0.8, list = FALSE)
trainData <- features[trainIndex, ]
trainResponse <- response[trainIndex]
testData <- features[-trainIndex, ]
testResponse <- response[-trainIndex]

carcass_adh_model <- randomForest(x = trainData, y = trainResponse, ntree = 500, importance = TRUE)
print(carcass_adh_model)

# Store predictions
predictions <- predict(carcass_adh_model, newdata = testData)

# Calculate performance metrics
mse <- mean((predictions - testResponse)^2)
rmse <- sqrt(mse)
r2 <- cor(predictions, testResponse)^2
mae <- mean(abs(predictions - testResponse))
rsq <- 1 - sum((testResponse - predictions)^2) / sum((testResponse - mean(testResponse))^2)

# Print performance metrics cutely
cat("Mean Squared Error:", mse, "\n")
cat("Root Mean Squared Error:", rmse, "\n")
cat("R-squared:", r2, "\n")

# Look at variable importance
importance <- importance(carcass_adh_model)
var_importance <- data.frame(Variable = rownames(importance), Importance = importance[, 1])
var_importance <- var_importance[order(var_importance$Importance, decreasing = TRUE), ]

# Filter to 10 most important variables
top_10_importance <- head(var_importance, 10)

top_10_importance_otus <- data.frame(
  OTU = c("c22b16cc6108c04f29fea3b6d4c81571", "d46e2205f0c6ecf67b51f83d111c509c",
          "04524511556ae43c41266c364d115459", "b3d2cc3b9f0184830bec2ef73bcc6166",
          "7339301f770131bd35e653eda2ac5106", "9908fffab7ed4f3bec44cda2f5084d49",
          "d12fa38dcafa4b438b14e45008e41f1a", "d114fb4c335125128be28401522dd41a",
          "69cb94bc7147e79c931c6c3fe7c5c4b1", "0e161099b2770fd8e02dca831d47f3bb"),
  Genus = c("Acinetobacter sp.", "Escherichia-Shigella sp.", "Ignatzschineria sp. 1", "Lactococcus garvieae", "Ignatzschineria sp. 2",
            "Enterococcus sp.", "Sporosarcina sp.", "Lactococcus lactis", "Wohlfahrtiimonas chitiniclastica", "Ignatzschineria sp. 3")
)

top_10_importance_with_genus <- merge(top_10_importance, top_10_importance_otus, by.x = "Variable", by.y = "OTU")
print(top_10_importance_with_genus)

# Plot importance
feature_importance_ADH <- ggplot(top_10_importance_with_genus, aes(x = reorder(Genus, Importance), y = Importance, fill = Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  ylab("Importance") +
  theme_pub() +
  scale_x_discrete(expand = c(0, 0)) +  # Remove space around x-axis
  scale_y_continuous(expand = c(0, 0), limits = c(0, 25)) +  # Remove space around y-axis
  scale_fill_gradient(low = "#5ec962", high = "darkgreen") +  # Gradient from light to dark green
  theme(text = element_text(family = "Econ Sans Cnd"),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_blank(),
        legend.position = "none")   # Remove legend

feature_importance_ADH

# with tree
tree <- phy_tree(all_carcasses)
top_10_otus <- top_10_importance_with_genus$Variable
pruned_tree <- ape::keep.tip(tree, top_10_otus)

# Get the order of tips in the tree
tree_order <- pruned_tree$tip.label

# Reorder your data frame to match the tree order
top_10_importance_with_genus <- top_10_importance_with_genus %>%
  mutate(Variable = factor(Variable, levels = tree_order))

# Create the tree plot
tree_plot <- ggtree(pruned_tree) +
  theme_tree2() +
  theme(legend.position = "none") +
  geom_tiplab(align = TRUE, size = 0)  # Adjust size as needed

tree_plot

feature_importance_ADH <- ggplot(top_10_importance_with_genus, aes(x = Variable, y = Importance, fill = Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  ylab("Importance") +
  theme_pub() +
  scale_x_discrete(expand = c(0, 0), labels = top_10_importance_with_genus$Genus) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 25)) +
  scale_fill_gradient(low = "#5ec962", high = "darkgreen") +
  theme(
    text = element_text(family = "Econ Sans Cnd"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold", hjust = 0.5),  # Center vertically
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_blank(),
    legend.position = "none",
  ) +
  geom_text(aes(label = Genus, y = 0), hjust = 1, nudge_y = -1, size = 4)  # Add genus labels


feature_importance_ADH

combined_importance_plot_ADH <- tree_plot + feature_importance_ADH +
  plot_layout(widths = c(1, 1)) &
  theme(plot.margin = margin(5.5, 5.5, 5.5, 0))

combined_importance_plot_ADH


# Preliminary model plotting
plot_data <- data.frame(Actual = testResponse, Predicted = predictions)

all_carcasses_adh_actual_versus_predicted <- ggplot(plot_data, aes(x = Actual, y = Predicted)) +
  geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlab("Actual ADH") +
  ylab("Predicted ADH") +
  theme_pub() +
  theme(aspect.ratio = 1,
    text = element_text(family = "Econ Sans Cnd"),
    axis.title = element_text(face = "bold", size = 20)
  ) +
  ylim(0, 2000) +  # Set y-axis limits
  xlim(0, 2000)

all_carcasses_adh_actual_versus_predicted

# Calculate residuals
residuals <- testResponse - predictions
residuals_data <- data.frame(Actual = testResponse, Residuals = residuals)
residual_density <- density(residuals_data$Residuals)

# Create the plot
all_carcasses_adh_residual_plot <- ggplot(residuals_data, aes(x = Actual, y = Residuals)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, color = "red") +
  xlab("Actual ADH") +
  ylab("Residuals") +
  ylim(-500, 500) +  # Set y-axis limits
  theme_pub() +
  theme(
    aspect.ratio = 1,
    text = element_text(family = "Econ Sans Cnd"),
    axis.title = element_text(face = "bold", size = 20))

print(all_carcasses_adh_residual_plot)

# Combine plots into one figure
carcass_adh_model_summary <- (all_carcasses_adh_actual_versus_predicted | all_carcasses_adh_residual_plot) / combined_importance_plot_ADH

print(carcass_adh_model_summary)

ggsave("carcass_adh_model_summary.png", width = 23, height = 23, units = "cm", dpi = 300, carcass_adh_model_summary)


