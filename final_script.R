# ==============================================

# Floral Scent Prediction: AdaSampling PU Learning Evaluation

# ==============================================

################# Set up Java for the `rcdk` package #############
# 1. Download Java here https://adoptium.net/download?link=https%3A%2F%2Fgithub.com%2Fadoptium%2Ftemurin25-binaries%2Freleases%2Fdownload%2Fjdk-25.0.2%252B10%2FOpenJDK25U-jdk_aarch64_mac_hotspot_25.0.2_10.pkg&vendor=Adoptium
# 2. Type /usr/libexec/java_home in your terminal
# 3. Run this script with your output as JAVA_HOME:

Sys.setenv(JAVA_HOME = "/Library/Java/JavaVirtualMachines/temurin-25.jdk/Contents/Home")
#install.packages("rJava")
library(rJava)
.jinit() # This initializes the Java Virtual Machine

################# Load Other Packages #############

#install.packages(c("AdaSampling", "fingerprint", "e1071", "caret", 
#"ggplot2", "dplyr", "rcdk", "pROC"))

library(AdaSampling)
library(fingerprint)
library(e1071)
library(caret)
library(ggplot2)
library(dplyr)
library(rcdk)
library(pROC)

################# Load the data #############

# Load from github directly
url <- "https://raw.githubusercontent.com/ARY2260/openpom/main/openpom/data/curated_datasets/curated_GS_LF_merged_4983.csv"
odor_data <- read.csv(url)

head(odor_data)

# The 'floral' scent will be the target
target_scent <- "floral"
(true_labels <- odor_data[[target_scent]])


################# Get the features from the SMILES strings #############


# Extract the correct column and ensure it is formatted as text
smiles_strings <- as.character(odor_data$nonStereoSMILES)

# Parse the SMILES strings into molecule objects
mols <- parse.smiles(smiles_strings)

# Generate 'Circular' (ECFP6) fingerprints - standard for scent prediction
fps <- lapply(mols, get.fingerprint, type='circular')

# Convert the list of fingerprints into a binary matrix for AdaSampling
# Each column is a specific chemical substructure bit
feature_matrix <- fp.to.matrix(fps)

################# Positive Unlabled Learning ##########################


# 5. Create the PU (Positive-Unlabeled) Setup
set.seed(42)
pos_indices <- which(true_labels == 1)

# Change 30% of known floral molecules from 1 to 0 (this is the "unlabeled" set)
hidden_indices <- sample(pos_indices, length(pos_indices) * 0.3)
pu_labels <- true_labels
pu_labels[hidden_indices] <- 0  

# AdaSampling requires named rows to track the molecules
rownames(feature_matrix) <- paste0("mol_", 1:nrow(feature_matrix))

# Separate the row names into Known Positives (1s) and Unlabeled (0s)
Ps_names <- rownames(feature_matrix)[which(pu_labels == 1)]
Ns_names <- rownames(feature_matrix)[which(pu_labels == 0)]

################# Run AdaSampling ##########################

# Explicitly name the columns
colnames(feature_matrix) <- paste0("Feature_", 1:ncol(feature_matrix))

# We pass the same matrix to train.mat and test.mat so we get scores for every molecule
svm_results <- adaSample(Ps = Ps_names, 
                         Ns = Ns_names, 
                         train.mat = feature_matrix, 
                         test.mat = feature_matrix, 
                         classifier = "svm")

knn_results <- adaSample(Ps = Ps_names, 
                         Ns = Ns_names, 
                         train.mat = feature_matrix, 
                         test.mat = feature_matrix, 
                         classifier = "knn")

logit_results <- adaSample(Ps = Ps_names, 
                         Ns = Ns_names, 
                         train.mat = feature_matrix, 
                         test.mat = feature_matrix, 
                         classifier = "logit")

################# Model Evaluation ##########################

svm_eval <- data.frame(
  Actual = as.factor(true_labels),
  PU_Input = pu_labels,
  Score = svm_results[,"P"] 
)

knn_eval <- data.frame(
  Actual = as.factor(true_labels),
  PU_Input = pu_labels,
  Score = knn_results[,"P"]
)

logit_eval <- data.frame(
  Actual = as.factor(true_labels),
  PU_Input = pu_labels,
  Score = logit_results[,"P"]
)

# Return a dataframe with 'P' (Positive) and 'N' (Negative) probabilities
head(svm_eval)
head(knn_eval)
head(logit_eval)

# Calculate all three ROC curves
# Note: Fixed 'eval_df' to 'svm_eval'
library(pROC)
svm_roc <- roc(response = svm_eval$Actual, predictor = svm_eval$Score, levels = c("0", "1"), direction = "<")
knn_roc <- roc(response = knn_eval$Actual, predictor = knn_eval$Score, levels = c("0", "1"), direction = "<")
logit_roc <- roc(response = logit_eval$Actual, predictor = logit_eval$Score, levels = c("0", "1"), direction = "<")

# Compare the models
auc_comparison <- data.frame(
  Model = c("SVM", "KNN", "Logistic Regression"),
  AUC = c(as.numeric(auc(svm_roc)), as.numeric(auc(knn_roc)), as.numeric(auc(logit_roc)))
)

print(auc_comparison)

################# Plot Results ##########################

# AUC plot
roc_plot <- ggroc(list(SVM = svm_roc, KNN = knn_roc, Logit = logit_roc), size = 1) +
  theme_minimal() +
  ggtitle("ROC Curves: AdaSampling Models for Floral Scent") +
  xlab("Specificity") + ylab("Sensitivity") +
  scale_color_manual(values = c("SVM" = "blue", "KNN" = "green", "Logit" = "red")) +
  theme(legend.title = element_blank())

print(roc_plot)

# Save plot
ggsave(filename = "ROC_Comparison_Plot.png", plot = roc_plot, width = 8, height = 6, dpi = 300)

# Distribution plot 
plot_distribution <- function(eval_data, model_name) {
  ggplot(eval_data, aes(x = Score, fill = Actual)) +
    geom_density(alpha = 0.6) +
    theme_minimal() +
    scale_fill_manual(values = c("0" = "gray60", "1" = "darkorchid"), 
                      labels = c("True Non-Floral", "True Floral")) +
    labs(title = paste(model_name, "- Probability Distribution"),
         x = "Predicted Probability of Being Floral",
         y = "Density",
         fill = "True Positive") +
    theme(legend.position = "bottom")
}

# Generate the individual plots
svm_dist <- plot_distribution(svm_eval, "SVM")
knn_dist <- plot_distribution(knn_eval, "KNN")
logit_dist <- plot_distribution(logit_eval, "Logistic Regression")

print(svm_dist)
print(knn_dist)
print(logit_dist)

# Save the distribution plots 
ggsave(filename = "SVM_Density_Distribution.png", plot = svm_dist, width = 8, height = 6, dpi = 300)
ggsave(filename = "KNN_Density_Distribution.png", plot = knn_dist, width = 8, height = 6, dpi = 300)
ggsave(filename = "Logit_Density_Distribution.png", plot = logit_dist, width = 8, height = 6, dpi = 300)


