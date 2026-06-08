if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}


suppressPackageStartupMessages({
  library(ranger)
})

# Example: use saved antibiotic random forest models on a CDS detection matrix.
#
# Required inputs:
# 1) antibiotic_random_forest_models.rds
#    A named list of fitted ranger models, one per antibiotic.
# 2) antibiotic_random_forest_input_feature_ids.csv
#    A table listing which CDS IDs were used for each antibiotic model.
# 3) A detection matrix with:
#    - rows = CDS IDs
#    - columns = sample IDs
#    - values = 0/1 for absence/presence
#
# Example detection matrix:
# validation_detection_matrix.rds

models_file <- "antibiotic_random_forest_models.rds"
feature_ids_file <- "antibiotic_random_forest_input_feature_ids.csv"
detection_matrix_file <- "validation_detection_matrix.rds"

models_list <- readRDS(models_file)
feature_ids_df <- read.csv(feature_ids_file, check.names = FALSE, stringsAsFactors = FALSE)
detection_mat <- readRDS(detection_matrix_file)
dim(detection_mat)

# Pick one antibiotic model to use.
antibiotic_name <- "tobramycin"

if (!antibiotic_name %in% names(models_list)) {
  stop("Antibiotic model not found: ", antibiotic_name)
}

model_fit <- models_list[[antibiotic_name]]
feature_ids <- feature_ids_df$CDS[feature_ids_df$antibiotic == antibiotic_name]
feature_ids <- unique(feature_ids[!is.na(feature_ids) & feature_ids != ""])

if (length(feature_ids) == 0) {
  stop("No feature IDs found for antibiotic: ", antibiotic_name)
}

# Keep only features that exist in the new detection matrix.
common_features <- intersect(feature_ids, row.names(detection_mat))

if (length(common_features) == 0) {
  stop("No model features were found in the detection matrix for: ", antibiotic_name)
}

# Build the sample-by-feature matrix expected by the model.
prediction_matrix <- t(as.matrix(detection_mat[common_features, , drop = FALSE]))
mode(prediction_matrix) <- "numeric"
prediction_df <- as.data.frame(prediction_matrix, stringsAsFactors = FALSE, check.names = FALSE)

# Add missing model features as 0 if they are absent from the new matrix.
missing_features <- setdiff(feature_ids, colnames(prediction_df))
if (length(missing_features) > 0) {
  for (feature_id in missing_features) {
    prediction_df[[feature_id]] <- 0
  }
}

# Reorder columns to match the model feature order.
prediction_df <- prediction_df[, feature_ids, drop = FALSE]

# Predict resistant probability.
prediction_out <- predict(model_fit, data = prediction_df)$predictions

result_df <- data.frame(
  sample_id = row.names(prediction_df),
  non_resistant_probability = prediction_out[, "non_resistant"],
  resistant_probability = prediction_out[, "resistant"],
  predicted_class = ifelse(
    prediction_out[, "resistant"] >= 0.5,
    "resistant",
    "non_resistant"
  ),
  stringsAsFactors = FALSE
)

print(head(result_df, 10))

write.csv(
   result_df,
   paste0(antibiotic_name, "_prediction_example.csv"),
   row.names = FALSE
)
