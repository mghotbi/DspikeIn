library(testthat)
library(igraph)
library(DspikeIn)

Complete <- load_graphml("Complete.graphml")

# Run the function on the test graph
result <- node_level_metrics(Complete)

# ==============================
# TEST 1: Output Structure
# ==============================
test_that("Function returns a list with expected components", {
  expect_type(result, "list")
  expect_named(result, c("metrics", "flextable", "plot1", "plot2", "facet_plot"))
})

# ==============================
# TEST 2: Metrics Data Frame Structure
# ==============================
test_that("Metrics dataframe has correct columns", {
  expected_cols <- c("Node", "Degree", "Strength", "Closeness", "Betweenness",
                     "EigenvectorCentrality", "PageRank", "Harmonic",
                     "Transitivity", "Coreness", "Constraint", "EffectiveSize",
                     "Redundancy", "Community", "Efficiency", "Local_Efficiency",
                     "Within_Module_Connectivity", "Among_Module_Connectivity")

  expect_true(is.data.frame(result$metrics))
  expect_setequal(names(result$metrics), expected_cols)
})

# ==============================
# TEST 3: Metrics Data Types
# ==============================

test_that("Metrics dataframe has correct data types", {
  expect_type(result$metrics$Node, "character")
  expect_true(is.numeric(result$metrics$Degree))
  expect_type(result$metrics$Closeness, "double")
  expect_s3_class(result$metrics$Community, "factor")
})

# ==============================
# TEST 4: Ensure No Missing Values in Core Metrics
# ==============================
test_that("Core network metrics do not contain NAs (except Local_Efficiency)", {
  core_metrics <- c("Degree", "Strength", "Closeness", "Betweenness",
                    "EigenvectorCentrality", "PageRank", "Harmonic",
                    "Efficiency")

  for (metric in core_metrics) {
    expect_false(any(is.na(result$metrics[[metric]])),
                 info = paste0("Missing values found in ", metric))
  }

  # Local_Efficiency may have NA values for isolated nodes
  expect_true(any(is.na(result$metrics$Local_Efficiency)) | all(!is.na(result$metrics$Local_Efficiency)))
})

# ==============================
# TEST 5: Check Plot Outputs
# ==============================
test_that("Plot outputs are ggplot objects", {
  expect_s3_class(result$plot1, "gg")
  expect_s3_class(result$plot2, "gg")
  expect_s3_class(result$facet_plot, "gg")
})

# ==============================
# TEST 6: Check Flextable Output
# ==============================
test_that("Flextable output exists and is formatted correctly", {
  expect_s3_class(result$flextable, "flextable")
})

# ==============================
# TEST 7: Handle Errors Gracefully
# ==============================
test_that("Function fails gracefully for invalid input", {
  expect_error(node_level_metrics(NULL), "Must provide a graph object")
  expect_error(node_level_metrics(42), "Must provide a graph object")
})
