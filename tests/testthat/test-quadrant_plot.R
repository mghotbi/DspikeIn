testthat::skip_on_bioc()

test_that("quadrant_plot generates a valid ggplot object using Complete.graphml", {
  library(testthat)
  library(igraph)
  library(ggplot2)
  library(DspikeIn)
  
  g <- load_graphml("Complete.graphml")
  suppressWarnings({
    result <- node_level_metrics(g)
  })
  #  Compute node-level metrics
  result <- node_level_metrics(g)
  metrics <- result$metrics # Extract metrics dataframe

  # Ensure the function returns a ggplot object
  plot <- quadrant_plot(metrics, x_metric = "Degree", y_metric = "Efficiency")
  expect_s3_class(plot, "ggplot")

  expect_true("Degree" %in% names(metrics))
  expect_true("Efficiency" %in% names(metrics))

  #  Test error handling for missing metric cols
  expect_error(quadrant_plot(metrics, x_metric = "InvalidMetric", y_metric = "Efficiency"),
    regexp = "Error: Specified x_metric or y_metric not found in metrics data frame."
  )

  expect_error(quadrant_plot(metrics, x_metric = "Degree", y_metric = "InvalidMetric"),
    regexp = "Error: Specified x_metric or y_metric not found in metrics data frame."
  )

  plot2 <- quadrant_plot(metrics, x_metric = "Degree", y_metric = "Redundancy")
  expect_s3_class(plot2, "ggplot")

  plot3 <- quadrant_plot(metrics, x_metric = "Among_Module_Connectivity", y_metric = "Betweenness")
  expect_s3_class(plot3, "ggplot")

  plot4 <- quadrant_plot(metrics, x_metric = "Within_Module_Connectivity", y_metric = "Degree")
  expect_s3_class(plot4, "ggplot")
})
