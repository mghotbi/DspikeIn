library(testthat)
library(igraph)
library(ggplot2)
library(DspikeIn)

test_that("simulate_network_robustness handles negative weights correctly", {

  #  Load GraphML from DspikeIn
  graph <- load_graphml("Complete.graphml")

  #  Assign positive weights (betweenness should work)
  E(graph)$weight <- abs(runif(ecount(graph), 0.1, 1))
  result_betweenness <- simulate_network_robustness(graph = graph, steps = 10, removal_strategy = "betweenness")

  expect_s3_class(result_betweenness$results, "data.frame")
  expect_s3_class(result_betweenness$plot, "ggplot")

  E(graph)$weight <- runif(ecount(graph), -1, 1)

  expect_warning(
    simulate_network_robustness(graph, steps = 10, removal_strategy = "betweenness"),
    "Adjusting edge weights to be positive"
  )
})
