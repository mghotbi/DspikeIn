library(testthat)
library(phyloseq)
library(flextable)

test_that("conclusion function works correctly with embedded dataset", {

  # Define output file paths
  output_docx <- "merged_data.docx"
  output_csv <- "merged_data.csv"

  if (file.exists(output_docx)) file.remove(output_docx)
  if (file.exists(output_csv)) file.remove(output_csv)

  data("physeq_16SOTU", package = "DspikeIn")

  # Define test parameters
  merged_spiked_species <- c("Tetragenococcus_halophilus")
  max_passed_range <- 20

  result_physeq <- conclusion(physeq_16SOTU, merged_spiked_species, max_passed_range)

  expect_type(result_physeq, "list")

  expect_true(all(c("summary_stats", "full_report", "phy_tree") %in% names(result_physeq)))

  #  Validate output types
  expect_s3_class(result_physeq$summary_stats, "flextable")
  expect_s3_class(result_physeq$full_report, "data.frame")

  required_cols <- c("Sample", "Total_Reads", "Spiked_Reads", "Percentage", "Result")
  expect_true(all(required_cols %in% colnames(result_physeq$full_report)))

  # **Fix: Properly remove 'spiked.volume' column**
  physeq_invalid <- physeq_16SOTU
  metadata_invalid <- phyloseq::sample_data(physeq_invalid)  # Extract metadata
  metadata_invalid$spiked.volume <- NULL  # Remove required column
  phyloseq::sample_data(physeq_invalid) <- metadata_invalid  # Reassign modified metadata

  # Ensure function throws error when 'spiked.volume' is missing
  expect_error(
    conclusion(physeq_invalid, merged_spiked_species, max_passed_range),
    regexp = "'spiked.volume' column not found in metadata."
  )

  # Clean up files after running the test
  if (file.exists(output_docx)) file.remove(output_docx)
  if (file.exists(output_csv)) file.remove(output_csv)

})
