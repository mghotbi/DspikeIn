# CLR (Centered Log-Ratio): Requires pseudocount addition to handle zeros.
# ALR (Additive Log-Ratio): Also requires pseudocount addition and a reference taxon.
# Hellinger: Square root of relative abundances.
# Log10: Requires pseudocount addition to handle zeros.
# Logp1: Logarithm with base 10 plus 1.
# Relative Abundance: Proportions of total counts.
# Z-transform: Standardizes to zero mean and unit variance.

normalize_phyloseq <- function(physeq, method = 'clr', pseudocount = 1) {
  # Convert OTU table to matrix
  otu_table <- as(otu_table(physeq), "matrix")
  
  # Normalize OTU table based on the selected method
  normalized_otu <- switch(method,
                           'clr' = {
                             # CLR transformation from compositions package
                             otu_table <- otu_table + pseudocount
                             compositions::clr(otu_table)
                           },
                           'alr' = {
                             # ALR transformation from compositions package
                             otu_table <- otu_table + pseudocount
                             compositions::alr(otu_table)
                           },
                           'hellinger' = {
                             # Hellinger transformation from vegan package
                             vegan::decostand(otu_table, method = 'hellinger')
                           },
                           'log10' = {
                             # Log10 transformation with pseudocount (base R)
                             log10(otu_table + pseudocount)
                           },
                           'logp1' = {
                             # Natural log transformation with pseudocount (base R)
                             log(otu_table + pseudocount)
                           },
                           'relabundance' = {
                             # Relative abundance (proportion) transformation (manual calculation)
                             sweep(otu_table, 2, colSums(otu_table), FUN = "/")
                           },
                           'Z' = {
                             # Z-transform from microbiome package
                             microbiome::transform(physeq, 'Z')
                           },
                           stop("Unknown method")
  )
  
  # If Z-transform method is used, return the normalized phyloseq object directly
  if (method == 'Z') {
    return(normalized_otu)  # This method already returns a phyloseq object
  } else {
    # Update the OTU table in the original phyloseq object with the normalized data
    otu_table(physeq) <- otu_table(normalized_otu, taxa_are_rows(physeq))
    return(physeq)
  }
}

# Example usage
#ps <- physeq_ITS_adj_scaled_Absolute_OTU

# Normalize the phyloseq object using different methods
#ps_normalized_clr <- normalize_phyloseq(ps, method = 'clr')
#ps_normalized_alr <- normalize_phyloseq(ps, method = 'alr')
#ps_normalized_hellinger <- normalize_phyloseq(ps, method = 'hellinger')
#ps_normalized_log10 <- normalize_phyloseq(ps, method = 'log10')
#ps_normalized_logp1 <- normalize_phyloseq(ps, method = 'logp1')
#ps_normalized_relabundance <- normalize_phyloseq(ps, method = 'relabundance')
#ps_normalized_Z <- normalize_phyloseq(ps, method = 'Z')
