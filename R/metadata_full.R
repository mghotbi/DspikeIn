#' @title Metadata for Microbiome Samples
#'
#' @description
#' This dataset contains detailed metadata associated with microbiome sequencing samples.
#' It includes host information, sequencing metrics, environmental classifications, and
#' experimental factors relevant for downstream analyses.
#'
#' @format A data frame with 312 rows and 46 columns:
#' \describe{
#'   \item{sample.id}{Unique identifier for each sample.}
#'   \item{Total_Reads_total}{Total number of sequencing reads before processing.}
#'   \item{Total_Reads_spiked}{Total number of sequencing reads after spiking.}
#'   \item{Percentage}{Relative abundance percentage of the target taxa.}
#'   \item{Result}{Experimental result classification.}
#'   \item{beta.distances}{Beta diversity distances between samples.}
#'   \item{Observed}{Observed microbial richness in the sample.}
#'   \item{Chao1}{Chao1 diversity estimator.}
#'   \item{ACE}{Abundance-based Coverage Estimator (ACE).}
#'   \item{se.ACE}{Standard error of ACE.}
#'   \item{Shannon}{Shannon diversity index.}
#'   \item{Simpson}{Simpson diversity index.}
#'   \item{InvSimpson}{Inverse Simpson diversity index.}
#'   \item{X16S.biosample}{16S rRNA sequencing biosample identifier.}
#'   \item{dna.biosample}{DNA biosample identifier.}
#'   \item{data.type}{General data type classification.}
#'   \item{ampliconlibrary.quantification.ng.ul}{Library quantification (ng/\eqn{\mu}L).}
#'   \item{plate.ID}{Identifier for the sequencing plate.}
#'   \item{well.location}{Well location on the sequencing plate.}
#'   \item{Env.broad.scale}{Broad-scale environmental classification.}
#'   \item{Host.taxon}{Taxonomic classification of the host.}
#'   \item{Host.genus}{Genus of the host organism.}
#'   \item{Host.species}{Species of the host organism.}
#'   \item{Animal.type}{Classification of the host (e.g., Mammal, Bird).}
#'   \item{Animal.ecomode}{Ecological mode of the host organism.}
#'   \item{Clade.Order}{Taxonomic order of the host.}
#'   \item{Family}{Taxonomic family of the host.}
#'   \item{Diet}{General dietary classification of the host.}
#'   \item{Diet.Detailed}{Detailed dietary classification.}
#'   \item{Habitat}{Habitat description where the sample was collected.}
#'   \item{Metamorphosis}{Indicates whether the species undergoes metamorphosis.}
#'   \item{Reproduction}{Reproductive strategy of the host species.}
#'   \item{Ecoregion.III}{Ecoregion classification (level III).}
#'   \item{Ecoregion.IV}{Ecoregion classification (level IV).}
#'   \item{Site}{Sampling site identifier.}
#'   \item{sample.name}{Human-readable sample name.}
#'   \item{biosample.parent}{Parent biosample identifier.}
#'   \item{data.type.1}{Secondary data type classification.}
#'   \item{ampliconlibrary.quantification.ng.ul.1}{Second library quantification (ng/\eqn{\mu}L).}
#'   \item{plate.ID.1}{Secondary plate ID (if applicable).}
#'   \item{well.location.1}{Secondary well location (if applicable).}
#'   \item{sample.or.blank}{Indicates whether the sample is biological or a blank control.}
#'   \item{sample.spiked.blank}{Indicates whether the sample was spiked with a blank.}
#'   \item{spiked.volume}{Volume of spike added to the sample.}
#'   \item{swab.presence}{Indicates whether a swab was used for collection.}
#'   \item{MK.spike}{Molecular spike-in used in the sequencing process.}
#' }
#'
#' @return A data frame containing sample-level metadata for microbiome sequencing analysis.
#'
#' @usage data(metadata_full)
#' @keywords datasets
#'
#' @examples
#' data(metadata_full)
#' head(metadata_full)
#' summary(metadata_full)
"metadata_full"
