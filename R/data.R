#' Annotated phylogenetic tree for pL1000_200_1
#'
#' A rooted binary \code{phylo} object (200 tips, 199 internal nodes) with
#' per-node confidence scores stored in \code{node.label}.
#' Built from a small in vitro LARRY barcode sample (pL1000, 200 cells).
#'
#' @format An \code{ape::phylo} object.
#' @source MitoDrift inference on the pL1000_200_1 dataset.
"pL1000_tree_annot"

#' Mutation count data for pL1000_200_1
#'
#' Long-format data frame of mitochondrial variant read counts
#' (200 cells, 186 variants, 37 200 rows).
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{cell}{Cell barcode identifier.}
#'   \item{variant}{Variant identifier (e.g. \code{"10111_T_C"}).}
#'   \item{a}{Alternate allele read count.}
#'   \item{d}{Total read depth.}
#' }
#' @source pL1000_200_1 dataset.
"pL1000_mut_dat"
