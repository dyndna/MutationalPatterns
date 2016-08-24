#' Extract mutational signatures from 96 mutation matrix using NMF
#' 
#' Decomposes trinucleotide count matrix into signatures and contribution of those signatures to the spectra of the samples/vcf files
#' @param mut_matrix 96 mutation count matrix 
#' @param rank Number of signatures to extract
#' @param nrun Number of iterations, default = 200
#' @return Named list of mutation matrix, signatures and signature contribution
#' @importFrom NMF nmf
#' @importFrom NMF basis
#' @importFrom NMF coef
#'
#' @examples
#' # Make 96 trinucleodide mutation count matrix
#' tri_matrix = mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
#' # Extract the signatures
#' nmf_res = extract_signatures(tri_matrix, rank = 3)
#' # Provide signature names (optional)
#' colnames(nmf_res$signatures) = c("Signature A", "Signature B")
#' # Plot signatures
#' plot_96_profile(nmf_res$signatures)
#'
#' @seealso \code{\link{mut_matrix}}
#'          \code{\link{plot_96_profile}}
#'
#' @export

extract_signatures = function(mut_matrix, rank, nrun = 200)
{
  mut_matrix = as.matrix(mut_matrix)
  # add small pseudocount to avoid features with zero counts
  mut_matrix = mut_matrix + 0.0001
  # Check if rank_range is appropriate
  if (!(rank > 0 & rank == round(rank)))
    stop("Rank should be a positive integer")
  if (ncol(mut_matrix) < max(rank))
    stop("Rank should be smaller than the number of samples in the input matrix")
  # Calculate nmf
  print("Decomposing matrix using NMF...")
  res = nmf(mut_matrix, rank = rank, method = "brunet", nrun=nrun, seed = 123456)
  print(paste("Number of iterations:", nrun))
  # Find signatures and contribution of signatures
  signatures = basis(res)
  contribution = coef(res)
  # Reconstruct mutation matrix
  reconstructed = signatures %*% contribution
  return(list(signatures = signatures, contribution = contribution, reconstructed = reconstructed))
}
