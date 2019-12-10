# H.R
#
# Compute Shannon entropy.
#
# Shannon entropy is computed in a discrete density distribution as the
# negative sum over all i of p_i * log(p_i). log(p_i) is computed as log2(p_i),
# thus the result has a unit of bits.
#
# Parameters:
#   pmf (numeric) probability mass function: a vector of states and
#                 associated probabilities, or counts. Each element in pmf
#                 must be > 0 and they are all scaled for sum(pmf) == 1. The
#                 return value of table() is therefore valid input.
#
# Return value:
#       (numeric) Shannon entropy in bits.
#
# Details:
#    None of the states of pmf must have a probability of 0, since this implies
#    that the distribution as a whole is impossible. Add pseudocounts if
#    non-observed states need to be included in the pmf. pmf may be given as raw
#    counts, since the elements are scaled to sum to 1.
#
# Author:
#   Boris Steipe (https://orcid.org/0000-0002-1134-6758)
#
# Examples:
#   H(c(A=0.25, C=0.25, G=0.25, T=0.25))  # 2 bits entropy in a random
#                                         # nucleotide sequence
#   H(rep(1, 20))  # Random peptide: 4.322 bits of entropy
#   H(1)     # If all elements are the same, entropy is zero
#

H <- function(pmf) {
  if (any(pmf <= 0)) {
    stop("Input is not a vector of counts or probabilities.")
  }
  pmf <- pmf / sum(pmf)
  H <- -sum(pmf * (log(pmf) / log(2)))
  return(H)
}

# [END]
