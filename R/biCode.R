# biCode.R
#
# Make a 5 character code from a binomial name by concatenating
# the first three letter of the first word and the first
# two letters of the second word.
#
# There are many ways to do this; here we assemble the two parts
# in a loop over the input, this way the function is vectorized and can
# work on a vector of input names.
#
# Parameters:
#    s (character) Binomial species name(s)
#
# Return value: (character) a vector of five-character uppercase
#                           codes, one per line of input.
#
# Details:
#    The function first writes the string "xxxxx" - these letters are over-
#    written with the actual code characters in uppercasder IF those exist.
#
# Author:
#   Boris Steipe (https://orcid.org/0000-0002-1134-6758)
#
# Examples:
#    biCode("Saccharomyces cerevisiae")
#    biCode("Saccharomycetales")
#    biCode("Sa c")
#    biCode("Sa ")
#    biCode("")
#    biCode("Schizosaccharomyces pombe 972h-")
#    biCode(c("Bipolaris oryzae", "Ustilago maydis", "Aspergillus nidulans"))
#

biCode <- function(s) {
  b <- character()
  for (i in 1:length(s)) {
    b[i] <- "xxxxx"
    bc <- strsplit(s[i], "\\s+")[[1]]
    if (length(bc) > 1) {
      substr(b[i], 1, 3) <- toupper(substr(bc[1], 1, 3))
      substr(b[i], 4, 5) <- toupper(substr(bc[2], 1, 2))
    } else if (length(bc) > 0) {
      substr(b[i], 1, 3) <- toupper(substr(bc[1], 1, 3))
    }
  }
  return(b)
}


# [END]
