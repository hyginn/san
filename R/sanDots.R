# sanDots.R
#
# Plot a dotplot. Function is not called dotplot bercause sequinr has a
# dotplot() function.
#
# Description
#
# Parameters:
#   ...
#
# Return value:
#       ...
#
# Details:
#    ...
#
# Author:
#   Boris Steipe (https://orcid.org/0000-0002-1134-6758)
#
# Examples:
#   ...
#

sanDots <- function(A, B,        # sequences
                    MDM = BLOSUM62,
                    f,           # filter
                    palette,     # a function that returns color values
                    xlab = "",
                    ylab = ""
) {
  # Purpose:
  #     Create a dotplot to measure sequence similarity between
  #     two amino acid sequences
  # Version:  1.0
  # Date:     2016-09
  # Author:   Boris Steipe
  #
  # Parameters:
  #     A, B: strings that contain no letters that are not found in MDM
  #     MDM: Mutation Data Matrix. Deafaults to BLOSUM62 which must be
  #          defined. Load package biostrings from BioConductor and load
  #          the matrix with data(BLOSUM62).
  #     f: filter matrix to weight an average around the neighborhood of
  #        an amino acid pair. Default to the identity matrix if missing.
  #        Average over a window of length f if length(f) is 1.
  #     palette: rainbow(), cm.colors(), or another function that returns
  #              a palette of color hexcodes. If missing, make our own
  #              palette.
  # Value:
  #     none. creates a dotplot.


  if (missing(f)) {
    f <- matrix(1) # default
  } else if (length(f) == 1) {
    if (! f %% 2) {stop("Sorry: f must be odd.")}
    w <- f
    f <- matrix(numeric(w * w), nrow = w)
    for (i in 1:w) { f[i, i] <- 1 }  # identity matrix
  }

  if (missing(palette)) {
    palette <- colorRampPalette(c("#000000",  #black
                                  "#111111",
                                  "#222222",
                                  "#332222",  # grey
                                  "#CC7755",
                                  "#EE8844",
                                  "#FF4400"), # red
                                bias = 0.8)
  }
  A <- unlist(strsplit(A, ""))
  V <- unlist(strsplit(B, ""))
  lA <- length(A)
  lB <- length(B)

  m <- matrix(numeric(lA * lB), nrow = lA, ncol = lB)
  for (i in 1:lA) {
    for (j in 1:lB) {
      m[i, j] <- MDM[A[i], B[j]]
    }
  }
  m2 <- m
  wr <- floor((dim(f)[1] - 1) / 2)  # half-window size for rows
  wc <- floor((dim(f)[2] - 1) / 2)  # half-window size for columns

  for (i in (wr + 1):(lA - wr)) {
    for (j in (wc + 1):(lB - wc)) {
      # apply the filter to each value in m by weighting and summing
      # over its wr x wc neighborhood. Put the new value in m2
      m2[i, j] <- sum(f * m[(i-wr):(i+wr), (j-wc):(j+wc)])
    }
  }
  image(1:lA, 1:lB, m2,
        col = palette(24),
        ylim=c(lB,1), xlim=c(1,lA),
        xlab = xlab,
        ylab = ylab,
        axes = FALSE)
  box()

  # find good values for axis ticks and gridlines
  steps <- c(1, 2, 5, 10, 20, 50, 100, 200, 500,
             1000, 2000, 5000, 10000, 20000, 50000)
  gridStep <- sum(steps < max(lA, lB))

  # draw axes
  axis(1, at = c(1, seq(steps[gridStep - 3], lA, by=steps[gridStep - 3])))
  axis(2, at = c(1, seq(steps[gridStep - 3], lB, by=steps[gridStep - 3])))
  axis(3, at = c(1, seq(steps[gridStep - 3], lA, by=steps[gridStep - 3])))
  axis(4, at = c(1, seq(steps[gridStep - 2], lB, by=steps[gridStep - 3])))

  # draw grid with thin, transparent lines
  for (pos in seq(steps[gridStep - 2], lA, by = steps[gridStep - 2])) {
    abline(v=pos, col = "#FFFFFF44", lwd = 0.5)
  }
  for (pos in seq(steps[gridStep - 2], lB, by = steps[gridStep - 2])) {
    abline(h=pos, col = "#FFFFFF44", lwd = 0.5)
  }
}


# [END]
