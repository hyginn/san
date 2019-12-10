# san.R
#
# Purpose: Sequence ANalysis. A simple introduction to methods of sequence
#          analysis and demonstration of a basic workflow.
# Version:  0.1
# Version history:
# Date:     2019-12-03
# Author:   Boris Steipe
# License:  MIT
#
# ToDo:
# Notes:
#
# ==============================================================================

# WARNING: SIDE EFFECTS
# Executing this script will execute code it contains.

# ====  PARAMETERS  ============================================================
# Execute the lines below to set important parameters for your work with this
# file.


# ====  PACKAGES  ==============================================================
# Load all required packages.
#

# === CRAN packages:

if (! requireNamespace("seqinr", quietly=TRUE)) {
  install.packages("seqinr")
}
# Package information:
#  library(help = seqinr)       # basic information
#  browseVignettes("seqinr")    # available vignettes
#  data(package = "seqinr")     # available datasets

if (! requireNamespace("httr", quietly=TRUE)) {
  install.packages("httr")
}



# === Bioconductor packages:

if (! requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}

if (! requireNamespace("Biostrings", quietly=TRUE)) {
  BiocManager::install("Biostringhs")
}
library(Biostrings)

if (! requireNamespace("msa", quietly=TRUE)) {
  BiocManager::install("msa")
}


# ====  FUNCTIONS  =============================================================

source("./R/biCode.R")             # make five-letter codes from binomial names
source("./R/H.R")                  # Compute Shannon entropy
source("./R/sanDots.R")            # plot a Dotplot

source("./scripts/sanCreateDB.R")  # create a database of reference sequences



# ====  PROCESS  ===============================================================


# === Sequence download

library(seqinr) # need library() for data()

URL <- sprintf("https://www.uniprot.org/uniprot/%s.fasta",
               "P39678")
tmp <- httr::GET(URL)
stopifnot(tmp$status_code == 200)
tmp <- gsub("^>[^ ]+",
            ">Mbp1",
            httr::content(tmp, as = "text", encoding = "UTF-8"))

tCon <- textConnection(tmp)
mbp1 <- seqinr::read.fasta(tCon, seqtype = "AA", as.string = TRUE)
close(tCon)

mbp1AA <- unlist(strsplit(mbp1$Mbp1, ""))

seqinr::computePI(mbp1AA)
seqinr::pmw(mbp1AA)

data(aaindex)  # "attach" the dataset

aaIndDesc <- character(length(aaindex))
for (i in seq_along(aaindex)) {
  aaIndDesc[i] <- sprintf("%d: %s", i, aaindex[[i]]$D)
}

aaIndDesc[grep("Composition", aaIndDesc)]

AAref <- aaindex[[459]]$I
names(AAref) <- a(names(AAref))
AAref <- AAref[order(names(AAref))]

AAobs <- (100 * table(mbp1AA)) / length(mbp1AA)

# color the bars by type.
# define colors

aaIndDesc[grep("Hydrophobicity", aaIndDesc)]

hyPho <- aaindex[[544]]$I
names(hyPho) <- a(names(hyPho))
hyPho <- hyPho[order(names(hyPho))]
hyPho <- round((255 * ((hyPho - min(hyPho)) / (max(hyPho) - min(hyPho)))) + 1)

if (! requireNamespace("viridis")) {
  install.packages("viridis")
}

hyCol <- viridisLite::viridis(n = length(hyPho))
names(hyCol) <- names(hyPho)[order(hyPho, decreasing = TRUE)]
barplot(rep(1, 20), names.arg = names(hyCol), col = hyCol, cex.names = 0.5)

lR <- sort(log(AAobs / AAref), decreasing = TRUE)
barplot(lR,
        ylim = c(-3.5,2),
        col = hyCol[names(lR)],
        xlab = "Amino acid",
        ylab = "log(Observed sequence / Database average)",
        cex.names = 0.9,
        main = "Amino acid enrichment in MBP1_YEAST")
abline(h = log(1), col="#00000055")
abline(h = log(c(1/2, 2)), lty = 2, col="#00000055")
abline(v = 9.7, lty = 2, lwd = 2, col="#6600BB")

arrows(4, 1.8, 0, 1.8, length = 0.07)
text(5.5, 1.8, "Enriched", cex = 0.9)
arrows(20, 1.8, 24, 1.8, length = 0.07)
text(19.5, 1.8, "Depleted", pos = 2, cex = 0.9)

# regex search for AT-hook

(m <- regexec("RGRP", mbp1[[1]][1]))
regmatches(mbp1[[1]][1], m)
sprintf("%s%s%s",
        tolower(substr(mbp1[[1]][1], 165, 167)),
                substr(mbp1[[1]][1], 168, 171),
        tolower(substr(mbp1[[1]][1], 172, 177)))

#SMART consensus http://smart.embl.de/smart/do_annotation.pl?DOMAIN=SM00384
# PRH_PETCR"] <-   KRSRGRPRKVQNS
# CONSENSUS/80%    t+tRGRP.+..t.
# CONSENSUS/65%    t+tRGRPtKttst
# CONSENSUS/50%    pRsRGRP+Kssss
#
# using the following ambiguity codes:
# alcohol	      o   S,T
# aliphatic	    l   I,L,V
# any	          .   A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y
# aromatic	    a   F,H,W,Y
# charged	      c   D,E,H,K,R
# negative	    -   D,E
# polar	        p   C,D,E,H,K,N,Q,R,S,T
# positive	    +   H,K,R
# small	        s   A,C,D,G,N,P,S,T,V
# tiny	        u   A,G,S
# turnlike	    t   A,C,D,E,G,H,K,N,Q,R,S,T
# hydrophobic	  h   A,C,F,G,H,I,K,L,M,R,T,V,W,Y

tP <- "[ACDEGHKNQRST]"
pP <- "[HKR]"
(AT80 <- sprintf("%s%s%sRGRP.%s..%s.", tP, pP, tP, pP, tP))   # t+tRGRP.+..t.

(m <- regexec(AT80, mbp1[[1]][1]))  # not found

# === SMART family alignment

SM00384ali <- character()
SM00384ali["PRH_PETCR"]      <- "KRSRGRPRKVQNS"
SM00384ali["O43167"]         <- "KRKRGRPKKVNTL"
SM00384ali["AAC68776"]       <- "KKPRGRPKKYEDS"
SM00384ali["AAC2_DICDI"]     <- "KRKRGRPPKMDEE"
SM00384ali["O49694"]         <- "KKKRGRPRKYAAD"
SM00384ali["O82964"]         <- "KRGRGRPRKAHAM"
SM00384ali["P92954"]         <- "KRGRGRPPKQKTQ"
SM00384ali["P91155"]         <- "PRARGRPRKATNL"
SM00384ali["HMGC_HUMAN_1"]   <- "KRPRGRPRKWPQQ"
SM00384ali["P76546"]         <- "PKKRGRPAKYNEK"
SM00384ali["YN06_CAEEL"]     <- "KRRRGRPRKDEEA"
SM00384ali["O49276"]         <- "PRKRGRPRKSMED"
SM00384ali["Q50887"]         <- "PKKRGRPPKAKTE"
SM00384ali["O45912"]         <- "KKGRGRPAKNPSA"
SM00384ali["O43095"]         <- "TRKRGRPRKTIAS"
SM00384ali["Q22173"]         <- "KKTRGRPRKDRSQ"
SM00384ali["ORC1_SCHPO"]     <- "SRGRGRPRKYPLP"
SM00384ali["Q43877_1"]       <- "PRPRGRPPKDPNA"
SM00384ali["O59212"]         <- "EKGRGRPKKYSTR"
SM00384ali["Q06339"]         <- "AKKRGRPRKSVVA"
SM00384ali["O15026_1"]       <- "KRRRGRPPKARDL"
SM00384ali["JH0797"]         <- "ARKRGRPPKKIQL"
SM00384ali["O65795"]         <- "AKPRGRPAKAAKT"
SM00384ali["Q38778_1"]       <- "GRPRGRPAKAKDP"
SM00384ali["Q40725_3"]       <- "GRGRGRPPKKASS"
SM00384ali["Q01086"]         <- "KRGRGRPRRSVTD"
SM00384ali["Q22381"]         <- "QKRRGRPRKTDAA"
SM00384ali["CPD1_DROME"]     <- "IKKRGRPAKNKGS"
SM00384ali["AAD21618_3"]     <- "KRGRGRPPLHRSE"
SM00384ali["HMGC_HUMAN_2"]   <- "KRPRGRPKGSKNK"
SM00384ali["Q43877_3"]       <- "GRPRGRPKKIART"
SM00384ali["O15030_1"]       <- "ARARGRPRKTKPG"
SM00384ali["AAD21618_6"]     <- "KRKRGRPPLNKPK"
SM00384ali["O23142_1"]       <- "EKKRGRPPGSSSK"
SM00384ali["CPD1_DROME_1"]   <- "VKKRGRPSKASVG"
SM00384ali["Q40725_4"]       <- "KRGVGRPRKNATP"
SM00384ali["CPD1_DROME_2"]   <- "KRKAGRPKKHQPS"
SM00384ali["O49276_1"]       <- "GRARGRPPGVKNG"
SM00384ali["O49276_2"]       <- "RKKRGRPKKFDRI"
SM00384ali["CAB40849_3"]     <- "SRGAGRPPKAKSP"
SM00384ali["O49694_1"]       <- "KRNRGRPPGSGGT"
SM00384ali["CPD1_DROME_5"]   <- "PKKRGRPSLAAGK"
SM00384ali["O00536_1"]       <- "TGKRGRPRNTEKA"
SM00384ali["CPD1_DROME_6"]   <- "TKGRGRPKSSGGA"
SM00384ali["CPD1_DROME_7"]   <- "GGQRGRPPKASKI"
SM00384ali["Q40725_12"]      <- "KRGAGRPRKKRPL"
SM00384ali["O15030_2"]       <- "PKRRGRPPSKFFK"
SM00384ali["CPD1_DROME_8"]   <- "GRGLGRPKKRAVE"
SM00384ali["CAB42096_2"]     <- "KRRPGRPRKHKPE"
SM00384ali["CPD1_DROME_9"]   <- "TKPRSRPAKNIDD"
SM00384ali["O45912_7"]       <- "VLKRGRSVKQPKD"






# ==== Complexity



# Sliding window computation of sequence complexity. The window extends for
# lW residues down from its index.
lW <- 12  # size of Window
hPos <- as.numeric(rep(NA, length(mbp1AA)))  # Record positional entropy
for (i in 1:(length(mbp1AA) - lW)) {
  pmf <- table(mbp1AA[i:(i + lW - 1)])
  hPos[i] <- H(pmf)
}

range(hPos, na.rm = TRUE)

(iMin <- which(hPos == min(hPos, na.rm = TRUE)))
paste0(mbp1AA[iMin:(iMin + lW - 1)], collapse = "")

(iMax <- which(hPos == max(hPos, na.rm = TRUE)))
paste0(mbp1AA[iMax[1]:(iMax[1] + lW - 1)], collapse = "")

which(hPos <= 2.2)

paste0(mbp1AA[108:(111 + lW - 1)], collapse = "")
paste0(mbp1AA[235:(236 + lW - 1)], collapse = "")
paste0(mbp1AA[278:(279 + lW - 1)], collapse = "")
paste0(mbp1AA[295:(296 + lW - 1)], collapse = "")
paste0(mbp1AA[700:(703 + lW - 1)], collapse = "")

subseq <- function(vAA, iFirst, iLast, yText, yTop,
                   lwd = 0.75, col = "#CC0000", cex = 1.0) {
  x <- (iFirst + iLast) / 2
  s <- sprintf("%d-%s", iFirst, paste0(vAA[iFirst:iLast], collapse = ""))
  w <- (strwidth(s) * 1.05) / 2

  rect(x - w, yText - 0.25, x + w, yText - 0.15, col = "#FFFFFFCC", border = NA)

  text(x, yText - 0.2, s, adj = c(0.5, NA), cex = cex)

  segments(x - w,  yText - 0.25, x - w,  yText - 0.15, lwd = lwd, col = col)
  segments(x - w,  yText - 0.15, iFirst, yText - 0.05, lwd = lwd, col = col)
  segments(iFirst, yText - 0.05, iFirst, yTop  + 0.05, lwd = lwd, col = col)

  segments(x + w,  yText - 0.25, x + w,  yText - 0.15, lwd = lwd, col = col)
  segments(x + w,  yText - 0.15, iLast,  yText - 0.05, lwd = lwd, col = col)
  segments(iLast,  yText - 0.05, iLast,  yTop  + 0.05, lwd = lwd, col = col)
}

plot(6:(length(hPos) + 5), hPos,
     type = "l",
     ylim = c(1.5, 3.5),
     xlab = "Sequenz Position",
     ylab = "Entropie (bits)")
abline(h = 2.2, col = "#CCCCFF")

subseq(mbp1AA, 295, 307, hPos[295] - 0.25, hPos[295])
subseq(mbp1AA, 278, 290, hPos[278] - 0.25, hPos[278])
subseq(mbp1AA, 235, 246, hPos[235] - 0.1, hPos[235])
subseq(mbp1AA, 108, 122, hPos[108], hPos[108])
subseq(mbp1AA, 700, 714, hPos[700], hPos[700])

# === dotplot


myFilter <- matrix(numeric(25), nrow = 5)
myFilter[1, ] <- c( 1, 0, 0, 0, 0)
myFilter[2, ] <- c( 0, 1, 0, 0, 0)
myFilter[3, ] <- c( 0, 0, 1, 0, 0)
myFilter[4, ] <- c( 0, 0, 0, 1, 0)
myFilter[5, ] <- c( 0, 0, 0, 0, 1)


sanDots(mbp1AA[300:600], mbp1AA[300:600],
        xlab = "Mbp1", ylab = "Mbp1",
        f = myFilter)

# Multiple Sequence alignment
#
# Get MBP1 named orthologues from myDB
cat("\n\n")
sel <- grep("^MBP1", myDB$protein$name)
cat(sprintf("\n>%s\n%s", myDB$protein$name[sel], myDB$protein$sequence[sel]))
cat("\n\n")

# use eg. for TCoffee at EBI

# Locally computed MSA
sel <- grep("MBP1", myDB$protein$name)
MBP1set <- Biostrings::AAStringSet(myDB$protein$sequence[sel])

# To help us make sense of the alignment we need to add the names for
# the sequences. Names for a seqSet object are held in the ranges slot...

MBP1set@ranges@NAMES <- myDB$protein$name[sel]

MBP1set

# Let's run an alignment with "Muscle"
(msaM <-  msa::msaMuscle( MBP1set, order = "aligned"))

# ... or to see the whole thing (cf. ?MsaAAMultipleAlignment ... print method):
msa::print(msaM, show=c("alignment", "complete"), showConsensus=FALSE)


# You see that the alignment object has sequence strings with hyphens as
# indel-characters. The names are printed to the console. And you also see that
# the order has not been preserved, but the most similar sequences are now
# adjacent to each other.

# You probabaly realize that computing an MSA is not that hard. It's not
# entirely trivial to collect meaningful sequences via e.g. PSI-BLAST ... but
# then computing the alignment gives you a result quickly. But what does it
# mean? What information does the MSA contain?

# Let's have a first look at conserved vs. diverged regions of the MSA. msa
# provides the function msaConservationScore() which outputs a vector of scores.
# The scores are the sum of pairscores for the column: for example a perfectly
# conserved column of histidines would have the following score in our MSA
# of eleven sequences:
#   -  one (H, H) pair score is 8 in BLOSUM62;
#   -  there are (n^2 - n) / 2 pairs that can be formed between amino acids
#        in a column from n sequences;
#   -  therefore the column score is 8 * (11^2 - 11) / 2 == 440


data("BLOSUM62")  # fetch the BLOSUM62 package from the Biostrings package

msaMScores <- msa::msaConservationScore(msaM, substitutionMatrix = BLOSUM62)
plot(msaMScores,
     ylim = c(-200, 1200),
     type = "l", col = "#205C5E",
     xlab = "Alignment Position",
     ylab = "msa::msaConservationScore()")

# That plot shows the well-aligned regions (domains ?) of the sequence, but it
# could use some smoothing. Options for smoothing such plots include calculating
# averages in sliding windows ("moving average"), and lowess() smoothing. Here
# is a quick demo of a moving average smoothing, to illustrate the principle.

wRadius <- 15     # we take the mean of all values around a point +- wRadius
len <- length(msaMScores)
v <- msaMScores
for (i in (1 + wRadius):(len - wRadius)) {
  v[i] <- mean(msaMScores[(i - wRadius):(i + wRadius)]) # mean of values in
  # window around i
}
points(v, col = "#FFFFFF", type = "l", lwd = 4.5)
points(v, col = "#3DAEB2", type = "l", lwd = 3)


# You can set a threshold and use rle() to define ranges of values that fall
# above and below the threshold, and thus approximate domain boundaries:
thrsh <- 30
(highScoringRanges <- rle(v > thrsh))
(idx <- c(1, cumsum(highScoringRanges$lengths)))
for (i in seq_along(highScoringRanges$lengths)) {
  if (highScoringRanges$values[i] == TRUE) { # If this range is above threshold,
    rect(idx[i], thrsh, idx[i+1], max(v),    # ... draw a rectangle
         col = "#205C5E33")                  # ... with a transparent color.
    cat(sprintf("Possible domain from %d to %d\n", idx[i], idx[i+1]))
  }
}

# Getting this right requires a bit of fiddling with the window radius and
# threshold (experiment with that a bit), but once we are satisfied, we can use
# the boundaries to print the MSA alignments separately for domains.

rangeText <- function(txt, iFirst, iLast, yText, hBars,
                      lwd = 0.75, col = "#CC0000", cex = 0.7) {
  x <- (iFirst + iLast) / 2
  dH <- hBars/2
  w <- (strwidth(txt) * 1.05) / 2

  rect(x - w, yText - dH, x + w, yText + (2 * dH),
       col = "#FFFFFFCC", border = NA)

  text(x, yText + dH, txt, adj = c(0.5, NA), cex = cex)

  segments(iFirst,  yText - dH, iFirst,  yText + dH, lwd = lwd, col = col)
  segments(iLast,   yText - dH, iLast,   yText + dH, lwd = lwd, col = col)
  segments(iFirst,  yText, iLast, yText, lwd = lwd, col = col)

}


# Add domain boundaries:
(mbp1Aligned <- as.character(msaM)["MBP1_SACCE"])
rangeText("KilA-N",
          regexec("TGSIMK", mbp1Aligned)[[1]][1],
          regexec("FDFTQT", mbp1Aligned)[[1]][1] + 6,
          yText = 1100, hBars = 100, cex = 1.2)

rangeText("Ank4",
          regexec("NGDTAL", mbp1Aligned)[[1]][1],
          regexec("GALT-TI", mbp1Aligned)[[1]][1] + 7,
          yText = 800, hBars = 80, cex = 0.6)

rangeText("Ank3",
          regexec("TVIHHI", mbp1Aligned)[[1]][1],
          regexec("IKDFSP", mbp1Aligned)[[1]][1] + 6,
          yText = 900, hBars = 80, cex = 0.6)

rangeText("Ank2",
          regexec("QGQTPL", mbp1Aligned)[[1]][1],
          regexec("VFDIDS", mbp1Aligned)[[1]][1] + 6,
          yText = 1000, hBars = 80, cex = 0.6)

rangeText("Ank1",
          regexec("ELHTAF", mbp1Aligned)[[1]][1],
          regexec("GTSIRS", mbp1Aligned)[[1]][1] + 6,
          yText = 1100, hBars = 80, cex = 0.6)

rangeText("SbcCD ATPase",
          regexec("EQHDNE", mbp1Aligned)[[1]][1],
          regexec("-EQIIT", mbp1Aligned)[[1]][1] + 6,
          yText = 1100, hBars = 100, cex = 1.2)



regexec("PHSAPY", mbp1Aligned)[[1]][1]
regexec("SNKEGL", mbp1Aligned)[[1]][1]

# [END]
