###########################################################
## Description
###########################################################

organism2method_des <- function(organism) {
  if (organism == "hsapiens") {
    method_des <- "human genome (GRCh38.p13)"
  } else if (organism == "mmusculus") {
    method_des <- "mouse genome (GRCm39)"
  } else if (organism == "rnorvegicus") {
    method_des <-	"rat genome Rnor_6.0"
  }
  return(method_des)
}

spikein_ERCC2method <- function(spikein_ERCC) {
  if (spikein_ERCC==TRUE) {
    spikein_method <- "Reads were mapped to a composite genome made by concatenating the reference genome and 92 ERCC ExFold RNA Spike-In Mixes sequences (ThermoFisher)."
    } else {
    spikein_method <- ""
    }
  return(spikein_method)
}