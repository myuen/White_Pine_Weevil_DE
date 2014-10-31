classifyCtgs <- function(x) {
  if(x[2] <= pCutoff & x[4] <= pCutoff) {
    
    if(x[1] > lfc & x[3] > lfc) {
      "UpUp"
      
    } else if (x[1] > lfc & x[3] < -1 * lfc) {
      "UpDown"
      
    } else if (x[1] < -1 * lfc & x[3] > lfc) {
      "DownUp"
      
    } else if (x[1] < -1 * lfc & x[3] < -1 * lfc) {
      "DownDown"
      
    } else {
      "SigBoring" # Significantly boring
    }
    
  } else if (x[2] <= pCutoff & x[4] > pCutoff) {
    type <- "NS" # Not stat. sig in 1-dimension
    
  } else if (x[2] > pCutoff & x[4] <= pCutoff) {
    type <- "NS" # Not stat. sig in 1-dimension
    
  } else if (x[2] > pCutoff & x[4] > pCutoff) {
    type <- "NSNS" # Not stat. sig in 2-dimension
  }
}
