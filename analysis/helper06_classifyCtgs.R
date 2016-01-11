classifyCtgs <- function(x, lfcCutoff, pCutoff) {
  if(x[2] <= pCutoff & x[4] <= pCutoff) {
    
    if(x[1] > lfcCutoff & x[3] > lfcCutoff) {
      "UpUp"
      
    } else if (x[1] > lfcCutoff & x[3] < -1 * lfcCutoff) {
      "UpDown"
      
    } else if (x[1] < -1 * lfcCutoff & x[3] > lfcCutoff) {
      "DownUp"
      
    } else if (x[1] < -1 * lfcCutoff & x[3] < -1 * lfcCutoff) {
      "DownDown"
      
    } else {
      "SigBoring" # Significantly boring
    }
    
  } else if (x[2] <= pCutoff & x[4] > pCutoff) {
    type <- "NS" # Not stat. sig in one of the dimension
    
  } else if (x[2] > pCutoff & x[4] <= pCutoff) {
    type <- "NS" # Not stat. sig in one of the dimension
    
  } else if (x[2] > pCutoff & x[4] > pCutoff) {
    type <- "NSNS" # Not stat. sig in both dimension
  }
}
