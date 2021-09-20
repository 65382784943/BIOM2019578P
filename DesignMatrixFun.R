



#############################################################################################
# R code for paper: Estimating the Optimal Timing of Surgery From Observational Data        #
# Date:             2020-01-18                                                              #
# R version:        3.6.0 (2019-04-26) -- "Planting of a Tree"                              #
# OS:               macOS X Catalina                                                        #
# Platform:         x86_64-apple-darwin15.6.0 (64-bit)                                      #  
#############################################################################################


# Generate the linear spline design matrix A 

lsp <- function(knot = NULL,  data) {
  if (length(knot) == 0) {
    return(data)
  }
  else{
    temp <-
      matrix(data = NA,
             nrow = length(data),
             ncol = length(knot) + 1)
    temp[, 1] <- data
    for (i in 2:(length(knot) + 1)) {
      for (j in 1:length(data)) {
        if ((data[j] - knot[i - 1]) >= 0) {
          temp[j, i] <- (data[j] - knot[i - 1])
        }
        else
          temp[j, i] <- 0
      }
    }
    return(temp)
  }
}












