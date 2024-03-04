



plot.GFDrmst <- function(x,...){
  # are confidence intervals included in the $res output?
  if(!("lwr_conf" %in% colnames(x$res))){
    warning("Confidence intervals are not included.")
    return(NULL)
  }


  x_lim <- c(min(unlist(x$res[,"lwr_conf"])),
             max(unlist(x$res[,"upr_conf"])))
  L <- nrow(x$res)
  y_lim <- c(0.8, L + 0.2)

 # par(mar = c(4, 5, 2, 0.3))
  plot(x = x$res[,"lwr_conf"], y = x$res[,"upr_conf"],
       xlim = x_lim, ylim = y_lim,
       main = paste(paste(100 * (1 - x$alpha), "%", sep = ""),
                    "simultaneous confidence intervals"),
       #xlab = "",#
       xlab = expression(paste(bold(H)[l], bold(mu))),
       ylab = "", type = "n", yaxt = "n")

  axis(2, at = L:1, las = 2, labels = paste("(",sapply(x$res[,"hyp_matrix"], paste, collapse=", ", sep = ""),")", sep =""))
  for (i_c in 1:L) {
    segments(x$res[[i_c,"lwr_conf"]], L - i_c+1, x$res[[i_c,"upr_conf"]], L +1- i_c)
    segments(x$res[[i_c,"lwr_conf"]], L - i_c +1- 0.1, x$res[[i_c,"lwr_conf"]], L +1- i_c + 0.1)
    segments(x$res[[i_c,"upr_conf"]], L - i_c +1- 0.1, x$res[[i_c,"upr_conf"]], L +1- i_c + 0.1)
    segments(x$res[[i_c,"estimator"]], L - i_c +1- 0.1, x$res[[i_c,"estimator"]], L +1- i_c + 0.1)
  }
}

