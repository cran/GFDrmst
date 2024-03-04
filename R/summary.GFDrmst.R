
## Summaryfunktion
summary.GFDrmst <- function (object, digits = 8,...) {
  cat("#", rep("-", 4 + round((64 - nchar(object$method))/4)), object$method,
      rep("-", 4 + round((64 - nchar(object$method))/4)) ,"#","\n", "\n")
    cat("- Significance level:", object$alpha)
    cat("\n", "\n")
    cat("#--- Hypothesis matrices-----------------------------------------------------------#",
        "\n")
    print(object$res[,1])
    cat("\n")
    if(any(sapply(object$res[,2],length) > 1)){ # estimators are only printed if not visible in res
    cat("#--- Estimators--------------------------------------------------------------------#",
        "\n")
    print(object$res[,2])
    cat("\n")}
    cat("#--- Wald-type test statistics ----------------------------------------------------#",
        "\n")
    print(object$test_stat)
    cat("\n")
    cat("#--- Overall p-values -------------------------------------------------------------#",
        "\n")
    print(round(object$p.value, digits))
    cat("\n")
    cat("#--- Overall results --------------------------------------------------------------#",
        "\n")
    temp_res <- object$res
    for (i in seq_len(nrow(object$res))) {
      for (j in c(2:(ncol(object$res) - 1))) {
        if (object$res[i, j] != "") {
          temp_res[[i, j]] <- round(as.numeric(object$res[[i,j]]), digits)
        }
      }
    }
    print(temp_res)
    cat("#----------------------------------------------------------------------------------#",
        "\n")
}
