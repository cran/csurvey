if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "capture.output", "sh", "xm.red",
    "grid_id", "freq", "prob", "x", "muhat", "fit_type", 
    "lwr", "upp", "domain", "block", "..common_vars", "ynm", "amat",
    "nd", "w", "v1", "xvs2", "xm.red", "ne", "zeros_ps", "small_cells", "obs_cells",
    "zeros", "Nhat", "Ds", "amat_0", "Nd",
    "x1nm", "x2nm", "data", "ci", "transpose", "main", "categ", "categnm", "surface", "type",
    "col", "cex.main", "xlab", "ylab", "zlab", "zlim", "box", "axes", "th", "ltheta",
    "ticktype", "nticks", "palette", "NCOL"
  ))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
