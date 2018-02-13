#  ------------------------------------------------------------------------
# 01. invisible_datamatrix : upgrade of check_datamatrix
# 02. invisible_logical    : single logical variable
# 03. invisible_PosReal    : single positive real number
# 04. invisible_PosIntMM   : single positive integer number with [min,max]
#  ------------------------------------------------------------------------



# 01. invisible_datamatrix ------------------------------------------------
#' @keywords internal
#' @noRd
invisible_datamatrix <- function(A, fname){
  cond1 = (is.matrix(A)||(inherits(A, "Matrix")))
  cond2 = (!any(is.infinite(A)))
  cond3 = (!any(is.na(A)))

  if (cond1&&cond2&&cond3){
    return(1)
  } else {
    stop(paste("* CovTools::",fname," : input matrix X is invalid.", sep=""))
  }
}

# 02. invisible_logical ---------------------------------------------------
#' @keywords internal
#' @noRd
invisible_logical <- function(x, fname, parname){
  cond1 = (length(as.vector(x))==1)
  cond2 = (is.logical(x))
  if (cond1&&cond2){
    return(1)
  } else {
    stop(paste("* CovTools::",fname," : an input ",parname," should be a logical variable.",sep=""))
  }
}


# 03. invisible_PosReal ---------------------------------------------------
#' @keywords internal
#' @noRd
invisible_PosReal <- function(x, fname, parname){
  cond1 = ((length(as.vector(x))==1)&&(is.numeric(x)))
  cond2 = ((x>0)&&(!is.na(x))&&(!is.infinite(x)))
  if (cond1&&cond2){
    return(1)
  } else {
    stop(paste("* CovTools::",fname," : an input ",parname," should be a positive real number.",sep=""))
  }
}

# 04. invisible_PosIntMM --------------------------------------------------
#' @keywords internal
#' @noRd
invisible_PosIntMM <- function(x, fname, parname, minvalue, maxvalue){
  cond1 = (length(as.vector(x))==1)
  cond2 = ((!is.na(x))&&(!is.infinite(x)))
  cond3 = ((x>=minvalue)&&(x<=maxvalue))

  if (cond1&&cond2&&cond3){
    return(1)
  } else {
    stop(paste("* CovTools::",fname," : an input ",parname," should be a positive integer number in [",minvalue,",",maxvalue,"].",sep=""))
  }
}
