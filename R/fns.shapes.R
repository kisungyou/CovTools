# port : 'estcov' but only mean part --------------------------------------
# For simplicity, let's use 2 Procrustes methods as well as AIRM
#' @keywords internal
#' @noRd
port_estcov <- function(S, method=c("Procrustes.SS","Procrustes.Full","AIRM")){
  M = dim(S)[3]
  if (missing(method)){
    stop("* port from shapes : method is not defined.")
  }
  if (!(method %in% c("Procrustes.SS","Procrustes.Full","AIRM"))){
    stop("* port from shapes : method is not valid.")
  }
  method  = match.arg(method)
  weights = rep(1, times=M)
  mcovmat = switch(method,
                   Procrustes.SS = shapes::estSS(S, weights),
                   Procrustes.Full=shapes::estShape(S, weights),
                   AIRM = shapes::estLogRiem2(S, weights)
                   )
  return(mcovmat)
}

# port : 'distcov' --------------------------------------------------------
#' @keywords internal
#' @noRd
port_distcov <- function(A,B,method=c("Procrustes.Full","AIRM")){
  if (missing(method)){
    stop("* port from shapes : method is not defined.")
  }
  if (!(method %in% c("Procrustes.SS","Procrustes.Full","AIRM"))){
    stop("* port from shapes : method is not valid.")
  }
  method = match.arg(method)
  output = switch(method,
                  Procrustes.Full=shapes::distProcrustesFull(A,B),
                  AIRM = shapes::distRiemPennec(A,B)
                  )
  return(output)
}
