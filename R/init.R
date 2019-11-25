## messages with start and end

.onAttach <- function(...){
  ## Retrieve Year Information
  date <- date()
  x <- regexpr("[0-9]{4}", date)
  this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)

  # Retrieve Current Version
  this.version = packageVersion("CovTools")

  ## Print on Screen
  packageStartupMessage("**--------------------------------------------------------**")
  packageStartupMessage("** CovTools")
  packageStartupMessage("**  - Geometric & Statistical Tools for Covariance Matrices.")
  packageStartupMessage("** Version    : ",this.version,"      (",this.year,")",sep="")
  packageStartupMessage("** Author     : Kyoungjae Lee, Lizhen Lin, and Kisung You")
  packageStartupMessage("** Maintainer : Kisung You (kyoustat@gmail.com)")
  packageStartupMessage("**")
  packageStartupMessage("** Please share any bugs or suggestions to the maintainer.")
  packageStartupMessage("**--------------------------------------------------------**")
}

.onUnload <- function(libpath) {
  library.dynam.unload("CovTools", libpath)
}
