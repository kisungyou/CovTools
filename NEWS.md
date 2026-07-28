# CovTools 0.5.7

* Replaced the external generalized-eigenvalue dependency with internal
  symmetric-definite matrix-pair helpers based on a Cholesky reduction.
* Removed an obsolete native generalized-eigenvalue wrapper.
* Corrected positive-definite hard-threshold selection and now verify every
  returned covariance estimate against the documented eigenvalue tolerance.
* Corrected Power-Euclidean distances and means for negative powers, fixed
  Procrustes size-and-shape distances, and covered scalar and one-slice inputs.
* Corrected bandwidth selection and matrix construction in the Banerjee, An,
  and Lee banded precision estimators.
* Corrected graphical-lasso BIC scoring and return the sparse positive-definite
  ADMM estimate.
* Centered covariance diagonality tests, restored the documented Cai test
  result field, and added minimum-sample and degenerate-data checks.
* Sequential cross-validation and BIC selection no longer create socket
  workers; parallel paths now clean up workers without changing the caller's
  global foreach backend.

# CovTools 0.5.4

* Corrected CRAN errors.

# CovTools 0.5.3

* Fix an error in `CovEst.2010OAS`.

# CovTools 0.5.2  

* Function list is moved to README file.
* Following functions are added - `CovEst.2003LW`, `CovEst.2010OAS`, `CovEst.2010RBLW`.
    
# CovTools 0.5.1

* Functions for hypothesis tests are separated for clarity.
* Initialize the following documentation:
  - `NEWS` for keeping record of updates.
  - `README` to briefly introduce the method.
