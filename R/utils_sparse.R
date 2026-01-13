suppressPackageStartupMessages({
  library(methods)
  library(Matrix)
  library(reticulate)
})

is_py_none <- function(obj) {
  if (is.null(obj)) {
    return(TRUE)
  }
  inherits(obj, "python.builtin.object") && reticulate::py_is_null_xptr(obj)
}

ensure_py_csc <- function(py_matrix, sparse_module = NULL) {
  if (is_py_none(py_matrix)) {
    return(NULL)
  }

  try({
    if (reticulate::py_has_attr(py_matrix, "tocsc")) {
      return(py_matrix$tocsc())
    }
  }, silent = TRUE)

  if (!is.null(sparse_module)) {
    try({
      return(sparse_module$csr_matrix(py_matrix)$tocsc())
    }, silent = TRUE)
  }

  NULL
}

#' Build a dgCMatrix from Python scipy sparse matrix
#'
#' @param py_matrix Python object exposing CSC/CSR attributes (data, indices, indptr, shape).
#' @return dgCMatrix or NULL when conversion fails
build_dgCMatrix_from_py <- function(py_matrix, sparse_module = NULL) {
  py_matrix <- ensure_py_csc(py_matrix, sparse_module)
  if (is.null(py_matrix)) {
    return(NULL)
  }

  indptr <- reticulate::py_to_r(py_matrix$indptr)
  indices <- reticulate::py_to_r(py_matrix$indices)
  data <- reticulate::py_to_r(py_matrix$data)
  shape <- reticulate::py_to_r(py_matrix$shape)

  methods::new(
    "dgCMatrix",
    Dim = as.integer(shape),
    p = as.integer(indptr),
    i = as.integer(indices),
    x = as.numeric(data)
  )
}
