#' Computes the Bray-Curtis dissimilarity between the rows of two matrices
#'
#' @param mat1 the matrix with the first set of samples
#' @param mat2 (optional) the matric with the second set of samples
#'
#' @return a distance matrix whose dimensions are `nrow(mat1)` x `nrow(mat2)`
#' @export
#'
BC <- function(mat1, mat2 = NULL) {
  if (is.null(mat2))
    mat2 <- mat1

  # we harmonize columns
  mat <-
    bind_rows(
      mat1 %>% as.data.frame(),
      mat2 %>% as.data.frame()
    ) %>%
    as.matrix()
  # we remove the NAs
  mat[is.na(mat)] = 0

  # re-convert to matrices
  m1 <-
    mat[1:nrow(mat1), ] %>%
    matrix(nrow = nrow(mat1), ncol = ncol(mat))
  m2 <-
    mat[nrow(mat1) + (1:nrow(mat2)), ] %>%
    matrix(nrow = nrow(mat2), ncol = ncol(mat))

  # we compute the distance
  # BC_ij = 1 - 2 * C_ij / (S_i + S_j)
  # where C_ij = sum of the minimum values between i and j
  dims <- c(nrow(m1), ncol(mat), nrow(m2))
  M1 <- array(m1, dim = dims)
  M2 <- array(rep(t(m2), each = nrow(m1)), dim = dims)
  C <- apply(pmin(M1, M2), c(1, 3), sum)
  S1 <- apply(M1, c(1,3), sum)
  S2 <- apply(M2, c(1,3), sum)
  D <- 1 - 2 * C / (S1 + S2)

  # we relabel the rows and columns
  rownames(D) <- rownames(mat1)
  colnames(D) <- rownames(mat2)
  D

}
