
#' Assign samples to Valencia clusters
#'
#' @param input a `matrix` or `data.frame` where each row is a sample for
#' which (sub)CST needs to be assigned
#' @param distance a `character` specifying which metric should be used to
#' measure (dis)similarity between sample composition and Valencia clusters
#'
#' @return a `list` with 3 components
#' - `distances`: the distance matrix between each sample (rows of `input`) and
#' Valencia clusters
#' - `assignment`: a `data.frame` with the CST and subCST assignment for each
#' sample (rows of `input`)
#' - `missing_taxa`: a `character` with the list of taxa in Valencia clusters
#' that are missing in the `input`.
#' @export
#'
#' @importFrom stringr str_c str_remove
#'
assign_to_Valencia_clusters <- function(input, distance = "YC"){

  v_clusters <- get_Valencia_clusters()

  if (any(input < 0)) stop("Data to assign must be (positive) proportions or counts\n")
  if (any(input > 1)) {
    warning("Assuming count data (as opposed to proportions) were provided.\n")
    input <- input / rowSums(input)
  }

  if (length(intersect(colnames(v_clusters), colnames(input))) < ncol(v_clusters)) {
    missing_taxa <- colnames(v_clusters)[!(colnames(v_clusters) %in% colnames(input))]
    missing_prop <- v_clusters[, missing_taxa] %>% rowSums()
    warning(
      str_c(
        "`input` does not have counts/proportions for all taxa represented in Valencia.\n",
        "Data for ", length(missing_taxa),"/", ncol(v_clusters)," taxa are missing.\n",
        "These missing taxa represent around ",round(mean(missing_prop)*100),"% (",
        str_c(round(range(missing_prop)*100), collapse = "-"),"%) of Valencia cluster composition.\n",
        "Their list is in the `missing_taxa` component of the returned value."
            )
    )
  }

  if (distance == "YC") {
    distances <- 1 - YC::YC(mat1 = input, mat2 = v_clusters)
  } else if (distance == "BC") {
    distances <- BC(mat1 = input, mat2 = v_clusters)
  } else {
    stop("Not yet implemented for other distances.\n")
  }

  subCST <- rownames(v_clusters)[apply(distances, 1, which.min)]
  CST <- subCST %>% stringr::str_remove(.,"-[A-C]") %>% stringr::str_remove(.,"[0-4]")

  list(
    distances = distances,
    assignment = data.frame(CST = CST, subCST = subCST) %>% set_rownames(rownames(input)),
    missing_taxa = missing_taxa
    )
}


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
