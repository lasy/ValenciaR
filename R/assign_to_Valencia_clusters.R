
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
        "Data for ", length(missing_taxa),"/", nrow(v_clusters)," taxa are missing.\n",
        "These missing taxa represent around ",round(mean(missing_prop)*100),"% (",
        str_c(round(range(missing_prop)*100), collapse = "-"),"%) of Valencia cluster composition.\n",
        "Their list is in the `missing_taxa` component of the returned value."
            )
    )
  }

  if (distance == "YC") {
    distances <- 1 - YC::YC(mat1 = input, mat2 = v_clusters)
  }
  else {
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
