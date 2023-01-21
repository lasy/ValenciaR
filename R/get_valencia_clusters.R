#' Retrieve Valencia centroids
#'
#' @return \code{matrix} composition of clusters (rows are subCST, columns are taxa)
#' @export
#'
#' @importFrom magrittr %>%
#' @examples get_Valencia_clusters()
#'
get_Valencia_clusters <- function() {
  clusters <-
    readr::read_csv(
    file = "https://raw.githubusercontent.com/ravel-lab/VALENCIA/master/CST_centroids_012920.csv",
    show_col_types = FALSE
  ) %>%
    as.data.frame()

  rownames(clusters) <- clusters$sub_CST
  clusters <- clusters %>% dplyr::select(-sub_CST) %>% as.matrix()
  clusters
}
