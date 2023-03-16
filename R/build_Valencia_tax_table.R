#' Builds a pseudo taxonomic table from Valencia taxa names
#'
#' @param taxa_vec a \code{character} vector with the names of the taxa as formatted by Valencia
#'
#' @return a \code{data.frame} with the pseudo taxonomic table for the Valencia
#' taxonomic labels in \code{taxa_vec}
#' @export
#'
#' @importFrom dplyr tibble as_tibble full_join mutate case_when
#' @importFrom stringr str_detect str_remove str_count
#' @importFrom magrittr set_colnames
build_Valencia_tax_table <- function(taxa_vec){

  tibble(
    valencia_taxa_label = taxa_vec
  ) %>%
    full_join(
      matrix(NA_character_, ncol = length(taxonomic_levels)) %>%
        set_colnames(taxonomic_levels) %>%
        as_tibble(),
      by = character()
    ) %>%
    mutate(
      n_underscores = str_count(valencia_taxa_label, "_"),
      Species =
        case_when(
          str_detect(valencia_taxa_label, "^[gfocpkd]_") ~ Species,
          str_detect(valencia_taxa_label, "rhamnosus") ~ "rhamnosus",
          str_detect(valencia_taxa_label,  "^[A-Z][A-Za-z]*_") ~
            str_remove(valencia_taxa_label, "^[A-Z][A-Za-z]*_"),
          !str_detect(valencia_taxa_label, "_") ~
            valencia_taxa_label,
          TRUE ~
            valencia_taxa_label
        ),
      Genus =
        case_when(
          str_detect(valencia_taxa_label,"^g_") ~
            str_remove(valencia_taxa_label, "^g_"),
          str_detect(valencia_taxa_label,  "^[A-Z][A-Za-z]*_[a-z_0-9]*") ~
            str_remove(valencia_taxa_label, "_[a-z_0-9]*"),
          TRUE ~ Genus
        ) %>%
        str_replace(., "\\.","/"),
      Family = ifelse(str_detect(valencia_taxa_label,"^f_"),
                      str_remove(valencia_taxa_label, "f_"),
                      NA_character_),
      Order = ifelse(str_detect(valencia_taxa_label,"^o_"),
                     str_remove(valencia_taxa_label, "o_"),
                     NA_character_),
      Class = ifelse(str_detect(valencia_taxa_label,"^c_"),
                     str_remove(valencia_taxa_label, "c_"),
                     NA_character_),
      Phylum = ifelse(str_detect(valencia_taxa_label,"^p_"),
                      str_remove(valencia_taxa_label, "p_"),
                      NA_character_),
      Kingdom = ifelse(str_detect(valencia_taxa_label,"^k_"),
                       str_remove(valencia_taxa_label, "k_"),
                       NA_character_),
      Domain = ifelse(str_detect(valencia_taxa_label,"^d_"),
                      str_remove(valencia_taxa_label, "d_"),
                      NA_character_)
    )

}
