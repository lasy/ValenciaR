#' "Converts" a taxonomic count/proportion table to match the Valencia cluster taxonomy.
#'
#' @param input a table of ASV/OTU/species counts or proportions (samples are rows, features are columns).
#' @param tax_table the taxonomic table associated to the columns of \code{input}.
#'
#' @return a \code{list} with two elements. The first element is the
#' \code{conversion_table}, which is the \code{tax_table} augmented with the
#' Valencia label for each columns of \code{input}.
#' The second element is the \code{converted_input}: the counts/proportions from
#' the original \code{input} aggregated for each Valencia taxonomic label.
#' @export
convert_to_Valencia_taxonomy <- function(input, tax_table){


  # TODO: input must be a matrix or a data.frame

  if (any(is.na(input))) {
    warning("Missing values in the `input` table are replaced by 0.\n")
    input[is.na(input)] <- 0
  }

  if (!all(taxonomic_levels %in% colnames(tax_table))) {
    stop(
      paste0("`tax_table` must have the following taxonomic levels: ",
             paste0(taxonomic_levels, collapse = ", "))
    )
  }

  v_clusters <- get_Valencia_clusters()
  v_tax_table <- build_Valencia_tax_table(colnames(v_clusters))
  conversion_table <- match_taxonomies(tax_table, v_tax_table)
  converted_input <- .convert(input, conversion_table)
  list(conversion_table = conversion_table, converted_input = converted_input)
}






#' Matches taxonomic tables
#'
#' @param tax_table a \code{data.frame} with the original taxonomic assignments
#' @param v_tax_table a \code{data.frame} with the
#'
#' @return a \code{data.frame} containing the "conversion table" between the
#' original taxa and the Valencia taxonomic labels.
#' @export
#' @importFrom dplyr mutate case_when inner_join select bind_rows arrange row_number filter left_join
#' @importFrom stringr str_detect
#'
match_taxonomies <- function(tax_table, v_tax_table){

  # General principle:
  # we match the taxonomic label from our input data (given in tax_table)
  # to those used by valencia.

  # TODO: replace the ASV_key by something at the Genus or Species level

  tax_table <-  tax_table %>% mutate(tax_id = row_number())

  # 1. we do the manual matches
  modified_tax_table <-
    tax_table %>%
    mutate(
      Species =
        case_when(
          str_detect(Genus, "Lachnocurva") ~ "BVAB1",
          (Genus == "Sneathia") & str_detect(Species, "amnii") ~ "amnii",
          TRUE ~ Species
        ),
      Genus =
        case_when(
          str_detect(Genus, "Lachnocurva") ~ NA_character_,
          str_detect(Genus, "Prevotella") & str_detect(Species, "corporis") ~ "Prevotella", # instead of "Prevotella_6",
          str_detect(Genus, "Prevotella") & str_detect(Species, "melaninogenica") ~ "Prevotella", # instead of "Prevotella_7",
          str_detect(Genus, "Prevotella") & str_detect(Species, "denticola") ~ "Prevotella", # instead of "Prevotella_7",
          str_detect(Genus, "Prevotella") & str_detect(Species, "copri") ~ "Prevotella", # instead of "Prevotella_9",
          str_detect(Genus, "Prevotella") & str_detect(Species, "histicola") ~ "Prevotella", # instead of "Prevotella_7",

          str_detect(Genus, "Rhizobium") ~ "Rhizobium", # removing the other genera
          str_detect(Genus, "Burkholderia") ~ "Burkholderia/Paraburkholderia", # removing the other genera
          TRUE ~ Genus
        )
    )

  # 2. we fix the species of some v_taxa for which we know different annotation strategies were used
  modified_tax_table <-
    modified_tax_table %>%
    mutate(
      Species =
        case_when(
          Genus == "Gardnerella" ~ "vaginalis",
          (Genus == "Atopobium") & ((Species == "")|is.na(Species)) ~ "vaginae",
          TRUE ~ Species
        )
    )

  # 3. we merge by Genus and Species,
  # keep the remaining ASVs and v_taxa for matching on higher taxonomic rank
  # and remove all unmatched v_taxa that have species-level info
  # (as we will only match on higher taxonomic rank afterwards)

  conversion_table <-
    inner_join(
      modified_tax_table %>% select(tax_id, Genus, Species),
      v_tax_table %>% select(valencia_taxa_label, Genus, Species),
      by = c("Genus", "Species")
    ) %>%
    select(tax_id, valencia_taxa_label)

  remaining_tax_table <-
    modified_tax_table %>% filter(!(tax_id %in% conversion_table$tax_id))
  remaining_v_tax_table <-
    v_tax_table %>%
    filter(
      !(valencia_taxa_label %in% conversion_table$valencia_taxa_label),
      is.na(Species)
    )

  # 4. Then, for each taxonomic rank from the Genus level,
  # we do a merge on that rank between the remaining tables

  for (tax_level in rev(taxonomic_levels[-1])[-1]) {
    # cat(tax_level,"\n")
    if (nrow(remaining_tax_table) > 0 & (nrow(remaining_v_tax_table) > 0)) {
      this_rank_merge <-
        inner_join(
          remaining_tax_table %>% select(tax_id, all_of(tax_level)),
          remaining_v_tax_table %>% select(valencia_taxa_label, all_of(tax_level)) ,
          by = tax_level
        )

      conversion_table <-
        conversion_table %>%
        bind_rows(., this_rank_merge %>% select(tax_id, valencia_taxa_label))

      remaining_tax_table <-
        remaining_tax_table %>% filter(!(tax_id %in% this_rank_merge$tax_id))
      remaining_v_tax_table <-
        remaining_v_tax_table %>% filter(!(valencia_taxa_label %in% this_rank_merge$valencia_taxa_label))
    }
  }

  # at the end, remains the unmatched Valencia taxonomy, and the unmatched ASVs

  # 5. all unmatched ASVs are grouped into d_Eukariota
  conversion_table <-
    conversion_table %>%
    bind_rows(
      .,
      remaining_tax_table %>%
        select(tax_id) %>%
        mutate(valencia_taxa_label = "d_Eukaryota")
    )

  conversion_table <-
    conversion_table %>%
    mutate(
      valencia_taxa_label =
        valencia_taxa_label %>% factor(., levels = v_tax_table$valencia_taxa_label)
      )  %>%
    arrange(valencia_taxa_label)

  conversion_table <-
    conversion_table %>%
    left_join(., tax_table, by = "tax_id")

  conversion_table
}



#' Aggregate the counts/proportions from the input matrix by taxonomic label
#'
#' @param input a `matrix` with counts or proportions (samples are rows, species are columns)
#' @param conversion_table a `data.frame` with at least two columns:
#'
#' @return a `matrix`
#' @importFrom dplyr mutate select arrange
#' @importFrom tidyr pivot_wider
.convert <- function(input, conversion_table) {

  conversion_table_wide <-
    conversion_table %>%
    select(tax_id, valencia_taxa_label) %>%
    mutate(value = 1) %>%
    pivot_wider(
      id_cols = tax_id,
      names_from = valencia_taxa_label, values_from = value,
      values_fill = 0
      ) %>%
    as.data.frame()
  conversion_matrix <-
    conversion_table_wide %>%
    arrange(tax_id) %>%
    select(-tax_id) %>%
    as.matrix()
  converted_mat <- as.matrix(input) %*% conversion_matrix
  converted_mat
}
