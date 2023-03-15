

#' Label topics to match the closest Valencia sub-CST
#'
#' @param lda_models a single LDA model (output of `topicmodels::LDA`) or
#' a list of LDA models (e.g., output of `alto::fit_lda_models`).
#' @param tax_table the taxonomic table for each row of the `beta` matrix
#' @param distance the distance used to compute the (dis)similarity between the
#' topic composition and the Valencia centroids.
#' @param max_distance the maximum distance allowed between a subCST and a topic
#' for the topic to be labeled as "other" ("-O") instead of with the subCST label.
#'
#' @return a single LDA model or a list of LDA models (depending on the input)
#' where the rownames of `beta` and the colnames of `gamma` match the label
#' of the closest Valencia sub-CST
#' @export
#'
label_topics <- function(lda_models, tax_table, distance, max_distance = 0.5){

  input <- lda_models
  if (!is.null(lda_models$beta)) {
    input <- list(unique_k = lda_models)
  }

  res <-
    purrr::map(
      .x = input,
      .f = .label_topics_1_model,
      tax_table = tax_table,
      distance = distance,
      max_distance = max_distance
    )

  if (!is.null(lda_models$beta)) {
    res <- res[[1]]
  }
  res
}


#' Label topics from a single model to match the closest Valencia sub-CST
#'
#' @param lda_model a single LDA model (output of `topicmodels::LDA`).
#' @param tax_table  the taxonomic table for each row of the `beta` matrix
#' @param distance the distance used to compute the (dis)similarity between the
#' topic composition and the Valencia centroids.
#' @param max_distance the maximum distance allowed between a subCST and a topic
#' for the topic to be labeled as "other" ("-O") instead of with the subCST label.
#'
#' @return a single LDA model where the rownames of `beta` and the colnames
#' of `gamma` match the label of the closest Valencia sub-CST
#'
#' @importFrom stringr str_replace str_c
#' @importFrom dplyr tibble mutate case_when group_by ungroup n
.label_topics_1_model <-
  function(lda_model, tax_table, distance, max_distance){

    beta_Valencia_tax <-
      convert_to_Valencia_taxonomy(
        lda_model$beta %>% exp(),
        tax_table = tax_table
      )

    clusters <-
      assign_to_Valencia_clusters(
        beta_Valencia_tax$converted_input,
        distance = distance
      )

    tmp <-
      tibble(
        topic_name =
          clusters$assignment$subCST %>%
          stringr::str_replace("I-[A-B]","I"),
        CST = clusters$assignment$CST,
        distance = clusters$distances %>% apply(., 1, min)
      ) %>%
      mutate(
        topic_name =
          case_when(
            distance > max_distance ~ str_c(CST, "-O"),
            TRUE ~ topic_name
          )
      )

    topic_names <- tmp$topic_name

    if (any(duplicated(topic_names))) {
      tmp <-
        tmp %>%
        group_by(topic_name) %>%
        mutate(
          n = rank(distance),
          N = n(),
          topic_name_suffix =
            case_when(
              N == 1 ~ "",
              N > 1 ~ str_c(".",letters[n])
            ),
          new_topic_name =
            str_c(topic_name, topic_name_suffix)
        )  %>%
        dplyr::ungroup()
      topic_names <-  tmp$new_topic_name
    }

    rownames(lda_model$beta) <- topic_names
    colnames(lda_model$gamma) <- topic_names

    lda_model

  }
