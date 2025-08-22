#' @importFrom magrittr %>%
#' @import dplyr
#' @import tidyr
NULL

#' Cluster Terms using DAVID's method
#'
#' This function performs clustering on enrichment results using an algorithm
#' inspired by DAVID's functional clustering method.
#'
#' @param enrichment_results A list of dataframes, each containing enrichment results.
#'        Each dataframe should include at least the columns 'Term', 'GeneID', and 'Padj'.
#' @param df_names Optional, a character vector of names for the enrichment result dataframes. Must
#'        match the length of `enrichment_results`. Default is `NULL`.
#' @param similarity_threshold A numeric value for the kappa score cutoff (0 < cutoff <= 1).
#' @param initial_group_membership Minimum number of terms to form an initial seed group.
#' @param final_group_membership Minimum number of terms for a final cluster.
#' @param multiple_linkage_threshold A numeric value for the merging threshold.
#'
#' @return A named list containing the clustering results.
#'
#' @export
david_cluster <- function(enrichment_results, df_names = NULL,
                          similarity_threshold = 0.5,
                          initial_group_membership = 3,
                          final_group_membership = 3,
                          multiple_linkage_threshold = 0.5) {

  if (is.null(df_names) || length(enrichment_results) != length(df_names)) {
    df_names <- as.character(seq_along(enrichment_results))
  }

  merged_df <- merge_enrichment_results(enrichment_results)
  term_vec <- merged_df$Term
  geneID_vec <- merged_df$GeneID

  cluster_result <- runDavidClustering(
    term_vec,
    geneID_vec,
    similarity_threshold,
    initial_group_membership,
    final_group_membership,
    multiple_linkage_threshold
  )

  cluster_options <- list(
    similarity_threshold = similarity_threshold,
    initial_group_membership = initial_group_membership,
    final_group_membership = final_group_membership,
    multiple_linkage_threshold = multiple_linkage_threshold
  )

  cluster_result$df_list <- enrichment_results
  cluster_result$merged_df <- merged_df
  cluster_result$cluster_options <- cluster_options
  cluster_result$df_names <- df_names
  cluster_result$final_clusters <- cluster_result$clusters
  cluster_result$cluster_df <- make_full_clusterdf(cluster_result$final_clusters, merged_df)

  return(cluster_result)
}
