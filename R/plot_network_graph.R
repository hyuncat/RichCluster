#' Plot Network Graph for a Cluster
#'
#' This function visualizes a single cluster as a network graph.
#'
#' @param cluster_result The result from the clustering function.
#' @param cluster_num The cluster number to plot.
#' @param distance_matrix The distance matrix used for clustering.
#' @param valuetype_list A list of value types (e.g., "Pvalue_1", "Padj_1") to use for node coloring.
#'
#' @return A plot object.
#'
#' @import igraph
#' @import fields
#' @import viridis
#' @export
plot_network_graph <- function(cluster_result, cluster_num, distance_matrix, valuetype_list) {

  term_names <- cluster_result$cluster_df$Term[cluster_result$cluster_df$Cluster == cluster_num]

  value_list <- list()
  for (valuetype in valuetype_list) {
    raw_values <- cluster_result$cluster_df[[valuetype]][cluster_result$cluster_df$Cluster == cluster_num]
    log_values <- -10 * log10(raw_values)
    value_list[[length(value_list) + 1]] <- log_values
  }

  num_datasets <- length(value_list)
  all_pvals <- unlist(value_list)
  global_min <- min(all_pvals, na.rm = TRUE)
  global_max <- max(all_pvals, na.rm = TRUE)

  node_colors_list <- list()

  for (i in 1:num_datasets) {
    pval <- value_list[[i]]
    color_for_nodes <- rep("grey", length(pval))

    if (any(!is.na(pval))) {
      norm_pval <- (pval - 0) / (global_max - 0)

      palette_colors <- viridis(100)
      palette_indices <- round(norm_pval * 99) + 1
      palette_indices <- pmin(100, pmax(1, palette_indices))

      color_for_nodes <- palette_colors[palette_indices]
      color_for_nodes[is.na(color_for_nodes)] <- "grey"
    }

    node_colors_list[[i]] <- color_for_nodes
  }

  node_colors <- list()
  for (i in 1:length(term_names)) {
    node_color_pieces <- c()
    for (j in 1:num_datasets) {
      node_color_pieces <- c(node_color_pieces, node_colors_list[[j]][i])
    }
    node_colors[[i]] <- node_color_pieces
  }

  vertex_pie <- list()
  for (i in 1:length(term_names)) {
    proportions <- c()
    for (j in 1:num_datasets) {
      proportions <- c(proportions, value_list[[j]][i])
    }
    proportions[is.na(proportions)] <- 0
    if(sum(proportions) == 0) {
        proportions <- rep(1, length(proportions))
    }
    vertex_pie[[i]] <- proportions
  }

  subset_matrix <- distance_matrix[term_names, term_names]
  subset_matrix[subset_matrix < 0] <- NA
  subset_matrix <- 1 / subset_matrix
  subset_matrix[is.infinite(subset_matrix)] <- NA
  subset_matrix[is.na(subset_matrix)] <- 0

  g <- graph_from_adjacency_matrix(subset_matrix, mode = "undirected", weighted = TRUE)

  pvalues <- cluster_result$merged_df$Pvalue
  term_indices_pvalues <- match(term_names, cluster_result$merged_df$Term)
  term_pvalues <- pvalues[term_indices_pvalues]
  max_index <- which.max(term_pvalues)
  representative_term <- term_names[max_index]

  short_term_names <- ifelse(nchar(term_names) > 10, paste0(substr(term_names, 1, 10), "..."), term_names)

  plot(g, vertex.shape = "pie", vertex.pie = vertex_pie,
       vertex.pie.color = node_colors, vertex.label.dist = -3,
       vertex.label.degree = -90, vertex.label = short_term_names,
       vertex.label.cex = 0.7, vertex.size = 30, main = representative_term)

  par(mar = c(5, 4, 4, 6))
  image.plot(
    legend.only = TRUE,
    col = viridis(100),
    zlim = c(0, max(all_pvals, na.rm = TRUE)),
    legend.lab = "-10log10(pvalue)",
    smallplot = c(0.8, 0.85, 0.2, 0.8),
    legend.line = 3,
    legend.mar = 5
  )
}
