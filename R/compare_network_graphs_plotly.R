#' Compare Network Graphs using Plotly
#'
#' This function creates a side-by-side comparison of network graphs for a single
#' cluster using different p-value types.
#'
#' @param cluster_result The result from the clustering function.
#' @param cluster_num The cluster number to plot.
#' @param pval_names A list of p-value names to compare.
#'
#' @return A plotly object.
#'
#' @import igraph
#' @import plotly
#' @import viridis
#' @export
compare_network_graphs_plotly <- function(cluster_result, cluster_num, pval_names) {

  # Get all pvals to find global scale
  all_pvals <- list()
  for (pval_name in pval_names) {
    all_pvals[[length(all_pvals) + 1]] <- -10 * log10(cluster_result$cluster_df[[pval_name]][cluster_result$cluster_df$Cluster == cluster_num])
  }

  global_min <- min(unlist(all_pvals), na.rm = TRUE)
  global_max <- max(unlist(all_pvals), na.rm = TRUE)

  # Generate layout once
  term_names <- cluster_result$cluster_df$Term[cluster_result$cluster_df$Cluster == cluster_num]
  subset_matrix <- cluster_result$distance_matrix[term_names, term_names]
  subset_matrix[subset_matrix < 0] <- 0
  g <- graph_from_adjacency_matrix(subset_matrix, mode = "undirected", weighted = TRUE)
  layout <- layout_with_fr(g)

  plots <- list()
  for (i in seq_along(pval_names)) {
    pval_name <- pval_names[[i]]
    show_colorbar <- (i == 1) # Show colorbar only for the first plot

    raw_pval <- cluster_result$cluster_df[[pval_name]][cluster_result$cluster_df$Cluster == cluster_num]
    pval <- -10 * log10(raw_pval)

    layout_df <- as.data.frame(layout)
    colnames(layout_df) <- c("x", "y")

    short_labels <- sapply(term_names, function(x) ifelse(nchar(x) > 10, paste0(substr(x, 1, 10), "..."), x))

    color_values <- rep(NA, length(pval))
    color_values[!is.na(pval)] <- pval[!is.na(pval)]

    hover_text <- paste0(
      "Term: ", term_names, "<br>",
      "P-value: ", signif(pval, 3)
    )

    edges <- as_edgelist(g)
    node_names <- V(g)$name
    edge_indices <- matrix(apply(edges, 2, function(col) match(col, node_names)), ncol = 2)

    edge_df <- data.frame(
      x = layout[edge_indices[,1], 1],
      y = layout[edge_indices[,1], 2],
      xend = layout[edge_indices[,2], 1],
      yend = layout[edge_indices[,2], 2]
    )

    p <- plot_ly() %>%
      add_segments(
        data = edge_df,
        x = ~x, xend = ~xend,
        y = ~y, yend = ~yend,
        line = list(color = 'lightgrey'),
        showlegend = FALSE
      ) %>%
      add_trace(
        data = layout_df,
        x = ~x, y = ~y,
        type = "scatter", mode = "markers+text",
        text = short_labels,
        textposition = "bottom center",
        hovertext = hover_text,
        hoverinfo = "text",
        marker = list(
          size = 20,
          color = color_values,
          colorscale = "Viridis",
          colorbar = if (show_colorbar) list(title = pval_name) else NULL,
          cmin = 0,
          cmax = global_max,
          showscale = show_colorbar,
          line = list(color = "black", width = 1)
        ),
        showlegend = FALSE
      ) %>%
      layout(
        title = list(text = paste0("Cluster ", cluster_num, " (", pval_name, ")"), y = 0.95),
        xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
        yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
      )
    plots[[length(plots) + 1]] <- p
  }

  subplot(plots, nrows = 1, margin = 0.05)
}
