library(dplyr)
library(tidyr)

load_cluster_data <- function(from_scratch=FALSE)
{
  if (from_scratch == FALSE)  {
    # load the cluster data from the files
    cluster_result <- readRDS('inst/extdata/cluster_result.rds')
  }
  else {
    # read and manually cluster
    rr1 <- read.delim(system.file("extdata", "HF36wk_vs_HF12wk.txt", package="richCluster"))
    rr2 <- read.delim(system.file("extdata", "HF36wk_vs_WT12wk.txt", package="richCluster"))
    
    enrichment_results <- list(rr1, rr2)
    rr_names <- c('hf36_vs_hf12', 'wt36_vs_wt12')
    
    cluster_result <- richCluster::cluster(
      enrichment_results, df_names=rr_names, min_terms=3, min_value=0.0001,
      distance_metric="kappa", distance_cutoff=0.5,
      linkage_method="average", linkage_cutoff=0.5
    )
    saveRDS(cluster_result, file = "inst/extdata/cluster_result.rds")
  }
  return(cluster_result)
}

cluster_result <- load_cluster_data(from_scratch=TRUE)

# WARD LINKAGE TEST
# ---
ward_cluster_result <- richCluster::cluster(
  cluster_result$df_list, df_names=cluster_result$df_names, min_terms=3, min_value=0.0001,
  distance_metric="kappa", distance_cutoff=0.5,
  linkage_method="ward", linkage_cutoff=0.5
)
print("Ward linkage clustering successful:")
print(head(ward_cluster_result$final_clusters))

# DAVID CLUSTER TEST
# ---
david_cluster_result <- richCluster::david_cluster(
  cluster_result$df_list, df_names=cluster_result$df_names
)
print("David clustering successful:")
print(head(david_cluster_result$final_clusters))

# ALL VISUALIZATION TESTS
# ---
plot_network_graph(
  cluster_result,
  cluster_num = 1,
  distance_matrix = cluster_result$distance_matrix,
  valuetype_list = c("Pvalue_1", "Padj_1")
)

compare_network_graphs_plotly(
  cluster_result,
  cluster_num = 1,
  pval_names = c("Pvalue_1", "Padj_1")
)

c_hmap <- richCluster::cluster_hmap(cluster_result)
c_hmap

# clusters 4, 6, 8 + terms from clusters x,y,z
clusters <- c("blood vessel development", "response to lipoprotein particle", "positive regulation of cell death")
terms <- c("myelination", "lipid oxidation")
t_hmap <- richCluster::term_hmap(cluster_result, clusters, terms, value_type="Padj")
t_hmap

c_bar <- richCluster::cluster_bar(cluster_result)
c_bar

t_bar <- richCluster::term_bar(cluster_result, 1)
t_bar

c_dot <- richCluster::cluster_dot(cluster_result)
c_dot

t_dot <- richCluster::term_dot(cluster_result, 1)
t_dot

cluster_df <- richCluster::export_df(cluster_result)

# write.csv(cluster_df, "~/Downloads/cluster_df.csv", row.names=FALSE)


# example of what workflow using the functions looks like
test_workflow <- function(cluster_result, min_terms=5) {
  
  # Testing: Choose one to comment out and test
  # result <- richCluster::all_clusters_hmap(full_clusterdf, "Padj")
  # result <- richCluster::cluster_correlation_hmap(final_clusters, distance_matrix, 3)
  result <- richCluster::cluster_network(final_clusters, distance_matrix, 1)
  # result <- richCluster::full_network(distance_matrix[1:30, 1:30])
  return(result)
  
}