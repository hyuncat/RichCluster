## richCluster
richCluster is a fast C++ agglomerative hierarchical clustering algorithm packaged into easily callable R functions, designed to help cluster biological 'terms' based on how similar of genes are expressed in their activation. 

Terms are clustered together based on how many genes are shared between them. We support two different types of similarity scores:
- Kappa score
- Jaccard index

As well as different linkage criteria for iteratively merging clusters together.
- Multiple linkage (from DAVID implementation)
- Single
- Complete
- Average
- Ward

## Installation
The package is currently under review in submission to CRAN, but for now users can install the package and try out the clustering and visualization features by installing through GitHub.

```shell
install.packages("devtools")
require(devtools)
install_github("hyuncat/richCluster")
```

## Usage
The basic flow is as follows. A demo can be followed along in the `tests.R` script in the `/scripts` subfolder. All demo data can be found in `/inst/extdata`.

<img src="https://github.com/user-attachments/assets/f5c857bc-8547-445f-be15-b9a1102f76f4" width="500" height="auto">

The user can supply a list of enrichment results, and simply run `cluster()` with certain options to return a `ClusterResult` named list object which is all you need to pass in order to the various visualization functions.

## Clustering
The goal of the richCluster package is to simplify the workflow of clustering biological terms across multiple enrichment results. 

- For instance, if a user has 3 different enrichment data corresponding to the same mouse at different stages of a disease and wants to directly compare those results side by side.

Our `cluster()` function accepts any arbitrary number of enrichment results in a list, and then will cluster different `$Term` rows together based on how many shared genes they have in the `$GeneID` column.

### Cluster input formatting
The absolute minimum a user must do is to provide a list of dataframes corresponding to the enrichment results they want to cluster.

- `enrichment_results` - A list of dataframes, each containing enrichment results. Each dataframe should include at least the columns 'Term', 'GeneID', and 'Padj'.

### Basic options
Some basic options to customize clustering outputs include:
- `df_names` - An optional character vector of names for the enrichment result dataframes. Must match the length of `enrichment_results`. Default is NULL.
- `min_terms` - Minimum number of terms each final cluster must include
- `min_value` - Minimum 'Pvalue' a term must have in order to be counted in final clustering

### Distance metric
We also allow user to specify which distance metric / cutoff score they want to use to cluster terms together.
- `distance_metric` - A string specifying the distance metric to use (e.g., "kappa").
- `distance_cutoff` - A numeric value for the distance cutoff (0 < cutoff <= 1).

Note that we technically are using a 'similarity' metric, so the cutoff is the *minimum* kappa score (for instance) that two terms must share in order to be clustered together. Hence a higher cutoff would lead to stricter clustering / smaller clusters.

### Linkage method
After an initial grouping of terms based on having a distance score above a certain cutoff, we merge our initial "seed groups" (pre-clusters) together based on the provided linkage method.

- `linkage_method` - A string specifying the linkage method to use (e.g., "average").
- `linkage_cutoff` - A numeric value between 0 and 1 for the minimum linkage cutoff needed to merge two groups together.

The default recommended should be `average`, which takes the average distance between all terms in the two clusters and uses that as the metric to merge.

The total list of supported metrics includes: 
- `"single"`
- `"complete"`
- `"average"`
- `"ward"` (recommended)

Again, a higher linkage_cutoff leads to stricter (smaller) clusters.

### Output
The output of the `cluster()` function is a `ClusterResult` which can be directly inputted into the visualizations or exported as a csv file with some additional options.

The name of each cluster is determined as the term in the cluster with the highest gene count.

### DAVID-style Clustering
For users who prefer a clustering method similar to the one used by the DAVID functional annotation tool, we provide the `david_cluster()` function. This function implements a clustering algorithm inspired by DAVID's method, which involves creating initial seeds and iteratively merging them.

Key parameters for `david_cluster()`:
- `similarity_threshold`: The kappa score cutoff for considering two terms similar.
- `initial_group_membership`: The minimum number of similar terms required to form an initial seed group.
- `final_group_membership`: The minimum number of terms required in a final cluster.
- `multiple_linkage_threshold`: The threshold for the Dice similarity coefficient used to merge seeds.

Example usage:
```r
david_result <- richCluster::david_cluster(
  enrichment_results,
  df_names = rr_names,
  similarity_threshold = 0.5,
  initial_group_membership = 3,
  final_group_membership = 3,
  multiple_linkage_threshold = 0.5
)
```

## Visualizations
We support two broad categories of cluster visualization:
1. Cluster-level: Aggregates pvalue across all terms in the cluster, to compare different clusters to each other
2. Term-level: Displays individual term pvalues *within* a cluster

Where comparisons are also done across the different enrichment results which originally went into clustering.

Using these two categories, we have three different types of plots which are currently supported:
1. Heatmaps
2. Bar plots
3. Dot plots

As well as an option to export as a dataframe (CSV).

### Network Graphs
`plot_network_graph` visualizes a single cluster as a network graph, where nodes are terms and edges represent the similarity between them. Node colors can represent p-values from multiple datasets.

<img src="https://i.imgur.com/example.png" width="500" height="auto">

`compare_network_graphs_plotly` provides an interactive, side-by-side comparison of network graphs for a single cluster using different p-value types.

<img src="https://i.imgur.com/example.png" width="500" height="auto">

### Heatmaps
`cluster_hmap` displays the -log10(pvalue) of all the different clusters across the user's supplied enrichment results.

<img src="https://github.com/user-attachments/assets/cd5197ba-a1bf-4f40-a2b8-73eca9f1af1b" width="500" height="auto">

`term_hmap` displays the -log10(pvalue) of the union of all terms in the specified clusters, as well as any additional terms specified.

<img src="https://github.com/user-attachments/assets/9ecd9686-d2dc-4548-b51e-61c8a7029817" width="500" height="auto">

### Bar plots
`cluster_bar` displays the -log10(pvalue) of all the different clusters across the user's supplied enrichment results.

<img src="https://github.com/user-attachments/assets/2d856a58-4a30-4870-b0ca-b3a2ea48a23d" width="600" height="auto">

`term_bar` displays the pvalue of all terms within a single cluster across the different input enrichment results.

<img src="https://github.com/user-attachments/assets/169c0a7b-0e0f-4f91-b9cf-8698bb2a13aa" width="500" height="auto">


### Dot plots
`cluster_dot` shows the enrichment data with the size of the cluster being associated with dot radius and the x-axis being the -log10(pvalue).

<img src="https://github.com/user-attachments/assets/1cbe0bfc-e9c7-4137-bac2-75fd5f1da7f0" width="500" height="auto">

`term_dot` shows the enrichment data with the number of genes in each term being associated with dot radius and the x-axis being the -log10(pvalue).

<img src="https://github.com/user-attachments/assets/dc76c127-f7f9-4f6e-9bbd-56c74d2c2107" width="500" height="auto">

### Export as CSV
To view the data in excel / as a dataframe, users can export the final clustered data as a dataframe with `export_df` and save their results to a csv as follows.

```r
cluster_df <- richCluster::export_df(cluster_result)
write.csv(cluster_df, "~/Downloads/cluster_df.csv", row.names=FALSE)
```

Otherwise, to retain access to the entire suite of visualization functions users can also download the entire cluster_result as a list object as follows:
```r
saveRDS(cluster_result, file = "~/Downloads/cluster_result.rds")
```

Then read it back in as an R list as such:
```r
cluster_result <- readRDS('~/Downloads/cluster_result.rds')
```

to avoid re-clustering large datasets across multiple work sessions.


## Troubleshooting
Note: If you are receiving errors, please try renaming your columns of interest to 'Term' and 'GeneID' across all genesets of interest, and make sure the formatting of the columns is consistent (eg, same Term/GeneID name spellings) across all datasets.

If that doesn't work, feel free to raise a GitHub issue or email @ssh2198@columbia.edu and I'll try to help you with your problem.
