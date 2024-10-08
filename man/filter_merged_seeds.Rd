% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_merged_seeds.R
\name{filter_merged_seeds}
\alias{filter_merged_seeds}
\title{Filter MergedSeeds by Number of Terms}
\usage{
filter_merged_seeds(merged_seeds, min_terms)
}
\arguments{
\item{merged_seeds}{A dataframe containing the merged seeds with column named `ClusterIndices`.}

\item{min_terms}{An integer specifying the minimum number of terms required in a cluster.}
}
\value{
The filtered data frame with clusters filtered to include only those with at least `min_terms` terms.
}
\description{
This function filters clusters from the MergedSeeds result of RichCluster,
keeping only clusters with greater than or equal to min_terms # of terms.
}
