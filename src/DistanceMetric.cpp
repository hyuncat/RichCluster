//
//  DistanceMetric.cpp
//  richCluster
//
//  Created by Sarah on 6/2/25.
//  edited by Junguk Hur on 8/22/2025
//

#include <stdio.h>
#include "DistanceMetric.h"
#include <unordered_set>
#include <string>
#include <stdexcept>
#include <algorithm>  // for std::set_intersection
#include <iterator>   // for std::inserter

double DistanceMetric::computeDistance(const std::unordered_set<std::string>& t1_genes,
                                       const std::unordered_set<std::string>& t2_genes,
                                       int totalGeneCount) {
  if (metric=="kappa")
    return getKappa(t1_genes, t2_genes, totalGeneCount);
  else if (metric=="jaccard")
    return getJaccard(t1_genes, t2_genes);
  else
    throw std::invalid_argument("unsupported distance metric: " + metric);
}


// the various distance metric computations
// kappa is the standard
double DistanceMetric::getKappa(const std::unordered_set<std::string>& t1_genes,
                                const std::unordered_set<std::string>& t2_genes,
                                int totalGeneCount) {
  
  // Calculate the intersection of t1_genes and t2_genes
  std::unordered_set<std::string> intersection;
  std::set_intersection(t1_genes.begin(), t1_genes.end(), t2_genes.begin(), t2_genes.end(), std::inserter(intersection, intersection.begin()));
  
  double common = static_cast<double>(intersection.size()); // Number of common genes
  
  if (common == 0) {
    return 0.0; // return 0 if no overlapping genes
  } 
  
  double t1_only = t1_genes.size() - common; // Genes unique to t1_genes
  double t2_only = t2_genes.size() - common; // Genes unique to t2_genes
  
  double unique = totalGeneCount - common - t1_only - t2_only; // Count of all genes not found in either term
  
  double relative_observed_agree = (common + unique) / totalGeneCount;
  double chance_yes = ((common + t1_only) / totalGeneCount) * ((common + t2_only) / totalGeneCount);
  double chance_no = ((unique + t1_only) / totalGeneCount) * ((unique + t2_only) / totalGeneCount);
  double chance_agree = chance_yes + chance_no;
  
  if (chance_agree == 1)
    return 0.0; // prevent divide by zero
  else
    return (relative_observed_agree - chance_agree) / (1 - chance_agree); // return kappa!
}

double DistanceMetric::getJaccard(const std::unordered_set<std::string>& t1_genes,
                                  const std::unordered_set<std::string>& t2_genes) {
  std::unordered_set<std::string> intersection; 
  double common = 0;
  for (const auto& gene : t1_genes) {
    if (t2_genes.count(gene)) {
      ++common;
    }
  }
  // double common = static_cast<double>(intersection.size()); // number of common genes
  double total = static_cast<double>(t1_genes.size()) + static_cast<double>(t2_genes.size());
  
  return common / total;
}
