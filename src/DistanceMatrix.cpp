//
//  DistanceMatrix.cpp
//  richCluster
//
//  Created by Sarah on 6/2/25.
//

#include <stdio.h>
#include "DistanceMatrix.h"
#include <Rcpp.h>

int DistanceMatrix::getDistanceIndex(int t1, int t2) const {
  int row_index = t1;
  int col_index = t2;
  return (row_index * n_terms) + col_index;
}

double DistanceMatrix::getDistance(int t1, int t2) const {
  int i = DistanceMatrix::getDistanceIndex(t1, t2);
  return distances[i];
} 

void DistanceMatrix::setDistance(double distance, int t1, int t2) {
  int i = DistanceMatrix::getDistanceIndex(t1, t2);
  distances[i] = distance;
}

// export utility to R
Rcpp::NumericMatrix DistanceMatrix::export_r() const {
  Rcpp::NumericMatrix dm(n_terms, n_terms);
  // unpack the vector 
  for (int i = 0; i < n_terms; ++i) {
    for (int j = 0; j < n_terms; ++j) {
      dm(i, j) = getDistance(i, j);
    }
  }
  Rcpp::List dimnames = Rcpp::List::create(terms, terms);
  dm.attr("dimnames") = dimnames;
  return dm;
}
