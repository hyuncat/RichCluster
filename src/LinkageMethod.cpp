//
//  LinkageMethod.cpp
//  richCluster
//
//  Created by Sarah on 6/2/25.
//

#include <stdio.h>
#include "LinkageMethod.h"

using Cluster = LinkageMethod::Cluster;

double LinkageMethod::computeLinkage(const Cluster& cluster1, const Cluster& cluster2) {
  if (method=="single")
    return LinkageMethod::single(cluster1, cluster2);
  else if (method=="complete")
    return LinkageMethod::complete(cluster1, cluster2);
  else if (method=="average")
    return LinkageMethod::average(cluster1, cluster2);
  else if (method=="ward")
    return LinkageMethod::ward(cluster1, cluster2);
  else
    return 0;
}

double LinkageMethod::single(const Cluster& cluster1, const Cluster& cluster2) {
  double minDist = 100; // a default stupid value
  for (auto i = cluster1.begin(); i!= cluster1.end(); ++i) {
    for (auto j = cluster2.begin(); j!= cluster2.end(); ++j) {
      if (i==j) 
        continue;
      // de-reference the two term iterators
      int t1 = *i;
      int t2 = *j;
      // get distance between them
      double dist = distFct(t1, t2);
      if (dist < minDist)
        minDist = dist;
    }
  }
  return minDist;
}

double LinkageMethod::complete(const Cluster& cluster1, const Cluster& cluster2) {
  double maxDist = 0; // a default stupid value
  for (auto i = cluster1.begin(); i!= cluster1.end(); ++i) {
    for (auto j = cluster2.begin(); j!= cluster2.end(); ++j) {
      if (i==j)
        continue;
      // de-reference the two term iterators
      int t1 = *i;
      int t2 = *j;
      // get distance between them
      double dist = distFct(t1, t2);
      if (dist > maxDist)
        maxDist = dist;
    } 
  }
  return maxDist;
}

double LinkageMethod::average(const Cluster& cluster1, const Cluster& cluster2) {
  double totalDist = 0;
  int n_terms = 0;
  for (auto i = cluster1.begin(); i!= cluster1.end(); ++i) {
    for (auto j = cluster2.begin(); j!= cluster2.end(); ++j) {
      if (i==j)
        continue;
      // de-reference the two term iterators
      int t1 = *i;
      int t2 = *j;
      // get distance between them
      double dist = distFct(t1, t2);
      totalDist += dist;
      n_terms ++;
    }
  }
  return totalDist/n_terms;
} 


double LinkageMethod::ward(const Cluster& cluster1, const Cluster& cluster2) {
  return LinkageMethod::average(cluster1, cluster2); // non-implemented for now
}
