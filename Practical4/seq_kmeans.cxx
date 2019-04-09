/** @file  seq_kmeans.c  (sequential version)
 *  @description
 *         Implementation of simple k-means clustering algorithm
 *         This program takes an array of N data objects, each with
 *         M coordinates and performs a k-means clustering given a
 *         user-provided value of the number of clusters (K)
 *         The clustering results are saved in 2 arrays:
 *         1 - clustersK][N] indicating the center coordinates
 *            of K clusters
 *         2 - membership[N] stores the cluster center ids, each
 *            corresponding to the cluster a data object is assigned
 */

#include "kmeans.h"

#include <iostream>
#include <valarray>

/** Find the nearest cluster in the sense of the Euclid norm
 * @param[in] object      current object
 * @param[in] clusters    set of clusters
 */
__inline static
size_t find_nearest_cluster(std::valarray<float> object,
    container clusters){

  // Index of the closest cluster
  size_t index = 0;

  // Current distance and minimum distance
  float dist, min_dist;

  // Compute the distance to first cluster center
  min_dist = (std::pow(object- clusters[0],2)).sum();
  for (size_t i=1; i<clusters.size(); i++) {
    // Compute distance to current cluster center
    dist = (std::pow(object- clusters[i],2)).sum();

    // Keep track of the min_dist and its index
    if (dist < min_dist) {
      min_dist = dist;
      index    = i;
    }
  }
  return(index);
}

/** Apply k-means clustering algorithm
 * @param[in] objects     the set of objects that we are trying to cluster
 * @param[in] numClusters target number of clusters
 * @param[in] threshold   value under which the result of k-mean is satisfactory
 * @param[out] membership  track to which cluster the objects are closest
 * @param[out] loop        number of iteration performed to stable solution
 */
void seq_kmeans(container & objects,
    int     numClusters,
    float   threshold,
    std::valarray<int> & membership,
    container & clusters,
    int    &loop) {

  // Current variation of the computed cluster
  float delta;

  // Various counters
  int   i, index;

  // Clusters holding previous and current states
  container newClusters;

  // Keep track of the added number of objects in current clusters
  std::valarray<int> newClusterSize;

  // pick first k-clusters elements of objects[]
  // as initial cluster centers
  container::iterator it = objects.begin();
  for (i=0; i<numClusters; i++,it++)
    clusters.push_back(*it);

  // Initialize re-clustering holder to null value
  newClusters = clusters;
  newClusterSize.resize(numClusters);

  // Initialize membership
  membership.resize(objects.size());
  membership= -1;


  // While the method did not converge, we iterate
  do {
    // Reset new clusters
    for(auto obj=newClusters.begin(); obj!=newClusters.end(); obj++){
      for(size_t i = 0; i<(*obj).size();i++){
        (*obj)[i]= 0.;
      }
    }
    newClusterSize = 0.;

    // Keep track of change during one iteration
    delta = 0.0;
    i=0;
    for(auto obj:objects) {
      // find the array index of nestest cluster center
      index = find_nearest_cluster(obj, clusters);

      // if membership changes, increase delta by 1
      if (membership[i] != index) delta += 1.0;

      // assign the membership to object i
      membership[i] = index;

      // update new cluster centers : sum of objects located within
      newClusterSize[index]++;
      newClusters[index]+=obj;
      i++;
    }
    // average the sum and replace old cluster centers with newClusters
    for (i=0; i<numClusters; i++) {
      if (newClusterSize[i] > 0)
        clusters[i] = newClusters[i] / float(newClusterSize[i]);
    }
    // Evaluate average data movemement
    delta /= objects.size();
  } while (delta > threshold && loop++ < 500);
  loop++;
}
