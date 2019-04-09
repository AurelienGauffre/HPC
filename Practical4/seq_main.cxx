#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

int      _debug;
#include "kmeans.h"

#include <string>
#include <list>
#include <vector>

#include <chrono>
#include <iostream>
static void usage(char *argv0, float threshold) {
  char *help =
    "Usage: %s [switches] -i filename -n num_clusters\n"
    "       -i filename    : file containing data to be clustered\n"
    "       -b             : input file is in binary format (default no)\n"
    "       -n num_clusters: number of clusters (K must > 1)\n"
    "       -t threshold   : threshold value (default %.4f)\n"
    "       -o             : output timing results (default no)\n"
    "       -d             : enable debug mode\n";
  fprintf(stderr, help, argv0, threshold);
  exit(-1);
}

int main(int argc, char **argv) {
  int     opt;
  extern char   *optarg;
  extern int     optind;
  int     isBinaryFile, is_output_timing;

  int     numClusters, numCoords, numObjs;
  std::string filename;
  container clusters;
  float   threshold;

  std::chrono::high_resolution_clock::time_point
    start, stop;
  std::chrono::duration<double> timing, io_timing, clustering_timing;

  int     loop_iterations;

  // default values
  _debug           = 0;
  threshold        = 0.001;
  numClusters      = 0;
  isBinaryFile     = 0;
  is_output_timing = 0;

  while ( (opt=getopt(argc,argv,"p:i:n:t:abdo"))!= EOF) {
    switch (opt) {
      case 'i': filename = std::string(optarg);
                break;
      case 'b': isBinaryFile = 1;
                break;
      case 't': threshold=atof(optarg);
                break;
      case 'n': numClusters = atoi(optarg);
                break;
      case 'o': is_output_timing = 1;
                break;
      case 'd': _debug = 1;
                break;
      case '?': usage(argv[0], threshold);
                break;
      default: usage(argv[0], threshold);
               break;
    }
  }

  if (filename.length()==0 || numClusters <= 1)
    usage(argv[0], threshold);

  if (is_output_timing)
    start = std::chrono::high_resolution_clock::now();

  // read data points from file
  container objects;
  readfile(filename, objects, isBinaryFile);
  if (is_output_timing) {
    stop      = std::chrono::high_resolution_clock::now();
    io_timing = stop - start;
    start     = stop;
  }

  /* membership: the cluster id for each data object */
  std::valarray<int> membership(objects.size());

  seq_kmeans(objects,
      numClusters, threshold,
      membership,
      clusters,
      loop_iterations);

  if (is_output_timing) {
    stop              = std::chrono::high_resolution_clock::now();
    clustering_timing = stop - start;
    start             = stop;
  }

  /* output: the coordinates of the cluster centres --------------*/
  writefile(filename, clusters, membership);

  /*---- output performance numbers -------------------------------*/

  if (is_output_timing) {
    stop      = std::chrono::high_resolution_clock::now();
    io_timing = stop - start;
    start     = stop;

    std::cout<<"Performing **** Regular Kmeans (sequential version) ****"<<std::endl;

    std::cout<<"Input file:     " << filename <<std::endl ;
    std::cout<<"numObjs       = " << objects.size()<<std::endl ;
    std::cout<<"numCoords     = " << objects[0].size()<<std::endl;
    std::cout<<"numClusters   = " << numClusters <<std::endl;
    std::cout<<"threshold     = " << threshold <<std::endl;

    std::cout<<"Loop iterations    = " << loop_iterations<<std::endl;

    std::cout<<"I/O time           = " << io_timing.count() << "sec" <<std::endl;
    std::cout<<"Computation timing = " << clustering_timing.count()  << "sec" <<std::endl;
  }


  return(0);
}

