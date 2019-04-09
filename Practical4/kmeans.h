#ifndef _H_KMEANS
#define _H_KMEANS

#include <vector>
#include <list>
#include <valarray>
#include <deque>

typedef std::deque<std::valarray<float> > container;

void readfile(std::string, container & ,int);
void writefile(std::string, container &, std::valarray<int> &);
void seq_kmeans(container & objects,
    int     numClusters, float   threshold,
    std::valarray<int> & membership,
    container & clusters,
    int    &loop_iterations);

double  wtime(void);

extern int _debug;

#endif
