/**
 *   File:         file_io.c
 *   Description
 *      This program reads point data from a file
 *      and write cluster output to files
 *   Input file format
 *      ascii  file: each line contains 1 data object
 *      binary file: first 4-byte integer is the number of data
 *      objects and 2nd integer is the number of features (or
 *      coordinates) of each object
 */

#include "kmeans.h"

// Store read objects in a list
#include <algorithm>

//
#include <fstream>
#include <iostream>
#include <sstream>

/** Read input file
 * @param[in]  filename pointer to the name of the file
 * @param[out] objects  hold the data read from files
 */
void readfile(std::string filename,
              container &objects, int isbinary)
{

  std::list<float> buffer_l;
  std::valarray<float> object;

  std::ifstream infile;
  if (isbinary == 1)
  {
    infile.open(filename, std::ios::binary);
    int nobj, ncoord;
    float val;
    infile.read(reinterpret_cast<char *>(&nobj), sizeof(int));
    infile.read(reinterpret_cast<char *>(&ncoord), sizeof(int));

    object.resize(ncoord);

    while (!infile.eof())
    {
      for (int k = 0; k < ncoord; k++)
      {
        infile.read(reinterpret_cast<char *>(&val), sizeof(float));
        object[k] = val;
      }
      objects.push_back(object);
    }
  }
  else
  {
    infile.open(filename);

    if (!infile)
    {
      std::cerr << "Could not open file";
      exit(1);
    }

    std::string line;
    float buffer_d;
    int count = 0;
    std::list<float>::iterator it;
    while (std::getline(infile, line))
    {
      std::istringstream strstream(line);
      if (count == 0)
      {
        // Ignore id of the element
        strstream >> buffer_d;
        while (strstream)
        {
          strstream >> buffer_d;
          buffer_l.push_back(buffer_d);
          count++;
        }
        object.resize(count);
        it = buffer_l.begin();
        for (int i = 0; i < count; i++, it++)
          object[i] = *it;

      }
      else
      {
        count = 0;
        // strstream >> buffer_d;
        while (strstream)
        {
          strstream >> object[count];
          count++;
        }
      }

      buffer_l.clear();
      objects.push_back(object);
    }
  }
  infile.close();
}

void writefile(std::string filename,
               container &clusters,
               std::valarray<int> &membership)
{

  std::ofstream outfile;

  // Coordinates of the cluster centers
  std::string outname;
  outname = filename + ".clusters_centers";
  std::cout << "Writing coordinates of K= " << clusters.size()
            << " clusters centers to file " << outname << std::endl;

  outfile.open(outname);
  int k = 0;
  for (auto cluster : clusters)
  {
    outfile << k << " ";
    for (auto coord : cluster)
    {
      outfile << coord << " ";
    }
    outfile << std::endl;
    k++;
  }
  outfile.close();

  // Write closest cluster center to each data points
  outname = filename + ".membership";
  std::cout << "Writing membership of N= " << membership.size() << " data objects to file " << outname << std::endl;

  outfile.open(outname);
  k = 0;
  for (auto id : membership)
  {
    outfile << k << " " << id << std::endl;
    k++;
  }
  outfile.close();
}
