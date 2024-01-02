/**
 * Copyright (c) 2023 John Robinson.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */


// DBSCANtest.cpp : This file contains the a test of the multithreaded DBSCAN code in MTDBSCAN.hpp. 

//#define DEBUG_PRINT  // enable this to add debug print messages.  must be defined befor including MTDBSCAN.hpp
//#define COMPARISON_TEST // enables a full comparison of clusters between the multithreaded DBSCAN and single threaded DBSCAN
//#define MEASURE_KDTREE_TIME

#include "MTDBSCAN.hpp"
#ifdef COMPARISON_TEST
#include "dbscan.h"
#endif

// set these to the type to be used for the test.  
typedef int64_t dkey_t; // dkey_t for the point coordinate type
typedef int64_t dval_t; // dval_t for the value associated with that point

// get the current time as a double.  Used for evaluating performance.
double timeNow() {
  struct timespec now;
  if (timespec_get(&now, TIME_UTC) != TIME_UTC) printf("timespec_get failure\n");
  double dblTime = (double)now.tv_sec + (double)now.tv_nsec * 1e-9;
  return dblTime;
}


// The test case generated in main() creates artificially clustered data around 1600 random points 
// in n-space.  Each artificial cluster has 4000 points.  When DBSCAN is run on this data the
// checkClusters function should print out the following result:
// Cluster Count = 1600 Max = 4000  Min = 4000 Average = 4000
// This is followed by a quick check to see that all values were returned
// This test has been run with in64_t, float, and double for the dkey_t type and dimensions 1 to 5.
// Note: the randomInterval class used in this fuction is defined in KdTree.h

//#define EXTREME_CASE1

int main(int argc, const char* argv[])
{
  const int numDimensions = 3;            // dimensions of the point
  int64_t clusterSpan = 1000;             // +/- range of the random numbers used to generate a cluster
  int64_t numClusters = 1600;             // number of clusters to generate
  int64_t clusterSpacing = 10000000000L;  // +/- range of random numbers used to generate the cluster center points
  int64_t numPointsPer = 4000;            // number of points per cluster
  dkey_t searchRad = (dkey_t)(clusterSpan * sqrt(numDimensions) / sqrt(3));  // size of the cluster search window.
  int numThreads = -1;                    // number of threads.  -1 indicates use the processor number of cores
  int rngSeed = 100;                      // random number generator seed

  // override the defaults with command lint switches.
  for (int k = 1; k < argc; k++) {
    if (strcmp(argv[k], "-t") == 0)
      numThreads = strtol(argv[++k], nullptr, 10);
    else if (strcmp(argv[k], "-r") == 0)
      searchRad = (dkey_t)strtoll(argv[++k], nullptr, 10);
    else if (strcmp(argv[k], "-c") == 0)
      numClusters = strtoll(argv[++k], nullptr, 10);
    else if (strcmp(argv[k], "-p") == 0)
      numPointsPer = strtoll(argv[++k], nullptr, 10);
    else if (strcmp(argv[k], "-s") == 0)
      clusterSpan = strtoll(argv[++k], nullptr, 10);
    else if (strcmp(argv[k], "-g") == 0)
      rngSeed = strtol(argv[++k], nullptr, 10);
    else if (strcmp(argv[k], "-h") == 0) {
      std::cout <<
        "-t [int]   Set number of threads\n" <<
        "-r [int]   Set size of cluster search radius\n" <<
        "-c [int]   Set number of clusters\n" <<
        "-p [int]   Number of points per cluster\n" <<
        "-s [int]   Span of the points in a cluster\n" <<
        "-g [int]   Data random number generator seed\n" <<
        "-h         Print this message" <<
        std::endl;
      exit(0);
    } else {
      std::cout << argv[k] << " is not valid.  Use -h for help" << std::endl;
      exit(1);
    }
  }

 // std::cout << "threads = " << numThreads << " searchRad = " << searchRad << 
 //              " numClusters = " << numClusters << " numPointsPer = " << numPointsPer << 
 //              " clusterSpan = " << clusterSpan << " rngSeed = " << rngSeed << std::endl;
  std::cout << "Creating artificially clustered data...." << std::endl;

  // Start by generating the center points of the clusters
  auto riLarge = RandomIntervalPDB(-clusterSpacing, clusterSpacing, rngSeed);
  int64_t* clusterCenters = new int64_t[numClusters*numDimensions];
  for (int64_t i = 0; i < numClusters; i++) {
    for (int64_t k = 0; k < numDimensions; k++)
      clusterCenters[i * numDimensions + k] = riLarge();
  }

  // calculate the total number of points and create a coordinate array to hold them
  size_t numPoints = numClusters * numPointsPer;
  auto coordinates = new std::vector<dkey_t>[numPoints];
  // create a list that will be used later  check that all the right values have been returned.
  auto riSmall = RandomIntervalPDB(-clusterSpan, clusterSpan, riLarge() % 1000);
  int64_t idx = 0;
  bool* indices = new bool[numPoints];
  for (int64_t j = 0; j < numPointsPer; j++) {
    for (int64_t i = 0; i < numClusters; i++) {
      for (int64_t k = 0; k < numDimensions; k++) {
        dkey_t tmp = riSmall();
        coordinates[idx].push_back(clusterCenters[i * numDimensions + k] + tmp);
      }
      indices[idx] = false;
      idx++;
    }
  }

  std::cout << "Adding data to MTDBSCAN..." << std::endl;
  auto mtdbscan = new MTDBSCAN< dkey_t, dval_t, numDimensions>;
  // Add each pair points and values to the MTDBSCAN object
  for (size_t i = 0; i < numPoints; ++i) {
    mtdbscan->addPointWithValue(coordinates[i], i);
  }
  // get the search range to about cluster distance window
  std::vector<dkey_t> window(numDimensions);
  for (int i = 0; i < numDimensions; i++) {
    window[i] = searchRad;
  }

  std::cout << "Running NTDBSCAN clustering algorithm..." << std::endl;
  mtdbscan->setWindow(window);
  mtdbscan->setNumThreads(numThreads);
  double startTime = timeNow();
  mtdbscan->build();
  double endTime = timeNow();

  double totalTime = (endTime - startTime);
#ifdef MEASURE_KDTREE_TIME
  std::cout << "Total MTDBSCAN cluster time: " << std::fixed << std::setprecision(2) << totalTime << " seconds" << 
    ",  KdTree time: " << mtdbscan->getKdTreeTime() << " seconds." << std::endl << std::endl;
#else
  std::cout << "Total MTDBSCAN cluster time: " << std::fixed << std::setprecision(2) << totalTime << " seconds" << std::endl << std::endl;
#endif
  if (!mtdbscan->checkClusters(numPoints)) {
    exit(1);
  }

  std::cout << "Checking to see that all values ended up in some cluster" << std::endl;
  for (size_t i = 0; i < mtdbscan->getNumClusters(); i++) {
    for (dval_t val : *mtdbscan->getClusterValueList(i)) {
      indices[val] = true;
    }
  }

  bool passed = true;
  for (size_t i = 0; i < numPoints; i++) passed &= indices[i];
  if (passed) {
    std::cout << "Passed" << std::endl;
  } else {
    std::cout << "Failed to find all the values in the clusters." << std::endl;
  }

#ifdef COMPARISON_TEST

  // the comparison test runs the same data through the single threaded dbscan code and 
  // checks the the contents of every cluster in each cluster set matches.

  std::cout << "Adding data to DBSCAN..." << std::endl;
  auto dbscan = new DBSCAN< dkey_t, dval_t, numDimensions>;
  // Add each pair points and values to the DBSCAN object
  for (size_t i = 0; i < numPoints; ++i) {
    dbscan->addPointWithValue(coordinates[i], i);
  }
  // get the search range to about cluster distance window
  for (int i = 0; i < numDimensions; i++) {
    window[i] = searchRad;
  }

  std::cout << "Running DBSCAN clustering algorithm..." << std::endl;
  dbscan->setWindow(window);
  startTime = timeNow();
  dbscan->build();
  endTime = timeNow();

  totalTime = (endTime - startTime);
  std::cout << "Total DBSCAN cluster time = " << std::fixed << std::setprecision(2) << totalTime << " seconds" << std::endl << std::endl;

  if (!dbscan->checkClusters(numPoints)) {
    exit(1);
  }

  // first check to see if the number of resulting clusters is the same
  if (mtdbscan->getNumClusters() != dbscan->getNumClusters()) {
    std::cout << "ERROR: The number of clusters in MTDBSCAN and DBSCAN do not match." << std::endl;
    exit(2);
  }

  std::cout << "Sorting the contents of all the clusters in MTDBSCAN" << std::endl;
  for (size_t i = 0; i < mtdbscan->getNumClusters(); i++) {
    mtdbscan->getClusterValueList(i)->sort();
  }
  std::cout << "Sorting the clusters in MTDBSCAN by first value" << std::endl;
  mtdbscan->sortClustersByFirstValue();

  std::cout << "Sorting the contents of all the clusters in DBSCAN" << std::endl;
  for (size_t i = 0; i < dbscan->getNumClusters(); i++) {
    dbscan->getClusterValueList(i)->sort();
  }
  std::cout << "Sorting the clusters in DBSCAN by first value" << std::endl;
  dbscan->sortClustersByFirstValue();


  // comparing exact contents
  // first sort all the contents of each cluster in both cluster results
  // then, for each cluster in the MTDBSCAN, find a cluster in DBSCAN with the same first element
  // once a match has been found, check that the entire contents of both match.
  
  for (size_t i = 0; i < mtdbscan->getNumClusters(); i++) {
    size_t j = i;
    if (*mtdbscan->getClusterValueList(i)->begin() != *dbscan->getClusterValueList(j)->begin()) {
      std::cout << "ERROR: No matching cluster found for MTDBSCAN cluster " << mtdbscan->getClusterID(i) << std::endl;
      exit(2);
    }
    if (mtdbscan->getClusterValueList(i)->size() != dbscan->getClusterValueList(j)->size()) {
      std::cout << "ERROR: Clusters with matching first elements are not the same size. MTDBSCAN[" << mtdbscan->getClusterID(i) 
        << "] DBSCAN[" << j << "]" << std::endl;
      exit(2);
    }
    auto pdbscanIt = mtdbscan->getClusterValueList(i)->begin();
    auto dbscanIt =  dbscan->getClusterValueList(j)->begin();
    for (; pdbscanIt != mtdbscan->getClusterValueList(i)->end(); pdbscanIt++, dbscanIt++) {
      if (*pdbscanIt != *dbscanIt) {
        std::cout << "ERROR: Clusters with matching first elements do not have matching contents. MTDBSCAN[" 
            << mtdbscan->getClusterID(i) << "] DBSCAN[" << j << "]" << std::endl;
        exit(2);
      }
    }

    auto MTDBSCANmin = mtdbscan->getClusterMinCorner(i);
    auto MTDBSCANmax = mtdbscan->getClusterMaxCorner(i);
    auto dbscanmin = dbscan->getClusterMinCorner(j);
    auto dbscanmax = dbscan->getClusterMaxCorner(j);
    bool eq = true;
    for (int k = 0; k < numDimensions; k++) {
      if (MTDBSCANmin->at(k) != dbscanmin->at(k)) eq = false;
      if (MTDBSCANmax->at(k) != dbscanmax->at(k)) eq = false;
    }

    if (!eq) {
      std::cout << "ERROR: Clusters with matching first elements do not have matching bounds. MTDBSCAN["
        << mtdbscan->getClusterID(i) << "] DBSCAN[" << j << "]" << std::endl;
      exit(2);
    }
  }

  std::cout << "The contents of MTDBSCAN and DBSCAN are identical." << std::endl << std::endl;

  delete dbscan;

#endif

  delete mtdbscan;

  std::cout << "Starting sorting test of a small 5 cluster case." << std::endl;
  auto dbscan1 = new MTDBSCAN<dkey_t, dval_t, 1>();

  std::vector<dkey_t> smallc{ 1, 101, 201, 301, 401, 2, 102, 202, 302, 3, 103, 203, 4, 104, 5 };
  for (size_t i = 0; i < smallc.size(); i++) {
    dbscan1->addPointWithValue(std::vector<dkey_t>{smallc[i]}, i);
  }
  std::vector<dkey_t> window1{ 2 };
  dbscan1->setWindow(window1);
  dbscan1->build();
  dbscan1->sortClustersBySize();
  std::cout << "Sorted cluster sized are:";
  for (size_t i = 0; i < dbscan1->getNumClusters(); i++) {
    std::cout << " " << dbscan1->getClusterSize(i);
  }
  std::cout << std::endl;

  delete[] indices;
  delete dbscan1;
  exit(0);

}
