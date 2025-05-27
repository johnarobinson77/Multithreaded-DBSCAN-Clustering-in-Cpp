/**
 * Copyright (c) 2024 John Robinson.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */


// DBSCANtest.cpp : This file contains a test of the multithreaded DBSCAN code in MTDBSCAN.hpp. 

//#define COMPARISON_TEST // enables a full comparison of clusters between the multithreaded DBSCAN and single threaded DBSCAN
//#define MEASURE_INDIVIDUAL_TIMES // enables printing of the time to build the k-d tree  which is also included in the total time.
#include <random>
#include <chrono>
#include <thread>
#include <iomanip>
#include "MTDBSCAN.hpp"
#include "avlSet.hpp"
#include "ParseArgs.hpp"
#include <stdlib.h>
#ifdef COMPARISON_TEST
#include "dbscan.h"
#endif

// a slight rewrite of the Romdomer class from
// https://stackoverflow.com/questions/13445688/how-to-generate-a-random-number-in-c/53887645#53887645
class RandomIntervalPDB {
  // random seed by default
  std::mt19937 gen_;
  std::uniform_int_distribution<int64_t> dist_;

public:
  RandomIntervalPDB(int64_t min, int64_t max, unsigned int seed = std::random_device{}())
    : gen_{ seed }, dist_{ min, max } {}

  // if you want predictable numbers
  void SetSeed(unsigned int seed) { gen_.seed(seed); }

  int64_t operator()() { return dist_(gen_); }
};


// set these to the type to be used for the test.  
typedef int64_t dkey_t; // dkey_t for the point coordinate type
typedef int64_t dval_t; // dval_t for the value associated with that point


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
  // command line arg modifiable values for the testing the args
  int64_t clusterSpan = 1000;             // +/- range of the random numbers used to generate a cluster
  int64_t numClusters = 1600;             // number of clusters to generate
  int64_t clusterSpacing = 10000000000L;  // +/- range of random numbers used to generate the cluster center points
  int64_t numPointsPer = 4000;            // number of points per cluster
  int64_t numThreads = -1;                // number of threads.  -1 indicates use the processor number of cores
  int64_t rngSeed = 100;                  // random number generator seed

  const int numDimensions = 3;            // dimensions of the point
  dkey_t searchRad = (dkey_t)(clusterSpan * sqrt(numDimensions) / sqrt(3));  // size of the cluster search window.
  
  ParseArgs parseArgs;
  parseArgs.addArg("-c", "Set number of artificial clusters", &numClusters);
  parseArgs.addArg("-p", "Set Number of points per cluster", &numPointsPer);
  parseArgs.addArg("-s", "Set span of the points in a cluster", &clusterSpan);
  parseArgs.addArg("-r", "Set size of cluster search window", &searchRad);
  parseArgs.addArg("-t", "Set number of threads", &numThreads);
  parseArgs.addCheck([](int64_t* const value) {
    if (*value < 0 || *value == 0)
      throw std::range_error("Number of threads must be > 0 or -1 which indicates use hardware thread count");
  });
  parseArgs.addArg("-g", "Set random number generator seed", &rngSeed);
  parseArgs.addHelpPretext("DBSCANtest <args> is a test of the multithreaded DBSCAN code in MTDBSCAN.hpp");
  if (!parseArgs.parse(argc, argv)) exit(1);
  parseArgs.clear();

  int exitCode = 0;
  
  if (numThreads == -1) numThreads = std::thread::hardware_concurrency();

  //std::cout << parseArgs.getValuesString() << std::endl;
  std::cout << "number of clusters = " << numClusters << " number of points per cluster = " << numPointsPer << " cluster span = " << clusterSpan <<
               "\nsearch window = " << searchRad << " threads = " << numThreads << " random number seed = " << rngSeed << std::endl << std::endl;
  std::cout << "Generating artificially clustered data...." << std::endl;

  // Start by generating the center points of the clusters
  auto riLarge = RandomIntervalPDB(-clusterSpacing, clusterSpacing, static_cast<unsigned int>(rngSeed));
  int64_t* clusterCenters = new int64_t[numClusters*numDimensions];
  for (int64_t i = 0; i < numClusters; i++) {
    for (int64_t k = 0; k < numDimensions; k++)
      clusterCenters[i * numDimensions + k] = riLarge();
  }

  // calculate the total number of points and create a coordinate array to hold them
  size_t numPoints = numClusters * numPointsPer;
  auto coordinates = new std::vector<dkey_t>[numPoints];
  auto riSmall = RandomIntervalPDB(-clusterSpan, clusterSpan, riLarge() % 1000);
  int64_t idx = 0;
  for (int64_t j = 0; j < numPointsPer; j++) {
    for (int64_t i = 0; i < numClusters; i++) {
      for (int64_t k = 0; k < numDimensions; k++) {
        dkey_t tmp = riSmall();
        coordinates[idx].push_back(clusterCenters[i * numDimensions + k] + tmp);
      }
      idx++;
    }
  }

  std::cout << "Adding data to MTDBSCAN..." << std::endl;
  auto mtdbscan = new MTDBSCAN< dkey_t, dval_t, numDimensions>;
  // Add each pair points and values to the MTDBSCAN object
  for (size_t i = 0; i < numPoints; ++i) {
    mtdbscan->addPointWithValue(coordinates[i], static_cast<dval_t>(i));
  }
  // get the search range to about cluster distance window
  // dynamic allocator of window is done to make valgrind happy
  auto windowPtr = new std::vector<dkey_t>(numDimensions);
  std::vector<dkey_t>& window = *windowPtr;
  for (int i = 0; i < numDimensions; i++) {
    window[i] = searchRad;
  }

  std::chrono::steady_clock::time_point startTime;
  std::chrono::steady_clock::time_point endTime;

  std::cout << "Running MTDBSCAN clustering algorithm..." << std::endl;
  mtdbscan->setWindow(window);
  mtdbscan->setNumThreads(numThreads);
  startTime = std::chrono::steady_clock::now();;
  mtdbscan->build();
  endTime = std::chrono::steady_clock::now();

  double totalMtdbscanTime = std::chrono::duration<double>(endTime - startTime).count();
#ifdef MEASURE_INDIVIDUAL_TIMES
  std::cout << "Total MTDBSCAN cluster time: " << std::fixed << std::setprecision(4) << totalMtdbscanTime << " seconds.\n" <<
    "KdTree time: " << mtdbscan->getKdTreeTime() << " seconds" << 
    ",  Search time: " << mtdbscan->getSearchTime() << " seconds" << 
    ",  Merge time: " << mtdbscan->getMergeTime() << " seconds" <<
    ",  copy time: " << mtdbscan->getCopyTime() << " seconds." << 
    ",  Cleanup time: " << mtdbscan->getCleanupTime() << " seconds.\n" << std::endl;
#else
  std::cout << "Total MTDBSCAN cluster time: " << std::fixed << std::setprecision(2) << totalMtdbscanTime << " seconds.\n" << std::endl;
#endif

  if (!mtdbscan->checkClusters(numPoints)) {
    delete[] coordinates;
    delete[] clusterCenters;
    delete windowPtr;
    exit(1);
  }
  if (false) {
    delete[] coordinates;
    delete[] clusterCenters;
    delete windowPtr;
    exit(1);
  }
  
  std::cout << "Checking to see that all values ended up in some cluster" << std::endl;
  
  // allocate a vector of bools, one for each point submitted to dbscan, to check that all points are represented in the clusters
  // dynamic allocation of boolSet done to make valgrind happy 
  auto boolSetPtr = new std::vector<bool>(numPoints);
  std::vector<bool>& boolSet = *boolSetPtr;
  // Using std::fill to set all values to false
  std::fill(boolSet.begin(), boolSet.end(), false);

  // set boolSet[val] true for each val or item in the clusters.  Make sure it is only set once.
  bool passed = true;
  for (size_t i = 0; i < mtdbscan->getNumClusters(); i++) {
    for (dval_t val : *mtdbscan->getClusterValueList(i)) {
      passed &= boolSet[val] == false;
      boolSet[val] = true;
      passed &= val < static_cast<dval_t>(numPoints) && val >= static_cast<dval_t>(0);
    }
  }
  // make sure every bool in boolSet is set, indicating all items are represented.
  if (passed && std::count(boolSet.begin(), boolSet.end(), true) == numPoints) {
    std::cout << "Passed" << std::endl;
  }
  else {
    std::cout << "ERROR: Failed to find the same set values in the clusters as generated." << std::endl;
    delete[] coordinates;
    delete[] clusterCenters;
    delete windowPtr;
    exit(1);
  }
  delete boolSetPtr;

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
  startTime = std::chrono::steady_clock::now();
  dbscan->build();
  endTime = std::chrono::steady_clock::now();

  double totalDbscanTime = std::chrono::duration<double>(endTime - startTime).count();
#ifdef MEASURE_INDIVIDUAL_TIMES
  std::cout << "Total DBSCAN cluster time: " << std::fixed << std::setprecision(4) << totalDbscanTime << " seconds" <<
    ",  KdTree time: " << dbscan->getKdTreeTime() << " seconds,  search time: " << totalDbscanTime - dbscan->getKdTreeTime() << " seconds." << std::endl << std::endl;
#else
  std::cout << "Total DBSCAN cluster time = " << std::fixed << std::setprecision(2) << totalDbscanTime << " seconds" << std::endl << std::endl;
#endif

  if (totalMtdbscanTime > totalDbscanTime)
    std::cout << "WARNING: MTDBSCAN time is greater than DBSCAN time." << std::endl;

  if (!dbscan->checkClusters(numPoints)) {
    exit(1);
  }
  if (false) {
    delete[] coordinates;
    delete[] clusterCenters;
    window.clear();
    exit(1);
  }

  // first check to see if the number of resulting clusters is the same
  if (mtdbscan->getNumClusters() != dbscan->getNumClusters()) {
    std::cout << "ERROR: The number of clusters in MTDBSCAN and DBSCAN do not match." << std::endl;
    exit(2);
  }

  std::cout << "Sorting the contents of all the clusters in MTDBSCAN" << std::endl;
  parallelFor((size_t)0, mtdbscan->getNumClusters(), [&mtdbscan](size_t i) {
    mtdbscan->getClusterValueList(i)->sort();
    });

  std::cout << "Sorting the clusters in MTDBSCAN by first value" << std::endl;
  mtdbscan->sortClustersByFirstValue();

  std::cout << "Sorting the contents of all the clusters in DBSCAN" << std::endl;
  parallelFor((size_t)0, dbscan->getNumClusters(), [&dbscan](size_t i) {
    dbscan->getClusterValueList(i)->sort();
    });
  std::cout << "Sorting the clusters in DBSCAN by first value" << std::endl;
  dbscan->sortClustersByFirstValue();


  // comparing exact contents
  // First sort all the contents of each cluster in both cluster results
  // Then, sort the clusters by the first element in each cluster
  // Now both cluster results should have the values in the sam order so do a linear compare.
  
  parallelFor((size_t)0, mtdbscan->getNumClusters(), [&mtdbscan, &dbscan](size_t i) {
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
    auto dbscanIt = dbscan->getClusterValueList(j)->begin();
    for (; pdbscanIt != mtdbscan->getClusterValueList(i)->end(); pdbscanIt++, dbscanIt++) {
      if (*pdbscanIt != *dbscanIt) {
        std::cout << "ERROR: Clusters with matching first elements do not have matching contents. MTDBSCAN["
          << mtdbscan->getClusterID(i) << "] DBSCAN[" << j << "]" << std::endl;
        exit(2);
      }
    }

    // now make sure the bounds of each cluster match
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
    });

  std::cout << "The contents of MTDBSCAN and DBSCAN are identical." << std::endl << std::endl;

  delete dbscan;

#endif
  delete mtdbscan;
  bool doSortTest = true;
  if (doSortTest) {
    std::cout << "Starting sorting test of a small 5 cluster case." << std::endl;
    auto dbscan1 = new MTDBSCAN<dkey_t, dval_t, 1>();

    const std::vector<dkey_t> smallClusters{ 501, 101, 201, 301, 401, 502, 102, 202, 302, 503, 103, 203, 504, 104, 505 };
    std::vector<dkey_t> tc(1);
    for (size_t i = 0; i < smallClusters.size(); i++) {
      tc[0] = smallClusters[i];
      dbscan1->addPointWithValue(tc, i);
    }
    std::vector<dkey_t> window1{ 2 };
    dbscan1->setWindow(window1);
    dbscan1->buildKdTree();
    dbscan1->buildClusters();
    dbscan1->sortClustersBySize();
    std::cout << "Sorted cluster sized are:";
    for (size_t i = 0; i < dbscan1->getNumClusters(); i++) {
      std::cout << " " << dbscan1->getClusterSize(i);
    }
    std::cout << std::endl;
    // clear and rebuild the cluster
    dbscan1->clear();
    dbscan1->buildClusters();
    dbscan1->sortClustersBySize();
    auto center = dbscan1->getClusterCenter(0);
    std::cout << "Center of the largest cluster is " << center->at(0) << std::endl;
    auto maxBounds = dbscan1->getClusterMaxCorner(0);
    std::cout << "Max corner of the largest cluster is " << maxBounds->at(0) << std::endl;
    auto minBounds = dbscan1->getClusterMinCorner(0);
    std::cout << "Min corner of the largest cluster is " << minBounds->at(0) << std::endl;
    delete dbscan1;
  }

  delete[] coordinates;
  delete[] clusterCenters;
  delete windowPtr;

  exit(0);

}
