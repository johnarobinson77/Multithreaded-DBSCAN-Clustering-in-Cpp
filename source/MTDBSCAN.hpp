/**
 * Copyright (c) 2023 John Robinson.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

 //
 //  MTDBSCAN.hpp
 //
 //  Created by John Robinson on 11/18/23.
 //


#ifndef MTDBSCAN_HPP
#define MTDBSCAN_HPP

#include <iostream>
#include <algorithm>
#include <atomic>
#include <stdexcept>
#include "KdTreePDB.hpp"
#include "avlSet.h"


template <typename K, typename V, size_t N>
class MTDBSCAN {

  // This is a shortened name of an std::pair that contains a tuple in the fist spot and a value in the second
  // When the KdTree search functions are called, a list of these pairs are returned.
  typedef std::pair<K*, std::list<V>*> retPair_t;

  int64_t numThreads = -1;  //number of threads to be used throughout the process
  int64_t numValues = 0;    //number of point,value pairs added to the cluster. 

  // this struct maintains a bounding hyper rectangle
  struct bounds {
    std::vector<K> max;
    std::vector<K> min;
    void addToBounds(std::vector<K>& nv) {
      for (size_t i = 0; i < max.size(); i++) {
        max[i] = (nv[i] > max[i]) ? nv[i] : max[i];
        min[i] = (nv[i] < min[i]) ? nv[i] : min[i];
      }
    }
    void addToBounds(K nv[]) {
      for (size_t i = 0; i < max.size(); i++) {
        max[i] = (nv[i] > max[i]) ? nv[i] : max[i];
        min[i] = (nv[i] < min[i]) ? nv[i] : min[i];
      }
    }
    void addToBounds(bounds& addBounds) {
      addToBounds(addBounds.max);
      addToBounds(addBounds.min);
    }
    void resetBounds() {
      for (size_t i = 0; i < max.size(); i++) {
        max[i] = std::numeric_limits<K>::min();
        min[i] = std::numeric_limits<K>::max();
      }
    }
    bool equals(bounds& otherBounds) {
      bool eq = true;
      for (size_t i = 0; i < max.size(); i++) {
        if (max[i] != otherBounds.max[i]) eq = false;
        if (min[i] != otherBounds.min[i]) eq = false;
      }
      return eq;
    }
    void resetBounds(int n) {
      max.resize(n);
      min.resize(n);
      resetBounds();
    }
    bounds(size_t n) {
      max.resize(n);
      min.resize(n);
      resetBounds();
    }
    bounds() = delete;
  };
  typedef bounds bounds_t;

  // Each cluster must have a unique ID.  This counter is incremented each time a new cluster is
  // created.  Since clusters are created in multiple threads, the increment process must be mutex locked.
  uint32_t IDCounter = 1;
  std::mutex IDCounterMutex;


  // Object of this structure contain the list of values of a single cluster
  struct Cluster {
    std::list<V>* values = nullptr;  // the value list
    bounds_t clusterBounds{ N }; // the minimum corner of a bounding hyper rectangle
    std::string* tag = nullptr; // optional tag string
    avlSet<Cluster*>* overlaps = nullptr; // pointer to a cluster that this cluster overlaps with.
    uint32_t ID = 0;



    // constructor for cluster.
    Cluster(uint32_t& IDCounter, std::mutex& IDCounterMutex) {
      values = new std::list<V>();
      tag = nullptr;
      { // mutex lock this operation.
        std::lock_guard<std::mutex> guard(IDCounterMutex);
        ID = ++IDCounter;
      }
    }

    inline uint32_t getID() { return ID; }

    bool addToCluster(K* kdKey, std::list<V>* values) {
      this->values->splice(this->values->end(), *values);
      clusterBounds.addToBounds(kdKey);
      return true;
    }

    bool combine(Cluster* ocl) {
      values->splice(this->values->end(), *ocl->values);
      clusterBounds.addToBounds(ocl->clusterBounds);
      return true;
    }

    // destructor
    ~Cluster() {
      if (values != nullptr) delete values;
      if (tag != nullptr) delete tag;
      if (overlaps != nullptr) delete overlaps;
    }

  private:
    Cluster() {}
  };

  // this is the list of clusters.
  std::vector<Cluster*>* clusters;
  // Pointer to the temporary KdTree that will be used to accelerated the clustering.
  KdTreePDB< K, V, Cluster, N>* kdTree;
  // this vector contains the distance to the edge of the search window in each dimension.
  std::vector<K> clusterWindowRadius;

  bool buildCluster(std::list<Cluster*>* clusters, std::list<Cluster*>* overlaps, int64_t threadNum, bool oneTime = false) {
    std::vector<K> qLower(N);  //allocate 2 vectors for the upper and lower corners of the window
    std::vector<K> qUpper(N);
    auto keys = new std::list<K*>();  // holds the list of the points to be searched for a cluster.
    retPair_t retPair;
    std::list<retPair_t*> retPairs;  // list of returned pairs from the kdTree Search
    //  create a selector for each thread that is as disparate as possible from the other thread.
    int64_t numThreadsMSB = static_cast<int64_t>(ceill(log2((float)numThreads)))-1;
    uint64_t selectionBias;
    if (((static_cast<int64_t>(0x1) << numThreadsMSB) & threadNum) == 0) {
      selectionBias = threadNum;
    }
    else {
      selectionBias = (0xFFFFFFFFFFFFFFFF << numThreadsMSB) | threadNum;
    }
    Cluster* newCluster = new Cluster(IDCounter, IDCounterMutex);  // create a new cluster
    while (kdTree->pickValue(retPair, selectionBias, newCluster)) {  // pick an initial cluster point until no more points
      keys->clear();
      newCluster->addToCluster(retPair.first, retPair.second); // add the picked Point to the cluster to the cluster
      keys->push_back(retPair.first);              // add the point to the keys list
      while (!keys->empty()) {                     // search around each point in the key list until empty
        K* center = keys->front();
        keys->pop_front();
        for (size_t i = 0; i < N; i++) {           // build the search bounds
          qLower.at(i) = center[i] - clusterWindowRadius[i];
          qUpper.at(i) = center[i] + clusterWindowRadius[i];
        }
        retPairs.clear();
        kdTree->searchRegionAndRemove(retPairs, qLower, qUpper, newCluster); //search the tree for points within the window
        for (retPair_t* kit : retPairs) { // add the points to the keys list and cluster
          keys->push_back(kit->first);
          newCluster->addToCluster(kit->first, kit->second);
          delete kit;
        }
      }
      if (newCluster->values->size() != 0)  {  // if newCluster has non-zero length values, add it to the clusters list
        clusters->push_back(newCluster);
        if (newCluster->overlaps != nullptr) {  // if newCluster has overlaps, add it to the overlap list.
          overlaps->push_back(newCluster);
        }
      }
      newCluster = new Cluster(IDCounter, IDCounterMutex);  // create a new cluster
      if (oneTime) break;
    }
    delete keys;
    return true;
  }

  // overlapMerge() traces the clusters that are linked to initialCluster during the buildCluster 
  // phase.  Each time it finds a non-zero cluster in combines that data into initialCluster.  Then
  // looks for any cluster linked to that cluster by recurring on it's overlap list.
  // if the cluster being inspected is a zero-length cluster, that means that its's data and the data in
  // all the clusters below it have already been copied so it can be skipped.
  void overlapMerge(avlSet<Cluster*>* overlaps, Cluster* initialCluster, int depth) {
#ifdef DEBUG_PRINT
    int cnt = 0;
#endif
    std::vector<Cluster*> overlapClusters(overlaps->size());
    overlaps->getKeys(overlapClusters);
    for (Cluster* thisOverlap : overlapClusters) { // iterate through the overlaps.
    
#ifdef DEBUG_PRINT
      cnt++;
      int64_t sz = thisOverlap->values->size();
      std::cout << "At cnt " << cnt << " depth " << depth << ", cluster " << thisOverlap->getID() << " has length " << sz << "\n";
#endif
      if (thisOverlap->values->size() == 0)   // if this cluster is zero length, it's data has already been copied
        continue; // so move on.
      if (thisOverlap == initialCluster) // check to see if this is pointing to initialCLuster
        continue; // if so, move on
#ifdef DEBUG_PRINT
      std::cout << "Merging cluster " << thisOverlap->getID() << " into cluster " << initialCluster->getID() << "\n";
#endif
      initialCluster->combine(thisOverlap);  // combine this overlap, clusters data into initialCluster
      thisOverlap->values->clear();          // and zero it's length
      if (depth >= numThreads * 10) {          // maximum number of search recursions should not exceed numThreads x 10
        throw std::runtime_error("Internal ERROR: maximum possible number of linked overlapped clusters exceeded on cluster ");
      }
      if (thisOverlap->overlaps != nullptr)
        overlapMerge(thisOverlap->overlaps, initialCluster, depth + 1); // go to the next level
    }
    return;
  }


public:
  // base constructor
  MTDBSCAN() {
    clusters = new std::vector<Cluster*>;
    kdTree = new KdTreePDB<K, V, Cluster, N>();
    setNumThreads();
  }
  // constructor with window
  MTDBSCAN(std::vector<K>& window) {
    clusters = new std::vector<Cluster*>;
    kdTree = new KdTreePDB<K, V, Cluster, N>();
    setWindow(window);
  }

  ~MTDBSCAN() {
    for (size_t i = 0; i < clusters->size(); i++) {
      delete clusters->at(i);
    }
    delete clusters;
  }

  // Set the number of threads used for building the clusters.  The default is to use the hardware concurrency
  void setNumThreads(int64_t thd = -1) {
    if (thd == -1) numThreads = std::thread::hardware_concurrency();
    else numThreads = thd;
    kdTree->setNumThreads(numThreads);
  }

  // addPoint fuction adds a point or key and an associated value to the MTDBSCAN object
  size_t addPointWithValue(std::vector<K> point, V value) {
    numValues = kdTree->add(point, value);
    return numValues;
  }

  // setWindow sets the clustering widnow size
  void setWindow(std::vector<K>& window) {
    clusterWindowRadius = window;
  }

#ifdef MEASURE_KDTREE_TIME
  std::chrono::steady_clock::time_point kdtBefore;
  std::chrono::steady_clock::time_point kdtAfter;
  // get the current time as a double.  Used for evaluating performance.
  double getKdTreeTime() {
    return std::chrono::duration<double>(kdtAfter - kdtBefore).count();
  }
#endif

  // build() builds the clusters.  All points and values must be added first.
  // It will return false if called a second time after MTDBSCAN was constructed.
  bool build() {
#ifdef MEASURE_KDTREE_TIME
    kdtBefore = std::chrono::steady_clock::now();
#endif
    if (kdTree != nullptr && kdTree->buildTree() == false) return false;  // build a kdTree using default threads.
#ifdef MEASURE_KDTREE_TIME
    kdtAfter = std::chrono::steady_clock::now();
#endif

    bool returnStatus = true;
    std::vector<std::future<bool>> futures;
    std::vector<std::list<Cluster*>*> clustersPerThread;
    std::vector<std::list<Cluster*>*> overlapClusters;

    // Create a temporary cluster list for each thread.
    // Add one cluster to each befor starting multithreading.  This is
    // done to mitigate the performance pothole when the search region is large
    // relative to the the region of the data.
    for (int i = 0; i < numThreads; i++) {
      clustersPerThread.push_back(new std::list<Cluster*>);
      overlapClusters.push_back(new std::list<Cluster*>);
      //returnStatus = returnStatus && buildCluster( clustersPerThread[i], overlapClusters[i], i, true);
    }

     //Now start start the individual threads.
    for (int i = 0; i < numThreads; i++) {
      futures.push_back(std::async(std::launch::async,
        &MTDBSCAN<K, V, N>::buildCluster, this, clustersPerThread[i], overlapClusters[i], i, false));
    }

    // get the results back from each thread.
    for (int i = 0; i < futures.size(); i++) 
      returnStatus = returnStatus && futures[i].get();

    if (returnStatus == true) {


    // loop through all of the overlapsClusters an make sure there is a bidirectional link to every cluster
/*      for (int i = 0; i < numThreads; i++) {  // iterate through per thread cluster lists
        //std::cout << "cl = " << clustersPerThread[i]->size() << " ov = " << overlapClusters[i]->size() << std::endl;
        for (Cluster* cl : *overlapClusters[i]) { // iterate through the overlap clusters
          for (Cluster* cl1 : *cl->overlaps) {
            if (cl1->overlaps == nullptr) {
              cl1->overlaps = new std::set<Cluster*>;
            }
            cl1->overlaps->insert(cl);
          }
        }
      }
      */
    // if there were no errors in the threaded clustering,
    // combine any clusters that were found to have overlapping points during the clustering process
      //size_t clcnt = 0, ovcnt = 0;
      for (int i = 0; i < numThreads; i++) {  // iterate through per thread cluster lists
       // clcnt += clustersPerThread[i]->size();
        //ovcnt += overlapClusters[i]->size();
        //std::cout << "cl = " << clustersPerThread[i]->size() << " ov = " << overlapClusters[i]->size() << std::endl;
        for (Cluster* cl : *clustersPerThread[i]) { // iterate through the overlap clusters
          if (cl->overlaps != nullptr && cl->values->size() != 0) {  // Check to see if are non-zero-length overlaps
            try {
              overlapMerge(cl->overlaps, cl, 0);  // combine those overlaps into this cluster
            }
            catch (std::runtime_error& e) {
              std::cerr << e.what() << std::endl;
              exit(1);
            }
          }
        }
      }
      //std::cout << "cl = " << clcnt << " ov = " << ovcnt << std::endl;

      // Move the contents temporary per-thread cluster lists to the final cluster vector.
      // Exclude clusters and delete those clusters that have no values.
      for (int i = 0;  i < numThreads; i++) {
        for (Cluster* cl : *clustersPerThread[i]) {
          if (cl->overlaps != nullptr) {
            delete cl->overlaps;
            cl->overlaps = nullptr;
          }
          if (cl->values->size() != 0) clusters->push_back(cl);
          else delete cl;
        }
      }
    }
    else {
      std::cerr << "ERROR: an error occurred during the cluster build process." << std::endl;
    }

    // Clean up the temporary lists.
    for (int i = 0; i < numThreads; i++) {
      delete overlapClusters[i];
      delete clustersPerThread[i];
    }
    // and the temporary KD Tree
    delete kdTree;
    kdTree = nullptr;
    return returnStatus;
  }

  // sortClusterBySize sorts the cluster by number of values in the cluster 
  // from largest to smallest
  void sortClustersBySize(bool smallestFirst = false) {
    if (!smallestFirst) {
      parallelSort(clusters->begin(), clusters->end(),
        [](const Cluster* a, const Cluster* b) -> bool
        {
          return a->values->size() > b->values->size();
        });
    }
    else {
      parallelSort(clusters->begin(), clusters->end(),
        [](const Cluster* a, const Cluster* b) -> bool
        {
          return a->values->size() < b->values->size();
        });
    }
  }

 // sortClusterByFirstValue sorts the cluster by the first value in the cluster 
 // from largest to smallest.  This is used to accelerate comparison of cluster algorithms.
  void sortClustersByFirstValue() {
    parallelSort(clusters->begin(), clusters->end(),
      [](const Cluster* a, const Cluster* b) -> bool
      {
        return a->values->front() > b->values->front();
      });
  }

  // cluster setters and getters
  // getClusterSize returns the number of values in cluster clusterIdx
  size_t getClusterSize(size_t clusterIdx) {
    if (clusterIdx >= clusters->size()) return 0;
    return clusters->at(clusterIdx)->values->size();
  }
  // getClusterValueList returns the list of values in cluster clusterIdx
  std::list<V>* getClusterValueList(size_t clusterIdx) {
    if (clusterIdx >= clusters->size()) return nullptr;
    return (clusters->at(clusterIdx))->values;
  }
  // getClusterMaxCorner returns the a vector of the maximum corner of the
  // bounding hypercube of cluster clusterIdx
  std::vector<K>* getClusterMaxCorner(size_t clusterIdx) {
    if (clusterIdx >= clusters->size()) return nullptr;
    Cluster* cluster = clusters->at(clusterIdx);
    return &(cluster->clusterBounds.min);
  }
  // getClusterMinCorner returns the a vector of the minimum corner of the
  // bounding hypercube of cluster clusterIdx
  std::vector<K>* getClusterMinCorner(size_t clusterIdx) {
    if (clusterIdx >= clusters->size()) return nullptr;
    Cluster* cluster = clusters->at(clusterIdx);
    return &(cluster->clusterBounds.max);
  }
  // getClusterTag returns the tag string set by setClusterTag
  std::string* getClusterTag(size_t clusterIdx) {
    if (clusterIdx >= clusters->size()) return nullptr;
    return clusters[clusterIdx]->tag;
  }

  // getClusterID returns the ID assigned when cluster was created.
  uint32_t getClusterID(size_t clusterIdx) {
    if (clusterIdx >= clusters->size()) return 0;
    Cluster* cluster = clusters->at(clusterIdx);
    return cluster->getID();
  }

  // cluster tag setter
  bool setClusterTag(size_t clusterIdx, std::string& tagIn) {
    if (clusterIdx >= clusters->size()) return false;
    clusters[clusterIdx]->tag = new std::string(tagIn);
    return true;
  }


  // getNumClusters returns the number of clusters
  size_t getNumClusters() {
    return clusters->size();
  }

  // checkCluster does some basic check on the clusters and print some statistics
  bool checkClusters(const size_t numLocations, const std::string* tag = nullptr) {
    size_t max = std::numeric_limits<size_t>::min();
    size_t min = std::numeric_limits<size_t>::max();
    size_t avgSize = 0;
    size_t count = 0;
    bool rb = true;
    for (size_t i = 0; i < clusters->size(); i++) {
      Cluster* c = clusters->at(i);
      //      if (tag == nullptr || c.tag != nullptr && c.hasTag(tag)) {
      if (tag == nullptr) {
        size_t t = c->values->size();
        if (t == 0) {
          std::cout << "Cluster " << i << " has 0 entries" << std::endl;
        }
        if (t > max) max = t;
        if (t < min) min = t;
        avgSize += t;
        count++;
      }
    }
    if (numLocations != -1 && avgSize != numLocations) {
      std::cout << "Number of locations in all clusters not equal input locations." << std::endl;
      rb = false;
    }
    avgSize = (count > 0) ? avgSize / count : 0;
    min = (count > 0) ? min : 0;
    max = (count > 0) ? max : 0;
    std::string t_tag = (tag == nullptr ? "" : *tag + " ");
    std::cout << "Cluster " << t_tag << "Count = " << count << " Max = " << max <<
      "  Min = " << min << " Average = " << avgSize << std::endl;
    return rb;
  }
};

#endif // MTDBSCAN_HPP
