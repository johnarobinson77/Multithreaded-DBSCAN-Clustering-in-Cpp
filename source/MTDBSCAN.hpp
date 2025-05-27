/**
 * Copyright (c) 2023 John Robinson.
 * Copyright (c) 2024 John Robinson.
 * Copyright (c) 2025 John Robinson.
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
#include "avlSet.hpp"

template <typename K, typename V, size_t N>
class MTDBSCAN
{
  // This is a shortened name of an std::pair that contains a tuple in the fist spot and a value in the second
  // When the KdTree search functions are called, a list of these pairs are returned.
  typedef std::pair<K *, std::list<V> *> retPair_t;

  /// @brief number of point,value pairs added to the cluster.
  int64_t numThreads = -1; 
  /// @brief number of point,value pairs added to the cluster.
  int64_t numValues = 0; 
  
  // Each cluster must have a unique ID.  This counter is incremented each time a new cluster is
  // created.  Since clusters are created in multiple threads, the increment process must be mutex locked.
  /// @brief counter used to assign a unique ID
  uint32_t IDCounter = 1;
  /// @brief mutex used for making sure counter increments correctly in multithreading environment.
  std::mutex IDCounterMutex;

  /// @brief his struct maintains a bounding hyper rectangle
  struct Bounds
  {
    /// @brief vector containing the max bounds
    std::vector<K> max;
    /// @brief vector containing the min bounds
    std::vector<K> min;
    
    /// @brief replaces the min and max bounds with nv, component by component
    ///  if smaller or bigger respectively.  nv must have the same number of components as the bounds.
    /// @param nv vector of values to replace the bounds.
    void addToBounds(std::vector<K> &nv)
    {
      for (size_t i = 0; i < max.size(); i++)
      {
        max[i] = (nv[i] > max[i]) ? nv[i] : max[i];
        min[i] = (nv[i] < min[i]) ? nv[i] : min[i];
      }
    }
    
    /// @brief replaces the min and max bounds with nv, component by component
    ///  if smaller or bigger respectively.  nv must have the same number of components as the bounds.
    /// @param nv pointer to array of values to replace the bounds.
    void addToBounds(K nv[])
    {
      for (size_t i = 0; i < max.size(); i++)
      {
        max[i] = (nv[i] > max[i]) ? nv[i] : max[i];
        min[i] = (nv[i] < min[i]) ? nv[i] : min[i];
      }
    }
    
    /// @brief Replaces the min and max bounds with addBounds, component by component
    ///  if smaller or bigger respectively.  nv must have the same number of components as the bounds.
    /// @param nv pointer to array of values to replace the bounds.
    void addToBounds(Bounds &addBounds)
    {
      addToBounds(addBounds.max);
      addToBounds(addBounds.min);
    }

    /// @brief Sets the min bound to max and the max bounds to min so that the next addToBounds will replace
    ///  both.
    void resetBounds()
    {
      for (size_t i = 0; i < max.size(); i++)
      {
        max[i] = std::numeric_limits<K>::lowest();
        min[i] = std::numeric_limits<K>::max();
      }
    }
    
    /// @brief resizes the number of components in bounds and sets the min bound to max and the max
    ///  bounds to min so that the next addToBounds will replace both.
    /// @param newSize
    void resetBounds(int newSize)
    {
      max.resize(newSize);
      min.resize(newSize);
      resetBounds();
    }

    /// @brief equals checks to see that two bounds are equal
    /// @param otherBounds reference to the other bounds
    /// @return compares two bounds to check for equality.
    bool equals(Bounds &otherBounds)
    {
      bool eq = true;
      for (size_t i = 0; i < max.size(); i++)
      {
        if (max[i] != otherBounds.max[i])
          eq = false;
        if (min[i] != otherBounds.min[i])
          eq = false;
      }
      return eq;
    }

    /// @brief constructor for bounds
    /// @param size number of components in bounds
    Bounds(size_t size)
    {
      max.resize(size);
      min.resize(size);
      resetBounds();
    }

    Bounds() = delete;
  }; // bounds
  typedef Bounds bounds_t;


  /// @brief This struct sums the average value of the points associated with the values in a cluster and then
  /// divides by the number of values to approximate the midpoint.
  struct Average {
    std::vector<K> sum;
    uint64_t count;

    /// @brief add a point to be averaged
    /// @param in the point to be averaged in
    inline void addToAverage(std::vector<K>& in) {
      for (size_t i = 0;  i < sum.size(); i++) sum[i] += in[i];
      count++;
    }

    /// @brief add a point to be averaged
    /// @param in the point to be averaged in
    inline void addToAverage(K* in) {
      for (size_t i = 0;  i < sum.size(); i++) sum[i] += in[i];
      count++;
    }

    /// @brief add currently running Average
    /// @param in the Average to be averaged in
    inline void addToAverage(Average& in) {
      for (size_t i = 0;  i < sum.size(); i++) sum[i] += in.sum[i];
      count += in.count;
    }

    /// @brief reset the average 
    inline void reset() {
      for (size_t i = 0;  i < sum.size(); i++) sum[i] = 0;
      count = 0;
    }

    /// @brief return the average 
    /// @return a reference to the average of all the points added.  0 if not points have been added.
    inline std::vector<K>* get() {
      if (count == 0) reset();
      if (count != 1) {
        for (size_t i = 0;  i < sum.size(); i++) sum[i] /= count;
        count = 1;
      }
      return &sum;
    }

    Average(size_t n) {
      sum.resize(n);
      reset();
    }
  };

  // Objects of this struct contain the list of values of a single cluster as will as bounds and a tag.
  struct Cluster
  {

    /// @brief Pointer to the value list
    std::list<V> *values = nullptr; 
    /// @brief the min and max corner of a bounding hyper rectangle
    bounds_t clusterBounds{N};
    /// @brief the average of the points associated with the values in the cluster
    Average clusterAverage{N};  
    /// @brief pointer to an optional tag string
    std::string *tag = nullptr; 
    /// @brief pointer to a sent of pointers to clusters that this cluster overlaps with. Used only during build
    avlSet<Cluster *> *overlaps = nullptr;
    /// @brief ID of this cluster.  Only used during build.
    uint32_t ID = 0; 

    /// @brief Make sure the default constructor isn't called.
    Cluster() = delete;

    /// @brief COnstructor for a cluster.
    /// @param IDCounter
    /// @param IDCounterMutex
    Cluster(uint32_t &IDCounter, std::mutex &IDCounterMutex)
    {
      values = new std::list<V>();
      tag = nullptr;
      { // mutex lock this operation.
        std::lock_guard<std::mutex> guard(IDCounterMutex);
        ID = IDCounter++;
      }
    }

    /// @brief getID retruns the ID of the cluster
    /// @return the ID
    inline uint32_t getID() { return ID; }

    /// @brief addToCluster copies a llist of values to the cluster and upsdates the bounds.
    /// @param kdKey pointer to a point or tuple associated with the values in vals.  Used to update the bounds
    /// @param vals list of values to add to the cluster.
    /// @return
    bool addToCluster(K *kdKey, std::list<V, STLAdaptor<V>> &vals)
    {
      for (V val : vals)
      {
        values->push_back(val);
      }
      clusterBounds.addToBounds(kdKey);
      clusterAverage.addToAverage(kdKey);
      return true;
    }

    /// @brief Combine combines the cluster pointed to by ocl into this cluster
    /// @param ocl pointer to the other cluster
    /// @return true if no failure.
    bool combine(Cluster *ocl)
    {
      values->splice(this->values->end(), *ocl->values);
      clusterBounds.addToBounds(ocl->clusterBounds);
      clusterAverage.addToAverage(ocl->clusterAverage);
      return true;
    }

    /// @brief destructor.  I deletes the values list and tag and overlap list
    ~Cluster()
    {
      if (values != nullptr)
        delete values;
      if (tag != nullptr)
        delete tag;
      if (overlaps != nullptr)
        delete overlaps;
    }
  };

  /// @brief this is the list of clusters.
  std::vector<Cluster *> clusters;
  /// @brief Pointer to the temporary KdTree that will be used to accelerated the clustering.
  KdTreePDB<K, V, Cluster, N> *kdTree;
  /// @brief this vector contains the distance to the edge of the search window in each dimension.
  std::vector<K> clusterWindowRadius;

  /// @brief buildClusterList repeatedly builds clusters until all points in the kdTree have been used
  /// @param clusters The list of clusters built by a call to this function
  /// @param overlaps The list of clusters built by this function that have overlaps with other clusters
  /// @param threadNum A unique thread indicator.  Should be part of a series of numbers between 0 and thread - 1.
  /// @return true if no failures
  bool buildClusterList(std::list<Cluster *> *clusters, std::list<Cluster *> *overlaps, int64_t threadNum)
  {
    std::vector<K> qLower(N); // allocate 2 vectors for the upper and lower corners of the window
    std::vector<K> qUpper(N);
    SimpleList<K *> keys; // holds the list of the points to be searched for a cluster.
    //  create a selector for each thread that is as disparate as possible from the other thread.
    int64_t numThreadsMSB = static_cast<int64_t>(ceill(log2((float)numThreads))) - 1;
    uint64_t selectionBias;
    if (((static_cast<int64_t>(0x1) << numThreadsMSB) & threadNum) == 0)
    {
      selectionBias = threadNum;
    }
    else
    {
      selectionBias = (0xFFFFFFFFFFFFFFFF << numThreadsMSB) | threadNum;
    }
    Cluster *newCluster = new Cluster(IDCounter, IDCounterMutex); // create a new cluster
    keys.reset();
    while (kdTree->pickValue(keys, selectionBias, newCluster))
    { // pick an initial cluster point until no more points
      while (!keys.empty())
      { // search around each point in the key list until empty
        K *center;
        keys.pop_front(center);
        for (size_t i = 0; i < N; i++)
        { // build the search bounds
          qLower[i] = center[i] - clusterWindowRadius[i];
          qUpper[i] = center[i] + clusterWindowRadius[i];
        }
        kdTree->searchRegionAndRemove(keys, qLower, qUpper, newCluster); // search the tree for points within the window
      }
      if (newCluster->values->size() != 0)
      { // if newCluster has non-zero length values, add it to the clusters list
        clusters->push_back(newCluster);
        if (newCluster->overlaps != nullptr)
        { // if newCluster has overlaps, add it to the overlap list.
          overlaps->push_back(newCluster);
        }
      }
      else
      {
        delete newCluster;
      }
      newCluster = new Cluster(IDCounter, IDCounterMutex); // create a new cluster
      keys.reset();                                        // reset the keys FIFO.
    }
    if (newCluster->values->size() == 0)
      delete newCluster;
    return true;
  }

  /// @brief overlapMerge() traces the clusters that are linked to initialCluster during the buildCluster
  // phase.  Each time it finds a non-zero cluster in combines that data into initialCluster.  Then
  // looks for any cluster linked to that cluster by recurring on it's overlap list.
  // if the cluster being inspected is a zero-length cluster, that means that its's data and the data in
  // all the clusters below it have already been copied so it can be skipped.
  /// @param overlaps set containing any clusters that were found to overlap with initial cluster
  /// @param initialCluster the initial cluster.
  /// @param depth current dept of the recursion
  void overlapMerge(avlSet<Cluster *> *overlaps, Cluster *initialCluster, int depth)
  {
#ifdef DEBUG_PRINT
    int cnt = 0;
#endif
    std::vector<Cluster *> overlapClusters(overlaps->size());
    overlaps->getKeys(overlapClusters);
    for (Cluster *thisOverlap : overlapClusters)
    { // iterate through the overlaps.

#ifdef DEBUG_PRINT
      cnt++;
      int64_t sz = thisOverlap->values->size();
      std::cout << "At cnt " << cnt << " depth " << depth << ", cluster " << thisOverlap->getID() << " has length " << sz << "\n";
#endif
      if (thisOverlap->values->size() == 0) // if this cluster is zero length, it's data has already been copied
        continue;                           // so move on.
      if (thisOverlap == initialCluster)    // check to see if this is pointing to initialCLuster
        continue;                           // if so, move on
#ifdef DEBUG_PRINT
      std::cout << "Merging cluster " << thisOverlap->getID() << " into cluster " << initialCluster->getID() << "\n";
#endif
      initialCluster->combine(thisOverlap); // combine this overlap, clusters data into initialCluster
      thisOverlap->values->clear();         // and zero it's length
      if (depth >= numThreads * 10)
      { // maximum number of search recursions should not exceed numThreads x 10
        throw std::runtime_error("Internal ERROR: maximum possible number of linked overlapped clusters exceeded on cluster ");
      }
      if (thisOverlap->overlaps != nullptr)
        overlapMerge(thisOverlap->overlaps, initialCluster, depth + 1); // go to the next level
    }
    return;
  }

public:
  /// @brief default constructor
  MTDBSCAN() noexcept 
  {
    kdTree = new KdTreePDB<K, V, Cluster, N>();
    setNumThreads();
  }

  /// @brief constructor with window
  MTDBSCAN(std::vector<K> &window) noexcept 
  {
    clusters = new std::vector<Cluster *>;
    kdTree = new KdTreePDB<K, V, Cluster, N>();
    setWindow(window);
  }

  /// @brief destructor.  Deletes the entire list of clusters.
  ~MTDBSCAN() noexcept 
  {
    for (size_t i = 0; i < clusters.size(); i++)
    {
      delete clusters.at(i);
    }
    deleteKdTree();
  }

  /// @brief setNumThreads sets the number of threads sed to build the kdTree and then build the cluster list.
  /// @param thd number of threads.  The default is to use the hardware concurrency.
  void setNumThreads(int64_t thd = -1)  noexcept 
  {
    if (thd == -1)
      numThreads = std::thread::hardware_concurrency();
    else
      numThreads = thd;
    kdTree->setNumThreads(numThreads);
  }

  /// @brief addPointWithValue adds a point to be clustered and the value associated with that point
  /// @param point vector of a point.  should have the same dimensions as specified in MTDBSCAN
  /// @param value value associated with the point
  /// @return returns the number of points added so far.
  size_t addPointWithValue(std::vector<K> &point, V value) noexcept 
  {
    numValues = kdTree->add(point, value);
    return numValues;
  }

  /// @brief setWindow sets the search window used to build a cluster
  /// @param window is a vector distances from the center of some search point to include new points
  ///               It must have the same dimensions as as the N in MTDBSCAN
  void setWindow(std::vector<K> &window)
  {
    clusterWindowRadius = window;
  }

#ifdef MEASURE_INDIVIDUAL_TIMES
  std::chrono::steady_clock::time_point kdtBefore;
  std::chrono::steady_clock::time_point kdtAfter;
  std::chrono::steady_clock::time_point cstrBefore;
  std::chrono::steady_clock::time_point cstrAfter;
  std::chrono::steady_clock::time_point mrgBefore;
  std::chrono::steady_clock::time_point mrgAfter;
  std::chrono::steady_clock::time_point cpyBefore;
  std::chrono::steady_clock::time_point cpyAfter;
  std::chrono::steady_clock::time_point cleanupBefore;
  std::chrono::steady_clock::time_point cleanupAfter;
  // return kdTree build time..
  double getKdTreeTime()
  {
    return std::chrono::duration<double>(kdtAfter - kdtBefore).count();
  }
  // return time to search the kd tree for all clusters and overlap clusters.
  double getSearchTime()
  {
    return std::chrono::duration<double>(cstrAfter - cstrBefore).count();
  }
  // return the time to merge overlap clusters
  double getMergeTime()
  {
    return std::chrono::duration<double>(mrgAfter - mrgBefore).count();
  }
  // return the time to copy cluster lists to the final vector
  double getCopyTime()
  {
    return std::chrono::duration<double>(cpyAfter - cpyBefore).count();
  }
  // return the time to delete intermediate data.
  double getCleanupTime()
  {
    return std::chrono::duration<double>(cleanupAfter - cleanupBefore).count();
  }
#endif

  /// @brief buildKdTree builds a kdTree using the points added by addPointWithValue.
  /// @return true is no errors
  bool buildKdTree() 
  {
#ifdef MEASURE_INDIVIDUAL_TIMES
    kdtBefore = std::chrono::steady_clock::now();
#endif
    if (kdTree != nullptr && kdTree->isBuilt())
    {
      std::cout << "KdTree has already been built. \n";
      return false;
    }

    if (kdTree != nullptr && kdTree->buildTree() == false)
    {
      std::cout << "KdTree build error\n";
      return false;
    }
#ifdef MEASURE_INDIVIDUAL_TIMES
    kdtAfter = std::chrono::steady_clock::now();
#endif
    return true;
  }


  /// @brief deleteKdTree deletes the kdTree used to build a set of clusters.
  /// After the clusters are built it has no function.
  void deleteKdTree()
  {
    if (kdTree != nullptr)
      delete kdTree;
    kdTree = nullptr;
  }


  /// @brief builds a set of clusters using an exiting kdTree and the current window settings.
  /// @return true if no errors.
  bool buildClusters()
  {
    // Create a temporary cluster list for each thread.
#ifdef MEASURE_INDIVIDUAL_TIMES
    cstrBefore = std::chrono::steady_clock::now();
#endif
    bool returnStatus = true;
    std::vector<std::future<bool>> futures;
    std::vector<std::list<Cluster *> *> clustersPerThread;
    std::vector<std::list<Cluster *> *> overlapClusters;
    for (int i = 0; i < numThreads; i++)
    {
      clustersPerThread.push_back(new std::list<Cluster *>);
      overlapClusters.push_back(new std::list<Cluster *>);
    }

    // totalClusterCount is used accumulate total cluster count
    size_t totalClusterCount = 0ULL;

    // Start building the clusters running the DBSCAN search using the specified number of threads.
    for (int i = 0; i < numThreads; i++)
    {
      futures.push_back(std::async(std::launch::async,
                                   &MTDBSCAN<K, V, N>::buildClusterList, this, clustersPerThread[i], overlapClusters[i], i));
    }

    // get the results back from each thread.
    for (int i = 0; i < futures.size(); i++)
    {
      returnStatus = returnStatus && futures[i].get();
      totalClusterCount += clustersPerThread[i]->size();
    }
#ifdef MEASURE_INDIVIDUAL_TIMES
    cstrAfter = std::chrono::steady_clock::now();
#endif

    if (returnStatus == true)
    {
      // if there were no errors in the threaded clustering,
      // merge any clusters that were found to have overlapping points during the clustering process
#ifdef MEASURE_INDIVIDUAL_TIMES
      mrgBefore = std::chrono::steady_clock::now();
#endif
      for (int i = 0; i < numThreads; i++)
      { // iterate through per thread cluster lists
        for (Cluster *cl : *overlapClusters[i])
        { // iterate through the overlap clusters
          if (cl->overlaps != nullptr && cl->values->size() != 0)
          { // Check to see if are non-zero-length overlaps
            try
            {
              overlapMerge(cl->overlaps, cl, 0); // combine those overlaps into this cluster
            }
            catch (std::runtime_error &e)
            {
              std::cerr << e.what() << std::endl;
              exit(1);
            }
          }
        }
      }
#ifdef MEASURE_INDIVIDUAL_TIMES
      mrgAfter = std::chrono::steady_clock::now();
#endif

#ifdef MEASURE_INDIVIDUAL_TIMES
      cpyBefore = std::chrono::steady_clock::now();
#endif
      // first remove overlap data since it's no longer needed.
      // Of course it would only be the clusters on the overlapClusters lists
      for (int i = 0; i < numThreads; i++)
      { // iterate through per thread cluster lists
        for (Cluster *cl : *overlapClusters[i])
        { // iterate through the overlap clusters
          if (cl->overlaps != nullptr)
          {
            delete cl->overlaps;
            cl->overlaps = nullptr;
          }
        }
      }

      // Move the contents temporary per-thread cluster lists to the final cluster vector.
      // Exclude and delete those clusters that have no values.
      // Start by reserving enough space in the vector to include all clusters
      clusters.reserve(totalClusterCount);
      for (int i = 0; i < numThreads; i++)
      {
        size_t j = 0;
        for (Cluster *cl : *clustersPerThread[i])
        {
          if (cl->values->size() != 0)
            clusters.push_back(cl);
          else
            delete cl;
        }
      }
      // Since not all clusters were moved to the clusters vector, shrink it to the remaining number
      clusters.shrink_to_fit();
#ifdef MEASURE_INDIVIDUAL_TIMES
      cpyAfter = std::chrono::steady_clock::now();
#endif
    }
    else
    {
      std::cerr << "ERROR: an error occurred during the cluster build process.\n";
      return false;
    }
    // delete the per tread lists.
    for (int i = 0; i < numThreads; i++)
    {
      delete overlapClusters[i];
      delete clustersPerThread[i];
    }
    return true;
  }

  /// @brief build() builds the clusters.  All points and values must be added first.
  /// @return t will return false if called a second time after MTDBSCAN was constructed.
  bool build()
  {
    // First build the KdTree from the points provided with the addPointWithValue function
    if (buildKdTree() == false)
      return false;

    // build clusters.
    if (buildClusters() == false)
      return false;

#ifdef MEASURE_INDIVIDUAL_TIMES
    cleanupBefore = std::chrono::steady_clock::now();
#endif

    // delete temporary KD Tree
    deleteKdTree();
#ifdef MEASURE_INDIVIDUAL_TIMES
    cleanupAfter = std::chrono::steady_clock::now();
#endif

    return true;
  }

  /// @brief sortClusterBySize sorts the clusters by number of values in the cluster
  /// @param smallestFirst true is smallest to largest , false is largest to smallest (default).
  void sortClustersBySize(bool smallestFirst = false) noexcept
  {
    if (!smallestFirst)
    {
      parallelSort(clusters.begin(), clusters.end(), [](const Cluster *a, const Cluster *b) -> bool
                   { return a->values->size() > b->values->size(); }, 4);
    }
    else
    {
      parallelSort(clusters.begin(), clusters.end(), [](const Cluster *a, const Cluster *b) -> bool
                   { return a->values->size() < b->values->size(); }, 4);
    }
  }

  /// @brief sortClusterByFirstValue sorts the cluster by the first value in the cluster.  Used in test for comparing 2 Cluster sets
  void sortClustersByFirstValue() noexcept 
  {
    parallelSort(clusters.begin(), clusters.end(), [](const Cluster *a, const Cluster *b) -> bool
                 { return a->values->front() > b->values->front(); }, 4);
  }


  /// @brief getClusterSize returns the size of the cluster by index
  /// @param clusterIdx index of the cluster 
  /// @return size of the cluster 0r 0 if the cluster index is out of range.
  size_t getClusterSize(size_t clusterIdx) const noexcept
  {
    if (clusterIdx >= clusters.size())
      return 0;
    return clusters.at(clusterIdx)->values->size();
  }


  /// @brief getClusterValueList returns a pointer to the list of values in the cluster at clusterIdx
  /// @param clusterIdx index of the cluster
  /// @return pointer to a list of values in that cluster or nullptr if clusterIdx is out of range.
  std::list<V> *getClusterValueList(size_t clusterIdx) const noexcept
  {
    if (clusterIdx >= clusters.size())
      return nullptr;
    return (clusters.at(clusterIdx))->values;
  }


  /// @brief getClusterMaxCorner returns a pointer to a vector of the maximum corner of the
  /// bounding hypercube of cluster clusterIdx
  /// @param clusterIdx index of the cluster
  /// @return pointer to the max bounds vector or nullptr if clusterIdx is out of range.
  std::vector<K> *getClusterMaxCorner(size_t clusterIdx) const noexcept
  {
    if (clusterIdx >= clusters.size())
      return nullptr;
    Cluster *cluster = clusters.at(clusterIdx);
    return &(cluster->clusterBounds.min);
  }

  
  /// @brief getClusterMinCorner returns a pointer to a vector of the minimum corner of the
  /// bounding hypercube of cluster clusterIdx
  /// @param clusterIdx index of the cluster
  /// @return pointer to the min bounds vector or nullptr if clusterIdx is out of range.
  std::vector<K> *getClusterMinCorner(size_t clusterIdx) const noexcept
  {
    if (clusterIdx >= clusters.size())
      return nullptr;
    Cluster *cluster = clusters.at(clusterIdx);
    return &(cluster->clusterBounds.max);
  }

  
  /// @brief getClusterCenter returns a vector of the point that is the average of the 
  /// points associated with the cluster at clusterIdx
  /// @param clusterIdx index of the cluster
  /// @return pointer to the center (average) vector, or nullptr if clusterIdx is out of range.
  std::vector<K> *getClusterCenter(size_t clusterIdx) noexcept
  {
    if (clusterIdx >= clusters.size())
      return nullptr;
    Cluster *cluster = clusters.at(clusterIdx);
    return cluster->clusterAverage.get();
  }

  
  /// @brief getClusterID returns the ID assigned when cluster was created.  
  /// IDs are unique and greater than 0 but not sequential
  /// @param clusterIdx index of the cluster
  /// @return cluster ID or 0 if clusterIdx is outside the range
  uint32_t getClusterID(size_t clusterIdx) const noexcept
  {
    if (clusterIdx >= clusters.size())
      return 0;
    Cluster *cluster = clusters.at(clusterIdx);
    return cluster->getID();
  }

  
  /// @brief setClusterTag sets a tag in the cluster at clusterIdx in the forma of an std::string
  /// @param clusterIdx index of the cluster
  /// @param tagIn string to copy to the cluster tag
  /// @return true if successful, false if clusterIdx is outside the range.
  bool setClusterTag(size_t clusterIdx, std::string &tagIn) noexcept
  {
    if (clusterIdx >= clusters.size())
      return false;
    clusters[clusterIdx]->tag = new std::string(tagIn);
    return true;
  }

  
  /// @brief getClusterTag returns the tag string set by setClusterTag in the cluster at clusterIdx
  /// @param clusterIdx index of the cluster
  /// @return pointer to the tag string or nullptr if clusterIdx is outside the range.
  std::string *getClusterTag(size_t clusterIdx) const noexcept
  {
    if (clusterIdx >= clusters.size())
      return nullptr;
    return clusters[clusterIdx]->tag;
  }

  
  /// @brief getNumClusters returns the number of clusters
  /// @return the number of clusters
  inline size_t getNumClusters() const noexcept
  {
    return clusters.size();
  }

  /// @brief clear resets the state of the MTDBSCAN clusters to be empty.  It allows setting a new window
  /// and calling buildCluster using the data retained in the kdTree to get a new cluster.
  void clear()
  {
    for (auto cl : clusters) {
      delete cl;
    }
    clusters.clear();
    if (kdTree != nullptr)
      kdTree->reset();
    IDCounter = 0;
  }


  /// @brief checkCluster does some basic check on the clusters and print some statistics
  /// @param numValues total number of values the clusters should contain
  /// @param tag if a tag is provided statistics are only checked for clusters with that tag.
  /// @return true if no errors are found.
  bool checkClusters(const size_t numValues, const std::string *tag = nullptr) const noexcept
  {
    size_t max = std::numeric_limits<size_t>::min();
    size_t min = std::numeric_limits<size_t>::max();
    size_t avgSize = 0;
    size_t count = 0;
    bool rb = true;
    for (size_t i = 0; i < clusters.size(); i++)
    {
      Cluster *c = clusters.at(i);
      if (tag == nullptr || (c->tag != nullptr && c->tag == tag))
      //if (tag == nullptr)
      {
        size_t t = c->values->size();
        if (t == 0)
        {
          std::cout << "Cluster " << i << " has 0 entries" << std::endl;
        }
        if (t > max)
          max = t;
        if (t < min)
          min = t;
        avgSize += t;
        count++;
      }
    }
    if (numValues != -1 && avgSize != numValues)
    {
      std::cout << "Number of locations in all clusters not equal input locations." << std::endl;
      rb = false;
    }
    avgSize = (count > 0) ? avgSize / count : 0;
    min = (count > 0) ? min : 0;
    max = (count > 0) ? max : 0;
    std::string t_tag = (tag == nullptr ? "" : *tag + " ");
    std::cout << "Cluster " << t_tag << "Count = " << count << " Max = " << max << "  Min = " << min << " Average = " << avgSize << std::endl;
    return rb;
  }
};

#endif // MTDBSCAN_HPP
