# A Multithreaded DBSCAN Clustering Implementation Using a K-d Tree

The code included here implements the popular DBSCAN clustering algorithm [1] in C++. it is a multithreaded version of [DBSCAN-in-Cpp](https://github.com/johnarobinson77/DBSCAN-in-Cpp). Like the single-threaded version, this implementation uses a k-d tree to accelerate the cluster-building process. The k-d tree code is also presented here. Using a k-d tree improves the clustering performance in three ways.

First, the k-d tree code presented here is highly optimized to build and search the tree using multithreading, significantly reducing the time needed for window searches.

Second, it supports marking the nodes in the k-d tree in such a way that whole branches of the tree can be marked as already taken. This eliminates redundant visiting a node while building the clusters.

Finally, it is used to facilitate multithreading the clustering algorithm itself. In the single-threaded version, the process of building the k-d tree is multithread, but the process of creating the clusters is single-threaded. In this multithreaded version, the clustering algorithm is also multithreaded, and as described later, the k-d tree is used to overcome the issues arising from multi-threading.

Note that because of the overhead required to implement multithreading, the multithreaded version using a single thread is not as fast as the single-threaded version. If multithreading is not required, use [DBSCAN-in-Cpp](https://github.com/johnarobinson77/DBSCAN-in-Cpp).

## Usage

This implementation assumes a set of data has a position in some n-space, geometric or otherwise. The position is represented in a vector of K-type numbers and an object called a value here, which can be an object or pointer to an object. Each position and value is presented as a pair to the MTDBSCAN object. Also, a cluster window needs to be provided. After running the clustering procedure, a list of clusters will be available where each cluster will be a list of the values presented to the MTDBSCAN object earlier.

In the example below, that data is stored as a class in the Locations class. The usage steps are as follows:

1. Create an instance of the MTDBSCAN class,
2. Loop through the positions and values, presenting them in pairs to the MTDBSCAN object
3. Call the build method.
4. Call the checkCluster function to get the statistics of the clustering if desired..
5. Access the clusters list as needed.

This source code example below can be found in the main function in DBSCANtest.cpp

```c++
#include "MTDBSCAN.hpp"
	:
	:
 std::cout << "Adding data to DBSCAN..." << std::endl;
  auto mtdbscan = new MTDBSCAN< dkey_t, dval_t, numDimensions>;
  // Add each pair of points and values to the dbscan object
  for (size_t i = 0; i < numPoints; ++i) {
    mtdbscan->addPointWithValue(coordinates[i], i);
  }
  // get the search range about cluster distance window
  std::vector<dkey_t> window(numDimensions);
  for (int i = 0; i < numDimensions; i++) {
    window[i] = searchRad;
  }

  std::cout << "Running DBSCAN clustering algorithm..." << std::endl;
  mtdbscan->setWindow(window);
  timespec startTime = getTime();
  mtdbscan->build();
  timespec endTime = getTime();

  double totalTime = (endTime.tv_sec - startTime.tv_sec) +
    1.0e-9 * ((double)(endTime.tv_nsec - startTime.tv_nsec));
  std::cout << "Total cluster time = " << std::fixed << std::setprecision(2) << totalTime << " seconds" << std::endl << std::endl;

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

  mtdbscan->sortClustersBySize();
```

Note that this interface is exactly the same as the single-threaded version in [DBSCAN-in-Cpp](https://github.com/johnarobinson77/DBSCAN-in-Cpp) and therefore is a drop-in replacement.

## Description of Source Files

**MTDBSCAN.hpp** holds the MTDBSCAN class definition. The DBSCAN algorithm is performed in the build function there.

**KdTreePDB.hpp** contains the KdTree class and the KdNode class. The KdTreePDB class is an API wrapper around the KdNode class. The KdNode class implements the building of the k-d tree using the Knlogn algorithm. The algorithm description for building a k-d tree and original code can be found at [https://github.com/chezruss/kd-tree](https://github.com/chezruss/kd-tree). Much of the code in the KdNode class for building the k-d tree is taken directly from that code base. The search code in the KdNode class was rewritten to accelerate the particular problem of multithreaded DBSCAN,

**ParallelSort.hpp, and ParallelFor.hpp** contain code for sorting in parallel. It is used in the kdTree build code.

**DBSCANtest.cpp** provides a test and an example and test case for the DBSCAN implementation The default test case is 1600 artificially clustered data with 4000 points per cluster. That data is presented to MTDBSCAN. The results should, therefore, end up with a Clusters list that is 1600 in size, and clusters will have an average of 4000 values. The parameters of the clustering can be overridden using command line arguments. Start with -h.

## MTDBSCAN Class Methods

```c++
 // base constructor
  DBSCAN()

  // constructor with window
  DBSCAN(std::vector<K>& window)

  // addPointWithValuefuction adds a point or key and an associated value to the dbscan object
  size_t addPointWithValue(std::vector<K> point, V value)

  // setWindow sets the clustering window size
  void setWindow(std::vector<K>& window)

  // build builds the clusters.  All points and values must be added first.
  // It will return false if called a second time after dbscan was constructed.
  bool build()

  //sortClusterBySize sorts the cluster by the number of values in the cluster from largest to smallest
  void sortClustersBySize()

  // cluster setters and getters
  // getClusterSize returns the number of values in cluster clusterIdx
  size_t getClusterSize(size_t clusterIdx)

  // getClusterValueList returns the list of values in cluster clusterIdx
  std::list<V>* getClusterValueList(size_t clusterIdx) 

  // getClusterMaxCorner returns a vector of the maximum corner of the
  // bounding hypercube of cluster clusterIdx
  std::vector<K>* getClusterMaxCorner(size_t clusterIdx)

  // getClusterMinCorner returns a vector of the minimum corner of the
  // bounding hypercube of cluster clusterIdx
  std::vector<K>* getClusterMinCorner(size_t clusterIdx) 

  // getClusterTag returns the tag string set by setClusterTag
  std::string* getClusterTag(size_t clusterIdx)

  // cluster tag setter
  bool setClusterTag(size_t clusterIdx, std::string& tagIn)

  // getNumClusters returns the number of clusters
  size_t getNumClusters()

  // checCluster does some basic checks on the clusters and prints some statistics
  bool checkClusters(const size_t numLocations, const std::string* tag = nullptr)  

 // Sets the number of threads used for building the clusters.  The default is to use the hardware concurrency value for the current processor.
  void setNumThreads(int64_t thd = -1)
```

## The Multithreaded DBSCAN Algorithm

### A Brief Explanation of DBSCAN

A cluster in DBSCAN is a set of points that are all within a user-defined distance from each other. The size of the cluster can be much larger than that distance. But each point in that cluster must be within that distance from some other point in the cluster. In this implementation, that distance is set by the window() function.

The process of building a cluster is relatively straightforward.

1. Pick an arbitrary point from the data set and create a new cluster starting with that point.
2. For each point in the cluster, search for all the data set points within the user-defined distance and add them to the cluster. Include the new points in the cluster as you go.
3. Once the search around all the points in the cluster is complete, and no new points have been added to the cluster, that cluster is complete.
4. Request the next arbitrary point from the data set not already added to a cluster. If a point is returned, repeat from step 2. If no points are left, the clustering process is complete.

In the code provided here, the data set is stored in a k-d tree where each point in the original data set is represented by a node in the tree. The picking and searching operations described above are performed on the k-d tree. When a point is added to a cluster by the picking or searching operation, its node in the tree is marked as taken, and that point is ignored on all subsequent operations.

### Issues with Multi-threading the Clustering Process

At first, it might seem that you could start multiple threads, each building its own list of clusters. When all threads have returned, simply concentrate the per-thread lists into a single list of clusters. But consider the possibility that two threads may pick different points that eventually want to be placed in the same cluster. Since the picking process chooses arbitrary untaken points, there is nothing that guarantees this won't happen. The probability that it will happen within a multi-threaded cluster build process is near 100%.

The result of the multithreaded process is that some clusters have been split into two or more sub-clusters. Worse, which clusters have been split will change from run to run, depending on the timing of how the threads interact. This happens even when rerunning the process on the same data set.

### Finding Associated Subclusters

A solution to this is to merge the sub-clusters after the multi-threaded clustering process is complete. First, the code needs to identify which sub-clusters need to be merged. This is done as follows.

In the picking and search routines, when a point is found that needs to be added to the current cluster, the node in the k-d tree is marked as having been taken by moving the values from that node to the cluster and storing a pointer to the current cluster that the values have been moved to in the node. When a subsequent search finds the same node that wants to be added to its cluster, one of two things happens.

1. A movedTo pointer that is the same as the current cluster pointer indicates that the point is already in the current cluster. It is just ignored.
2. A movedTo pointer that is different from the current cluster pointer indicates that the current cluster overlaps with the cluster to which the values have been moved. The movedTo and current clusters will need to be merged in the future. To identify this, the movedTo pointer is added in an std::set in the current cluster called the overlap set.

### Merging Subclusters

After the threads have built their individual list of clusters, those lists are scanned to find the subclusters that need to be merged, indicated by a non-null overlap std::set. The set of subclusters that need to be merged can be simple or quite complex. It is helpful to think of the set of subclusters as a graph where the subcluster is a node in that graph, and the links are the pointers in the overlap std::set. A subcluster may not have a direct link to all the subclusters in the cluster set. Several links may need to be followed to find them all. In addition, this overlap graph will have loops where following the overlap pointers from subcluster to subcluster ends up back at the original cluster.

To find and merge these subclusters, the following post-clustering process is performed.

1. Scan the clusters in the per-thread cluster lists for a cluster with a non-null overlap std::set and non-zero length. We'll refer to this cluster as the initiating cluster. Note: a zero-length subcluster indicates that its values have already been moved to an initiating cluster.
2. Merge all the subclusters in the overlap graph into the initiating cluster by calling the overlapMerge() function with a pointer to the initiating cluster and a pointer to the overlap std::set. The overlapMerge() function iterates over the overlap set.
  1. If an overlap cluster is zero length or is the same as the initiating cluster, it is skipped over because it indicates a loop in the graph.
  2. If not, its values are added to the initiating cluster, and its length is zeroed. Then overlapMerge is called with its overlap set and a pointer to the initiating cluster.
3. When all subclusters have been merged into final clusters, scan the cluster list again and remove all zero-length clusters. Also, delete any non-null overlap std::stes as they are no longer needed. This "cleaning" function is done in this implementation as part of a move from the per-thread lists to the final consolidated cluster list.

### Other Performance Considerations

While overlapping subclusters are functionally resolved, they hinder performance and should be avoided as much as possible. A way to do that is for each thread to work on clusters that are geometrically as far apart as possible. The points in the most widely separated clusters will be on opposite sides of the k-d tree.

The picker function used to provide a seed point for a new cluster has a selectionBias parameter used to select a path through the k-d tree. The thread number is used to create a selectionBias that chooses a seed point that is as far away as possible from the other thread.

Another optimization is to mark each child branch in a node in the k-d tree with an index of a cluster that has used all of the points below. If a search for a cluster with the same index is performed, that branch can be skipped. However, if a search for a cluster with a different index is in progress, that branch still needs to be searched for only overlaps.

### Performance Pothole
There is one performance pothole with the current code that can cause extreme performance degradation.  This happens when the search region is large relative to the area containing the data set, such that the results will end up with just a few clusters.  The requirement that a thread must search for a point in the cluster it's forming that has been taken by other clusters can result in very large search times in the k-d tree.  This performance degradation is only present when running multithreaded.  The current code may slow down as much as 4x relative to single threading.  

Doing dbscan clustering with such a large search region is generally not very useful.  However, depending on the algorithm's use, it may not be avoidable.  There are several ways that this performance pothole is being mitigated, and more are being explored. 

## Results

The chart below shows how processing time diminishes with increasing thread.  This performance data was taken on a 64-core Graviton processor on AWS.  The performance flattens out above 32 threads, likely due to memory bandwidth and increased overlapping sub-clusters.  The data set is 1600 artificially created groups of 4000 3-dimensional points each.  The search region is chosen to identify those groups as clusters.  These are the default settings in DBSCANtest.cpp.

![Multithreaded DBSCAN Performance](https://github.com/johnarobinson77/Multithreaded-DBSCAN-Clustering-in-Cpp/blob/main/MTDBSCAN.PNG)

## References

[1] Ester, Martin; Kriegel, Hans-Peter; Sander, Jörg; Xu, Xiaowei (1996). Simoudis, Evangelos; Han, Jiawei; Fayyad, Usama M. (eds.). A density-based algorithm for discovering clusters in large spatial databases with noise. Proceedings of the Second International Conference on Knowledge Discovery and Data Mining (KDD-96). AAAI Press. pp. 226–231. CiteSeerX 10.1.1.121.9220. ISBN 1-57735-004-9.

[2] Russell A. Brown, Building a Balanced k-d tree in O(kn log n) Time, Journal of Computer Graphics Techniques (JCGT), vol. 4, no. 1, 50-68, 2015.
