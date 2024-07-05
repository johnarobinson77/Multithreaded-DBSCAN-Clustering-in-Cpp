# A Multithreaded DBSCAN Clustering Implementation Using a K-d Tree
\[ Special Note:  a previous version of this code had a performance pothole when the cluster size was large.  That pothole has been substantially mitigated in this version. \]

The code included here implements the popular DBSCAN clustering algorithm [1] in C++. it is a multithreaded version of [DBSCAN-in-Cpp](https://github.com/johnarobinson77/DBSCAN-in-Cpp). Like the single-threaded version, this implementation uses a k-d tree to accelerate the cluster-building process. The k-d tree code is also presented here. Using a k-d tree improves the clustering performance in three ways.

First, the k-d tree code presented here is highly optimized to build and search the tree using multithreading, significantly reducing the time needed for window searches.

Second, it supports marking the nodes in the k-d tree so that whole branches can be marked as already taken. This eliminates redundant visits to a node while building the clusters.

Finally, it facilitates the multithreading of the clustering algorithm itself. In the single-threaded version, the process of building the k-d tree is multithreaded, but the process of creating the clusters is single-threaded. In this multithreaded version, the clustering algorithm is also multithreaded, and as described later, the k-d tree is used to overcome the issues arising from multithreading.

Note that because of the overhead required to implement multithreading, the multithreaded version using a single thread is not as fast as the single-threaded version. If multithreading is not required, use [DBSCAN-in-Cpp](https://github.com/johnarobinson77/DBSCAN-in-Cpp).

## Usage
This implementation assumes a data set has a position in some n-space, geometric or otherwise. The position is represented in a vector of K-type numbers, and an object called a value here can be an object or a pointer to an object. Each position and value are added as a pair to the MTDBSCAN object. Also, a cluster window needs to be provided. After running the clustering procedure, a vector of clusters will be available, where each cluster will be a list of the values presented to the MTDBSCAN object earlier.

In the example below, that data is stored as a class in the Locations class. The usage steps are as follows:

1. Create an instance of the MTDBSCAN class,
2. Loop through the positions and values, presenting them in pairs to the MTDBSCAN object
3. Call the build method.
4. Call the checkCluster function to get the statistics of the clustering if desired.
5. Access the clusters vector as needed.

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
  //Get the search range about the cluster distance window
  std::vector<dkey_t> window(numDimensions);
  for (int i = 0; i < numDimensions; i++) {
    window[i] = searchRad;
  }

  std::cout << "Running DBSCAN clustering algorithm..." << std::endl;
  mtdbscan->setWindow(window);
  startTime = std::chrono::steady_clock::now();;
  mtdbscan->build();
  endTime = std::chrono::steady_clock::now();

  double totalMtdbscanTime = std::chrono::duration<double>(endTime - startTime).count();

  std::cout << "Total cluster time = " << std::fixed << std::setprecision(4) << totalMtdbscanTime << " seconds" << std::endl << std::endl;

  if (!mtdbscan->checkClusters(numPoints)) {
    exit(1);
  }

  std::cout << "Checking to see that all values ended up in some cluster" << std::endl;

  avlSet<dval_t> valMatch;
  bool passed = true;
  for (size_t i = 0; i < mtdbscan->getNumClusters(); i++) {
    for (dval_t val : *mtdbscan->getClusterValueList(i)) {
      passed &= valMatch.insert(val);
      passed &= val < static_cast<dval_t>(numPoints) && val >= static_cast<dval_t>(0);
    }
  }
  if (passed && valMatch.size() == numPoints) {
    std::cout << "Passed" << std::endl;
  }
  else {
    std::cout << "Failed to find the same set values in the clusters as generated." << std::endl;
  }
  valMatch.clear();

  mtdbscan->sortClustersBySize();
```

Note that this interface is the same as the single-threaded version in [DBSCAN-in-Cpp](https://github.com/johnarobinson77/DBSCAN-in-Cpp) and therefore is a drop-in replacement.

## Description of Source Files
### Implementation Files
The following files are needed to implement multithreaded DBSCAN in your code.

**MTDBSCAN.hpp** holds the MTDBSCAN class definition. The DBSCAN algorithm is performed in the build function there.

**KdTreePDB.hpp** contains the KdTree class and the KdNode class. The KdTreePDB class is an API wrapper around the KdNode class. The KdNode class implements the building of the k-d tree using the Knlogn algorithm. The algorithm description for building a k-d tree and original code can be found at [https://github.com/chezruss/kd-tree](https://github.com/chezruss/kd-tree). Much of the code in the KdNode class for building the k-d tree is taken directly from that code base. The search code in the KdNode class was rewritten to accelerate the particular problem of multithreaded DBSCAN,

**ParallelSort.hpp, and ParallelFor.hpp** contain code for sorting in parallel. It is used in the kdTree build code as well as the cluster size sort function.

**SimpleContainers.h**  provides set, list, and memory allocator classes that are faster and take less memory than the STL versions for this particular use case.

**avlSet.h** avlSet.h[3] is a slightly modified version of avlTree.cpp that can be found at [https://github.com/chezruss/avl-tree](https://github.com/chezruss/avl-tree).

### Test files
The following files are provided to test the multithreaded DBSCAN code and are not needed to implement this DBSCAN package implementation in your code.

**DBSCANtest.cpp** provides an example of usage and, a test case for the DBSCAN implementation. The default test case is 1600 artificially clustered data with 4000 points per cluster. That data is presented to MTDBSCAN. The results should, therefore, end up with a Clusters list that is 1600 in size, and clusters will have an average of 4000 values. The parameters of the clustering can be overridden using command line arguments. Start with -h.  This file is not necessary other than to test the MTDBSCAN class. Compile with -std=c++17.

**ParseArgs.h** is my attempt at providing a general command line parser.  It is used in DBSCANtest.  This file is not necessary other than to test the MTDBSCAN class.

## MTDBSCAN Class Methods

```c++
 // base constructor
  MTDBSCAN< K, V, N>()
  // where: K is the type of the components of the point associated with the value
  //        V is the type or class of the value
  //        N is the number of components of the point

  // constructor with window
  MTDBSCAN< K, V, N >(std::vector<K>& window)

  // addPointWithValuefuction adds a point or key and an associated value to the dbscan object
  size_t addPointWithValue(std::vector<K> point, V value)

  // setWindow sets the clustering window radius oe=r 1/2 the search kernal size.  It must have the same number of values as N.  
  void setWindow(std::vector<K>& window)

  // build builds the clusters.  All points and values must be added first.
  // It will return false if called a second time after dbscan was constructed.
  bool build()

  //sortClusterBySize sorts the cluster by the number of values in the cluster from largest to smallest.  To implement
  // noise cutoff, call this function, and don't use the values at the end of the vector beyond the cutoff.
  void sortClustersBySize()

  // cluster setters and getters
  // getClusterSize returns the number of values in cluster clusterIdx
  size_t getClusterSize(size_t clusterIdx)

  // getClusterValueList returns a pointer to the list of values in cluster clusterIdx
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
A cluster in DBSCAN is a set of points that are all within a user-defined distance from each other. The cluster size can be much larger than that distance, but each point must be within that distance from some other point in the cluster. In this implementation, that distance is set by the window() function.

The process of building a cluster is relatively straightforward.

1. Pick an arbitrary point from the data set and create a new cluster starting with that point.
2. For each point in the cluster, search for all the data set points within the user-defined distance and add them to the cluster. Include the new points in the cluster as you go.
3. Once the search around all the points in the cluster is complete, and no new points have been added to the cluster, that cluster is complete.
4. Request the next arbitrary point from the data set not already added to a cluster. If a point is returned, repeat from step 2. If no points are left, the clustering process is complete.

In the code provided here, the data set is stored in a k-d tree where each point in the original data set is represented by a node in the tree. The picking and searching operations described above are performed on the k-d tree. When a point is added to a cluster by the picking or searching operation, its node in the tree is marked as taken, and that point is ignored on all subsequent operations.

### Multi-threading the Clustering Process
At first, it might seem that you could start multiple threads, each building its own list of clusters. When all threads have returned, simply concentrate the per-thread lists into a single list of clusters. However, consider the possibility that two threads may pick different points that eventually want to be placed in the same cluster. Since the picking process chooses arbitrary untaken points, nothing prevents this from happening. 

The result of the multithreaded process is that some clusters have been split into two or more sub-clusters. Worse, which clusters have been split will change from run to run, depending on the timing of how the threads interact. This happens even when rerunning the process on the same data set.

### Finding Associated Subclusters
A solution to this is to identify and merge the sub-clusters after the multi-threaded clustering process is complete. First, the code needs to identify which sub-clusters need to be merged. This is done as follows.

In the picking and search routines, when a point is found that needs to be added to the current cluster, the node in the k-d tree is marked as having been taken by moving the values from that node to the cluster and storing a pointer to the current cluster that the values have been moved to in the node. When a subsequent search finds the same node that wants to be added to its cluster, one of two things happens.

1. A movedTo pointer that is the same as the current cluster pointer indicates that the point is already in the current cluster. It is just ignored.
2. A movedTo pointer that is different from the current cluster pointer indicates that the current cluster overlaps with the cluster to which the values have been moved. The movedTo and current clusters will need to be merged in the future. To identify this, the movedTo pointer is added to a set in the current cluster called the overlap set.

### Merging Subclusters
After the threads have built their list of clusters, those lists are scanned to find the subclusters that need to be merged, indicated by a non-null overlap set. The set of subclusters that need to be merged can be simple or quite complex. It is helpful to think of the set of subclusters as a graph where the subcluster is a node in that graph, and the links are the pointers in the overlap set. A subcluster may not have a direct link to all the subclusters in the cluster set. Several links may need to be followed to find them all. In addition, this overlap graph will have loops where following the overlap pointers from subcluster to subcluster ends up back at the original cluster.

The following post-clustering process is performed to find and merge these subclusters.

1. Scan the clusters in the per-thread cluster lists for a cluster with a non-null overlap set and non-zero length. We'll refer to this cluster as the initiating cluster. Note: a zero-length subcluster indicates that its values have already been moved to an initiating cluster.
2. Merge all the subclusters in the overlap graph into the initiating cluster by calling the overlapMerge() function with a pointer to the initiating cluster and a pointer to the overlap set. The overlapMerge() function iterates over the overlap set.
  1. If an overlap cluster is zero length or is the same as the initiating cluster, it is skipped over because it indicates a loop in the graph.
  2. If not, its values are added to the initiating cluster, and its length is zeroed. Then overlapMerge is called with its overlap set and a pointer to the initiating cluster.
3. At any recursion level, once all the clusters in the overlap set have been addressed, the code returns from that level.  Eventually, all of the nodes in that overlap set except the initiating cluster will have zero length. 
4. When all subclusters have been merged into final clusters, scan the cluster list again and remove all zero-length clusters. Also, delete any non-null overlap sets, as they are no longer needed. 

## Performance
### Limiting Subcluster Overlap
While overlapping subclusters are functionally required, they hinder performance and should be avoided. A way to do that is for each thread to work on clusters that are geometrically as far apart as possible. The points in the most widely separated clusters will be on opposite sides of the k-d tree.

The picker function that provides a seed point for a new cluster has a selectionBias parameter used to select a path through the k-d tree. The thread number is used to create a selectionBias that chooses a seed point that is as far away as possible from the other threads.

### Limiting k-d tree Search Depth
The nature of searching the k-d tree for points to add to the cluster is such that nodes in the tree will be visited multiple times unless something is done to mitigate those multiple visits.  If it is known at a node in the k-d tree that this node and all the nodes below on both child branches have already been added to this cluster, there is no need to search those branches again.

To achieve this, each cluster generated has a unique ID assigned at construction time.  Each k-d tree node contains a set called the AllVisited set.  During the recursive search, when a node is found to include the point in that node, the cluster ID is added to the node's AllVisited set but only if the two child nodes have that ID in their set.  On return to the parent node, the return code indicates whether it is an AllVisited node for this cluster.

On the downward traversal of subsequent searches, when the thread finds the cluster-ID it is working on in the AllVisited set, it returns without going down that tree branch.  This goes a long way to keep the k-d tree cluster search from being N^2 complexity. 

### Memory Allocation
For the purposes of building the clusters, the k-d tree is an intermediate data structure that is destroyed at the end of the cluster-building process.  The k-d tree will contain a k node for each point added to the MTDBSCAN.  So for a large data set, the temporary memory taken by the k-d tree can be quite large, but worse it is allocated in small chunks, which are the node's size.

During performance tuning, I found that the time to delete all the nodes in the k-d tree was substantial.  I believe that this was due to the memory allocator having to coalesce all the small per-node allocations.  On some systems, such as in MSVC, that process took place at the time the node destructor was called.  On other systems, such as g++ on Ubuntu, the coalescing seemed to happen in two parts, first when the node destructor was called and then again when a new large container, such as a std::vector, was allocated.

To avoid small memory allocation fragmentation, I wrote an allocator for the k-d tree nodes, referred to as SimpleAlloc in the code.  SimpleAlloc uses malloc to allocate large blocks of heap memory.  When a new k-d tree node is constructed, it gets a pointer from SimpleAlloc, which simply puts the next allocation after the previous one.  When SimpleAlloc finds insufficient memory in the current block, it mallocs a new large block and puts the allocation in that block.  

SimpleAlloc keeps a std::list of all the blocks it allocates.  When it is time to deallocate the memory for the k-d tree, free is called on each of the large blocks in the list, speeding up the process to the point where it is no longer a significant part of the processing time.

### Stepping Away from STL Containers
As mentioned in the Limiting k-d tree Search Depth section, each node has an allVisited set. Since this set is in the inner loop of the cluster search of the k-d tree, it must be efficient for this use case. I found that writing my own set that was optimized for small set sizes and the most basic operations was superior in performance and size to using std::set or std::unordered_set. It is referred to as SimpleSet in the code. It also needed to be compatible with SimpleAlloc, as described above.

Likewise, the set that keeps track of overlapped clusters uses avlSet[3].

## Performance Results
Two major processes contribute most to the execution time: building the k-d tree and searching the k-d tree to build the clusters.  Building the k-d tree is always done with the next power of 2 threads less than the requested number of threads, while the cluster search uses all the requested threads.  For instance, if 8 threads are requested, 8 will be used for both the k-d tree build and the cluster search.  If 6 threads are requested, 4 will be used for the k-d tree build and 6 for the cluster search.  That is why you will observe a bump in performance at 2,4,8,16 and 32 threads.

Except where noted, the data set presented to MTDBSCAN is 6.4 million points in 1600 clusters of 4000 points each.  The window size is set so the resulting clusters match the input data.

### Amdahl Model
The number of threads vs. time Figures shows the processing time and the result for fitting it to an Amdahl model.  The circles indicate the measured data, and the dashed line indicates the best fit to a 3-parameter Amdahl model.  The Amdahl model is described by the following equation:
 > t = t<sub>S</sub> + t<sub>1</sub>/q + t<sub>C</sub>(q-1)       
   where:  
     t is the total execution time  
     t<sub>S</sub> is the time to execute the non-parallelizable code  
     t<sub>1</sub> is the time to execute parallelizable code  
     t<sub>C</sub> is the time loss due to contention between threads.  
     q is the number of threads

A more complete description of the Amdahl model can be found in in section 6.10 in Reference [2] below. 

### Core i7 Results
Figure 1 shows the processing time on a 24-core Intel i7 processor running on Ubuntu under WSL. The processor actually has 8 cores with 2 hyper threads each plus 8 efficiency cores. The code was compiled with g++-11.  

![Intel i7 Multithreaded Performance](https://github.com/johnarobinson77/Multithreaded-DBSCAN-Clustering-in-Cpp/blob/main/hp_gpp_total.png)
**<p align="center">Figure 1. Intel i7 Performance</p>**

### Graviton Results
The remaining figures show the performance as measured on a 48-core Graviton processor on AWS.  It is a better example of using this solution on an enterprise-class processor.

Figure 2 shows the performance on the same 1600 cluster 4000 points per cluster data out to 48 threads.  Note the small value of t~S~ indicating that the threads are not interfering with each other for this data set.

![Figure 2. Graviton Multithreaded Performance](https://github.com/johnarobinson77/Multithreaded-DBSCAN-Clustering-in-Cpp/blob/main/gvt_gpp_total_v2.png)
 **<p align="center">Figure 2. Graviton Performance</p>**
### Extreme Cases    
The above performance figures show the performances for relatively large numbers of clusters.  However, the performance tends to suffer when there are large numbers of clusters with few values per cluster or small clusters with a large number of points per cluster.  Chart 3 shows the performance of different cluster sizes.  In each case, the number of points in the input data is the same, 6.4 million, but arranged into different numbers of clusters as indicated on the horizontal axis.

![Figure 3. Performance vs. Cluster Count](https://github.com/johnarobinson77/Multithreaded-DBSCAN-Clustering-in-Cpp/blob/main/ClusterCntVsTime.png)
**<p align="center">Figure 3. Performance vs Cluster Count</p>**

Figure 4 shows the performance vs number of threads for a 6.4 million point data set grouped into only 16 clusters.  When the number of threads exceeds the number of clusters, there is more overlapped searching of the k-d tree, and therefore more overlapping clusters formed.  Note that t~S~ is much larger, and the processing time increases with increasing thread count.

![Figure 4. 16 Cluster Performance](https://github.com/johnarobinson77/Multithreaded-DBSCAN-Clustering-in-Cpp/blob/main/gvt_gpp_total_X_v2.png)
**<p align="center">Figure 4. 16 Cluster Performance</p>**

## Thanks
Thanks to Russel Brown for providing the figures with the Amdahl Model curves. 

## References
[1] Ester, Martin; Kriegel, Hans-Peter; Sander, Jörg; Xu, Xiaowei (1996). Simoudis, Evangelos; Han, Jiawei; Fayyad, Usama M. (eds.). A density-based algorithm for discovering clusters in large spatial databases with noise. Proceedings of the Second International Conference on Knowledge Discovery and Data Mining (KDD-96). AAAI Press. pp. 226–231. CiteSeerX 10.1.1.121.9220. ISBN 1-57735-004-9.

[2] Russell A. Brown, Building a Balanced k-d Tree in O(kn log n) Time, [arXiv:1410.5420v45](https://arxiv.org/abs/1410.5420v45) [cs.DS] .

[3] Code for the avlSet was provided by Russell A. Brown from his repository at https://github.com/chezruss/avl-tree

