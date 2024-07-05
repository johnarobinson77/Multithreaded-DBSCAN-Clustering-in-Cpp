/*
 * Copyright (c) 2023 John Robinson.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

 //
 //  KdTreePDB.h
 //
 //  Created by John Robinson on 11/18/23.
 // 
 //  Note that most of the code contained here is original code written
 //  by the author.  However, the code specific to building the kdTree,
 //  specifically the createKdTree function in the kKdNode class and the
 //  functions it calls were taken from code written by Russel Brown and  
 //  can found at https://github.com/chezruss/kd-tree.
 //



#ifndef KDTREE_HPP
#define KDTREE_HPP

#include <vector>
#include <list>
#include <set>
#include <exception>
#include <future>
#include "parallelSort.hpp"
#include "avlSet.h"
#include "SimpleContainers.h"

std::atomic<size_t> nodesVisited = 0;

/*
 * This type is the signed equivalent of size_t and might be equivalent to intmax_t
 */
typedef int64_t signed_size_t;

enum SearchRet { // The following are codes that the search routines will return indicating the status of then child node visited.
  Undefined,                // indicates the values, ltChile and gtChild pointer are all null.
  NoResult,                 // indicates no result has been returned.
  AllTakenByID,             // indicates all nodes below have been taken by current cluster
  AllTaken,                 // indicates all nodes below have been taken by different clusters
};

// this is an ID that is used to indicate all taken when the nodes below do not match all IDs.
const uint32_t allTakenID = 0xFFFFFFFF;


// the KdTreePDB class provide a user interface wrapper around the kdNodePDB objects that 
// form the actual k-d tree.
template<typename K, typename V, typename M, size_t N>
class KdTreePDB {

  // this is ther return type
  typedef std::pair<K*, std::list<V>*> retPair_t;


  class KdNodePDB {

    typedef std::pair<K*, std::list<V>*> retPair_t;

  private:
    K tuple[N];
    KdNodePDB* ltChild;
    KdNodePDB* gtChild;
    std::list<V, STLAdaptor<V>>* values;  // list of values associated with the point in tuple
    std::mutex thisMutex;       // used to lock this node when it is being taken by a cluster or updating the allTakenSet
    M* movedTo = nullptr;       // pointer to the container that is holding the values that were here.
    SimpleSet allTakenSet;      // set of IDs of the clusters that have taken values from this node and all nodes below.

  public:
    // need to make sure allocator pointers are provided when construction occures so delete the default
    KdNodePDB() = delete;

    // constructor for entering the point only.
    KdNodePDB(V const value, STLAdaptor<V>* adaptorV, SimpleAlloc* alloc) { // Pass non-primitive types as 'V const&'
      for (size_t i = 0; i < N; i++) this->tuple[i] = (K)0;
      ltChild = gtChild = nullptr; // redundant
      // create the lest of values
      void* tp = alloc->allocate(sizeof(std::list<V, STLAdaptor<V>>), alignof(std::list<V, STLAdaptor<V>>));
      values = new(tp)  std::list<V, STLAdaptor<V>>(*adaptorV);
      values->push_back(value);
      // create the list for SimpleSet used for the allTakenSet
      allTakenSet.setup(alloc);
    }

    // constructor for entering the point and value
    KdNodePDB(std::vector<K>& tuple, V const value, STLAdaptor<V>* adaptorV, SimpleAlloc* alloc)
    { // Pass non-primitive types as 'V const&'
      for (size_t i = 0; i < N; i++) this->tuple[i] = tuple[i];
      ltChild = gtChild = nullptr; // redundant
      // create the lest of values
      void* tp = alloc->allocate(sizeof(std::list<V, STLAdaptor<V>>), 
        alignof(std::list<V, STLAdaptor<V>>));
      values = new(tp)  std::list<V, STLAdaptor<V>>(*adaptorV);
      values->push_back(value);
      // create the list for SimpleSet used for the allTakenSet
      allTakenSet.setup(alloc);
    }

  public:
    ~KdNodePDB() {
      std::cout << "Should not have been called.\n";
    }

  public:
    K const* getTuple() {
      return this->tuple;
    }

    /*
      * The superKeyCompare function compares two K arrays in all k dimensions,
      * and uses the sorting or partition coordinate as the most significant dimension.
      *
      * Calling parameters:
      *
      * a - a K*
      * b - a K*
      * p - the most significant dimension
      *
      * returns a K result of comparing two K arrays
      */
  private:
    inline
      static K superKeyCompare(K const* a, K const* b, size_t p) {
      // Typically, this first calculation of diff will be non-zero and bypass the 'for' loop.
      K diff = a[p] - b[p];
      for (size_t i = 1; diff == 0 && i < N; i++) {
        size_t r = i + p;
        // A fast alternative to the modulus operator for (i + p) < 2 * N.
        r = (r < N) ? r : r - N;
        diff = a[r] - b[r];
      }
      return diff;
    }

    /*
      * The removeDuplicates function checks the validity of the merge sort and
      * removes duplicates from the kdNodes array.
      *
      * Calling parameters:
      *
      * kdNodes - a KdNode** array that has been sorted via merge sort according to (x,y,z,w...) tuples
      * i - the leading dimension for the super key
      *
      * returns the end index of the reference array following removal of duplicate elements
      */
  private:
    inline
      static size_t removeDuplicates(KdNodePDB** kdNodes, size_t i, size_t size) {
      size_t end = 0;
      for (size_t j = 1; j < size; ++j) {
        K compare = superKeyCompare(kdNodes[j]->tuple, kdNodes[end]->tuple, i);
        if (compare < 0) {
          std::cout << "merge sort failure: superKeyCompare(kdNodes[" << j << "], kdNodes["
            << end << "], " << i << ") = " << compare << std::endl;
          exit(1);
        }
        else if (compare > 0) {
          // Keep the jth element of the kdNodes array.
          kdNodes[++end] = kdNodes[j];
        }
        else {
          // append vales to the end of the kept kdNode and delete the now unused node
          kdNodes[end]->values->splice(kdNodes[end]->values->end(), *(kdNodes[j]->values));
          //delete kdNodes[j];
        }
      }
      return end;
    }

    /*
      * The buildKdTree function builds a k-d tree by recursively partitioning
      * the reference arrays and adding KdNodes to the tree.  These arrays
      * are permuted cyclically for successive levels of the tree in
      * order that sorting occur in the order x, y, z, w...
      *
      * Calling parameters:
      *
      * reference - a KdNode*** array to recursively sort via its (x, y, z, w...) tuples array
      * temporary - a KdNode*** temporary array from which to copy sorted results;
      * start - start element of the reference array
      * end - end element of the reference array
      * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
      * depth - the depth in the tree
      *
      * returns: a KdNode pointer to the root of the k-d tree
      */
  private:
    static KdNodePDB* buildKdTree(KdNodePDB*** references,
      std::vector< std::vector<size_t> > const& permutation,
      size_t start, size_t end,
      signed_size_t maximumSubmitDepth, signed_size_t depth) {

      KdNodePDB* node = nullptr;

      // The partition permutes as x, y, z, w... and specifies the most significant key.
      size_t p = permutation.at(depth).at(permutation.at(0).size() - 1);

      // Obtain the reference array that corresponds to the most significant key.
      KdNodePDB** reference = references[permutation.at(depth).at(N)];

      if (end == start) {

        // Only one reference was passed to this function, so add it to the tree.
        node = reference[end];

      }
      else if (end == start + 1) {

        // Two references were passed to this function in sorted order, so store the start
        // element at this level of the tree and store the end element as the > child.
        node = reference[start];
        node->gtChild = reference[end];

      }
      else if (end == start + 2) {

        // Three references were passed to this function in sorted order, so
        // store the median element at this level of the tree, store the start
        // element as the < child and store the end element as the > child.
        node = reference[start + 1];
        node->ltChild = reference[start];
        node->gtChild = reference[end];

      }
      else if (end > start + 2) {

        // Four or more references were passed to this function, so the
        // median element of the reference array is chosen as the tuple
        // about which the other reference arrays will be partitioned
        // Avoid overflow when computing the median.
        size_t median = start + ((end - start) / 2);

        // Store the median element of the reference array in a new KdNode.
        node = reference[median];

        // Build both branches with child threads at as many levels of the tree
        // as possible.  Create the child threads as high in the tree as possible.
        // Are child threads available to build both branches of the tree?
        if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

          // No, child threads are not available, so one thread will be used.
          // Initialize startIndex=1 so that the 'for' loop that partitions the
          // reference arrays will partition a number of arrays equal to N.
          size_t startIndex = 1;

          // If depth < N-1, copy references[permute[N]] to references[permute[0]]
          // where permute is the permutation vector for this level of the tree.
          // Sort the two halves of references[permut[0]] with p+1 as the most
          // significant key of the super key. Use as the temporary array
          // references[permute[1]] because that array is not used for partitioning.
          // Partition a number of reference arrays equal to the tree depth because
          // those reference arrays are already sorted.
          if (depth < N - 1) {
            startIndex = N - depth;
            KdNodePDB** dst = references[permutation.at(depth).at(0)];
            for (size_t i = start; i <= end; ++i) {
              dst[i] = reference[i];
            }
            // Sort the lower half of references[permute[0]] with the current thread.
            size_t numThreads = 1;
            size_t pp1 = p + 1;
            if ((maximumSubmitDepth - depth) > 0) numThreads = (1LL << (maximumSubmitDepth - depth));
            parallelSort(dst + start, dst + median, [&](KdNodePDB* a, KdNodePDB* b) {
              return 0 > superKeyCompare(a->tuple, b->tuple, pp1);
              }, numThreads);

            // Sort the upper half of references[permute[0]] with the current thread.
            parallelSort(dst + median + 1, dst + end + 1, [&](KdNodePDB* a, KdNodePDB* b) {
              return 0 > superKeyCompare(a->tuple, b->tuple, pp1);
              }, numThreads);
          }

          // Partition the reference arrays specified by 'startIndex' in
          // a priori sorted order by comparing super keys.  Store the
          // result from references[permute[i]] in references[permute[i-1]]
          // where permute is the permutation vector for this level of the
          // tree, thus permuting the reference arrays. Skip the element
          // of references[permute[i]] that equals the tuple that is stored
          // in the new KdNode.
          K* tuple = node->tuple;
          for (size_t i = startIndex; i < N; ++i) {
            // Specify the source and destination reference arrays.
            KdNodePDB** src = references[permutation.at(depth).at(i)];
            KdNodePDB** dst = references[permutation.at(depth).at(i - 1)];

            // Fill the lower and upper halves of one reference array
            // in ascending order with the current thread.
            for (size_t j = start, lower = start - 1, upper = median; j <= end; ++j) {
              KdNodePDB* src_j = src[j];
              K compare = superKeyCompare(src_j->tuple, tuple, p);
              if (compare < 0) {
                dst[++lower] = src_j;
              }
              else if (compare > 0) {
                dst[++upper] = src_j;
              }
            }
          }

          // Recursively build the < branch of the tree with the current thread.
          node->ltChild = buildKdTree(references, permutation, start, median - 1,
            maximumSubmitDepth, depth + 1);

          // Then recursively build the > branch of the tree with the current thread.
          node->gtChild = buildKdTree(references, permutation, median + 1, end,
            maximumSubmitDepth, depth + 1);

        }
        else {

          // Yes, child threads are available, so two threads will be used.
          // Initialize endIndex=0 so that the 'for' loop that partitions the
          // reference arrays will partition a number of arrays equal to N.
          size_t startIndex = 1;

          // If depth < N-1, copy references[permute[N]] to references[permute[0]]
          // where permute is the permutation vector for this level of the tree.
          // Sort the two halves of references[permute[0]] with p+1 as the most
          // significant key of the super key. Use as the temporary array
          // references[permute[1]] because that array is not used for partitioning.
          // Partition a number of reference arrays equal to the tree depth because
          // those reference arrays are already sorted.
          if (depth < N - 1) {
            startIndex = N - depth;
            KdNodePDB** dst = references[permutation.at(depth).at(0)];
            // Copy and sort the lower half of references[permute[0]] with a child thread.
            std::future<void> copyFuture =
              std::async(std::launch::async, [&] {
              for (size_t i = start; i <= median - 1L; ++i) {
                dst[i] = reference[i];
              }
              size_t numThreads = 1;
              size_t pp1 = p + 1;
              if (maximumSubmitDepth > -1) numThreads = (1LL << (maximumSubmitDepth - depth + 1));
              parallelSort(dst + start, dst + median, [&](KdNodePDB* a, KdNodePDB* b) {
                return 0 > superKeyCompare(a->tuple, b->tuple, pp1);
                }, numThreads / 2);
                });

            // Copy and sort the upper half of references[permute[0]] with the current thread.
            for (size_t i = median + 1; i <= end; ++i) {
              dst[i] = reference[i];
            }
            int64_t numThreads = 1;
            size_t pp1 = p + 1;
            if (maximumSubmitDepth > -1) numThreads = (1LL << (maximumSubmitDepth - depth + 1));
            parallelSort(dst + median + 1, dst + end + 1, [&](KdNodePDB* a, KdNodePDB* b) {
              return 0 > superKeyCompare(a->tuple, b->tuple, pp1);
              }, numThreads / 2);

            // Wait for the child thread to finish execution.
            try {
              copyFuture.get();
            }
            catch (std::exception const& e) {
              std::cout << "caught exception " << e.what() << std::endl;
            }
          }

          // Create a copy of the node->tuple array so that the current thread
          // and the child thread do not contend for read access to this array.
          K* tuple = node->tuple;
          K point[N];
          for (size_t i = 0; i < N; ++i) {
            point[i] = tuple[i];
          }

          // Partition the reference arrays specified by 'startIndex' in
          // a priori sorted order by comparing super keys.  Store the
          // result from references[permute[i]]] in references[permute[i-1]]
          // where permute is the permutation vector for this level of the
          // tree, thus permuting the reference arrays. Skip the element
          // of references[permute[i]] that equals the tuple that is stored
          // in the new KdNode.
          for (size_t i = startIndex; i < N; ++i) {
            // Specify the source and destination reference arrays.
            KdNodePDB** src = references[permutation.at(depth).at(i)];
            KdNodePDB** dst = references[permutation.at(depth).at(i - 1)];

            // Two threads may be used to partition the reference arrays, analogous to
            // the use of two threads to merge the results for the merge sort algorithm.
            // Fill one reference array in ascending order with a child thread.
            std::future<void> partitionFuture =
              std::async(std::launch::async, [&] {
              for (size_t lower = start - 1, upper = median, j = start; j <= median; ++j) {
                KdNodePDB* src_j = src[j];
                K compare = superKeyCompare(src_j->tuple, point, p);
                if (compare < 0) {
                  dst[++lower] = src_j;
                }
                else if (compare > 0) {
                  dst[++upper] = src_j;
                }
              }
                });

            // Simultaneously fill the same reference array in descending order with the current thread.
            for (size_t lower = median, upper = end + 1, k = end; k > median; --k) {
              KdNodePDB* src_k = src[k];
              K compare = superKeyCompare(src_k->tuple, tuple, p);
              if (compare < 0) {
                dst[--lower] = src_k;
              }
              else if (compare > 0) {
                dst[--upper] = src_k;
              }
            }

            // Wait for the child thread to finish execution.
            try {
              partitionFuture.get();
            }
            catch (std::exception const& e) {
              std::cout << "caught exception " << e.what() << std::endl;
            }
          }

          // Recursively build the < branch of the tree with a child thread.
          // The recursive call to buildKdTree must be placed in a lambda
          // expression because buildKdTree is a template not a function.
          std::future<KdNodePDB*> buildFuture =
            std::async(std::launch::async, [&] {
            return buildKdTree(references, permutation, start, median - 1,
              maximumSubmitDepth, depth + 1);
              });

          // And simultaneously build the > branch of the tree with the current thread.
          node->gtChild = buildKdTree(references, permutation, median + 1, end,
            maximumSubmitDepth, depth + 1);

          // Wait for the child thread to finish execution.
          try {
            node->ltChild = buildFuture.get();
          }
          catch (std::exception const& e) {
            std::cout << "caught exception " << e.what() << std::endl;
          }
        }

      }
      else if (end < start) {

        // This is an illegal condition that should never occur, so test for it last.
        std::cout << "error has occurred at depth = " << depth << " : end = " << end
          << "  <  start = " << start << std::endl;
        exit(1);

      }

      // Return the pointer to the root of the k-d tree.
      return node;
    }


    /*
      * The swap function swaps two elements in a vector<size_t>.
      *
      * Calling parameters:
      *
      * a - the vector
      * i - the index of the first element
      * j - the index of the second element
      */
  private:
    inline
      static void swap(std::vector<size_t>& a, size_t i, size_t j) {
      size_t t = a.at(i);
      a.at(i) = a.at(j);
      a.at(j) = t;
    }

    /*
      * The createKdTree function performs the necessary initialization then calls the buildKdTree function.
      *
      * Calling parameters:
      *
      * kdNodes - a vector<KdNode*> wherein each KdNode contains a (x,y,z,w...) tuple
      * numDimensions - the number of dimensions
      * numThreads - the number of threads
      * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
      *
      * returns: a KdNode pointer to the root of the k-d tree
      */
  public:
    static KdNodePDB* createKdTree(std::vector<KdNodePDB*>& kdNodes,
      size_t numThreads, signed_size_t maximumSubmitDepth) {

      // Create the references arrays including one additional array for use in building the k-d tree.
      KdNodePDB*** references = new KdNodePDB * *[N + 1];

      // The first references array is the .data() array of the kdNodes vector.
      references[0] = kdNodes.data();

      // Allocate the remaining references arrays.
      for (size_t i = 1; i < N + 1; ++i) {
        references[i] = new KdNodePDB * [kdNodes.size()];
      }

      // Sort the first reference array using multiple threads. Importantly,
      // for compatibility with the 'permutation' vector initialized below,
      // use the first dimension (0) as the leading key of the super key.
      size_t numSortThreads = 1;
      if (maximumSubmitDepth >= 0) numSortThreads = (1ULL << (maximumSubmitDepth + 1));
      parallelSort(references[0], references[0] + kdNodes.size(),
        [&](KdNodePDB* a, KdNodePDB* b) {
          return (0 > superKeyCompare(a->tuple, b->tuple, 0));
        }, numSortThreads);

      // Remove references to duplicate coordinates via one pass through the first reference array.
      size_t end = removeDuplicates(references[0], 0, kdNodes.size());


      // Determine the maximum depth of the k-d tree, which is log2( kdNodes.size() ).
      size_t size = kdNodes.size();
      size_t maxDepth = 1;
      while (size > 0) {
        ++maxDepth;
        size >>= 1;
      }

      // It is unnecessary to compute either the permutation of the reference array or
      // the partition coordinate upon each recursive call of the buildKdTree function
      // because both depend only on the depth of recursion, so they may be pre-computed.
      // Create and initialize an 'indices' vector for the permutation calculation.
      // Because this vector is initialized with 0, 1, 2, 3, 0, 1, 2, 3, etc. (for
      // e.g. 4-dimensional data), the leading key of the super key will be 0 at the
      // first level of the nascent tree, consistent with having sorted the reference
      // array above using 0 as the leading key of the super key.
      std::vector<size_t> indices(N + 2);
      for (size_t i = 0; i < indices.size() - 1; ++i) {
        indices[i] = i;
      }

      // Create a 'permutation' vector from the 'indices' vector to specify permutation
      // of the reference arrays and of the partition coordinate.
      std::vector< std::vector<size_t> > permutation(maxDepth, std::vector<size_t>(N + 2));

      // Fill the permutation vector by calculating the permutation of the indices vector
      // and the the partition coordinate of the tuple at each depth in the tree.
      for (size_t i = 0; i < permutation.size(); ++i) {
        // The last entry of the indices vector contains the partition coordinate.
        indices.at(N + 1) = i % N;
        // Swap the first and second to the last elements of the indices vector.
        swap(indices, 0, N);
        // Copy the indices vector to one row of the permutation vector.
        permutation.at(i) = indices;
        // Swap the third and second to the last elements of the indices vector.
        swap(indices, N - 1, N);
      }

      // Build the k-d tree with multiple threads if possible.
      KdNodePDB* root = buildKdTree(references, permutation, 0, end,
        maximumSubmitDepth, 0);
      // Delete all but the first of the references arrays.
      for (size_t i = 1; i < N + 1; ++i) {
        delete[] references[i];
      }
      delete[] references;

      // Return the pointer to the root of the k-d tree.
      return root;
    }

// The following two macros are used as pairs to protect critical code from
// being accessed at the same time by different thread.  These macros 
// allow different means of preforming that locks to be experimented with.
#define LOCK_THIS_BLOCK thisMutex.lock()
#define UNLOCK_THIS_BLOCK thisMutex.unlock() 
//#define LOCK_THIS_BLOCK 
//#define UNLOCK_THIS_BLOCK 

    /*
     * The checkAllTakenStatus function checks the all taken status of this
     * node for the moveToID.
     *
     * Calling parameters:
     *
     * moveToID - the ID of the container currently being processed.
     *
     * return - a SearchReturn enum indicating whether all the nodes 
     *          in the tree below have already been taken by this cluster.
     */

  private:
    inline SearchRet checkAllTakenStatus(const uint32_t moveToID) {
      // figure out if there is any need to proceed further
      if (allTakenSet.size() > 0) {
        LOCK_THIS_BLOCK;
        if (allTakenSet.contains(moveToID)) {
          UNLOCK_THIS_BLOCK;
          return AllTakenByID;
        }
        UNLOCK_THIS_BLOCK;
        return AllTaken;
      }
      return NoResult;
    }


    /*
     * The regionSearch function searches the k-d tree to find the KdNodes that
     * lie within a hyper-rectangle defined by the query lower and upper bounds.
     *
     * Calling parameters:
     *
     * result:      reference to a list of where to put pointers to the tuples of the taken points.
     * queryLower:  the query lower bound vector
     * queryUpper:  the query upper bound vector
     * moveTo:      pointer to the cluster into which the found items will be placed
     * permutation: vector that specifies permutation of the partition coordinate
     * maximumSubmitDepth: the maximum tree depth at which a child task may be launched
     * depth:       the current depth in the k-d tree
     *
     * return - a SearchReturn enum indicating whether all the nodes in the tree below have
     *          already been taken by this cluster.
     */
  public:
    SearchRet regionSearchAndRemove(SimpleList<K*>& result,
      const std::vector<K>& queryLower, const std::vector<K>& queryUpper,
      M* const moveTo, const signed_size_t depth, const std::vector<int64_t>& permutation, size_t treeIdx) {

      // get a local copy of the cluster ID;
      const uint32_t moveToID = moveTo->getID();

      //nodesVisited++;

      bool ltFlag = false;
      bool gtFlag = false;

      // Look up the primary coordinate index.
      unsigned int p = (unsigned int)permutation[(unsigned int)depth];

      // These return codes indicate whether the status of the node that was just returned from
      // init the return result to NoResult
      SearchRet ltRetCode = Undefined;
      SearchRet gtRetCode = Undefined;
      bool inside = false;

      // Search the < branch of the k-d tree if the partition coordinate of the queryPlus is
      // <= the partition coordinate of the k-d node.  The < branch
      // must be searched when the cutoff distance equals the partition coordinate because the super
      // key may assign a point to either branch of the tree if the sorting or partition coordinate,
      // which forms the most significant portion of the super key, shows equality.
      if (queryLower[p] <= tuple[p]) {
        // but only search if the ltChild pointer is not null 
        if (ltChild != nullptr) {
          ltRetCode = ltChild->checkAllTakenStatus(moveToID);
          if (ltRetCode != AllTakenByID)
            ltFlag = true;// defer going doen the lt branch until after this nodes point have veen added to the result FIFO
        }
        else
          ltRetCode = AllTakenByID;
        // Search the > branch of the k-d tree if the partition coordinate of the queryPlus is
        // >= the partition coordinate of the k-d node.  The < branch
        // must be searched when the cutoff distance equals the partition coordinate because the super
        // key may assign a point to either branch of the tree if the sorting or partition coordinate,
        // which forms the most significant portion of the super key, shows equality.
        if (queryUpper[p] >= tuple[p]) {
          // Only search if the gtChild pointer is not null.  Otherwise report the return code as all taken.
          if (gtChild != nullptr) {
            gtRetCode = gtChild->checkAllTakenStatus(moveToID);
            if (gtRetCode != AllTakenByID)
              gtFlag = true; // defer going doen the gt branch until after this nodes point have veen added to the result FIFO
          }
          else
            gtRetCode = AllTakenByID;
          // Now check to see if the local tuple is inside the the hypercube.
          // the test to down both children is a necessary but not sufficient indication of being inside.
          // If the distance from the query node to the k-d node is within the query rectangle
          // in all k dimensions, add the tuple to the return list an add the values to the cluster.
          if (moveTo != movedTo) {
            inside = true;
            for (size_t i = 0; i < queryUpper.size(); i++) {
              if ((queryUpper[i] < tuple[i]) || (queryLower[i] > tuple[i])) {
                inside = false;
                break;
              }
            }
            
            if (inside) {
              // check if we need to grab a lock on this node
              LOCK_THIS_BLOCK;
              // check if the node had not already been taken after we have the lock
              if (movedTo == nullptr) {
                // if it has not, add this nodes data to the return list and cluster
                result.push_back(tuple);
                moveTo->addToCluster(tuple, *values);
                movedTo = moveTo;  // mark the value as taken by this cluster.
              }
              else {
                // this node was taken by another cluster, so indicate this cluster overlaps 
                // that cluster by adding it to the overlaps set.
                if (moveTo->overlaps == nullptr) {
                  moveTo->overlaps = new avlSet<M*>;
                }
                auto rslt = moveTo->overlaps->insert(movedTo);
              }
              UNLOCK_THIS_BLOCK;
            }
          }
        }
        else {
          // since we did not go down the gtChild path, get that nodes taken status.
          if (gtChild != nullptr)
            gtRetCode = gtChild->checkAllTakenStatus(moveToID);
          else
            gtRetCode = AllTakenByID;
        }
      }
      else { // Since we did not go down the < path, check the > path.
        if (queryUpper[p] >= tuple[p]) {
          // but only search if the gtChild pointer is not null and all nodes below are not taken for this cluster
          if (gtChild != nullptr) {
            gtRetCode = gtChild->checkAllTakenStatus(moveToID);
            if (gtRetCode != AllTakenByID)
              gtFlag = true; // defer going doen the lt branch until later
          }
          else
            gtRetCode = AllTakenByID;
          if (ltChild != nullptr)
            ltRetCode = ltChild->checkAllTakenStatus(moveToID);
          else
            ltRetCode = AllTakenByID;
        }
      }

      if (ltFlag) // now perform the search of the lt branch
        ltRetCode = ltChild->regionSearchAndRemove(result, queryLower, queryUpper, moveTo, depth + 1, permutation, treeIdx << 1);
      if (gtFlag) // now perform the search of the gt branch
        gtRetCode = gtChild->regionSearchAndRemove(result, queryLower, queryUpper, moveTo, depth + 1, permutation, (treeIdx << 1) + 1);

      // now figure out what to record in the allTakenSet and what code to return
      // this first condition covers the case where the current and all nodes below are have been taken by the same ID. 
      if (movedTo == moveTo && ltRetCode == AllTakenByID && gtRetCode == AllTakenByID) {
        LOCK_THIS_BLOCK;
        allTakenSet.insert(moveToID);
        UNLOCK_THIS_BLOCK;
        return AllTakenByID;
      }
      // this second condition covers the case where all the nodes below have been have been or would have been taken by this ID.
      if (inside && ltRetCode == AllTakenByID && gtRetCode == AllTakenByID) {
        LOCK_THIS_BLOCK;
        allTakenSet.insert(moveToID);
        UNLOCK_THIS_BLOCK;
        return AllTakenByID;
      }
      // this case handles the condition where  this node and the children have been taken by some cluster but not all the same ID. 
      if (movedTo != nullptr && ltRetCode != NoResult && gtRetCode != NoResult) {
        LOCK_THIS_BLOCK;
        allTakenSet.insert(allTakenID);
        UNLOCK_THIS_BLOCK;
        return AllTaken;
      }
      // if none of the above return nothing to change.
      return NoResult;
    }

    void allTakenCheck(bool& allTaken, M* moveTo) {
      if (movedTo != moveTo) {
        allTaken = false;
        return;
      }
      if (ltChild != nullptr && allTaken) {
        ltChild->allTakenCheck(allTaken, moveTo);
      }
      if (gtChild != nullptr && allTaken) {
        gtChild->allTakenCheck(allTaken, moveTo);
      }
     }

    /*
    The pickValue picks a value from the kdTree by descending to the lowest point in the tree that it finds
    a node where a value has not been taken yet.
    Parameters:
      result:     reference to a list of where to put a pointer to the tuple of the picked point.
      selector:   unsigned uint64_t where the bits are used to guide which child node to descend
      moveTo:     pointer to cluster where that will contain the value picked
      depth:      integer counting the number of levels descended.

      return - a SearchReturn enum indicating whether all the nodes in the tree below have
                already been taken by this cluster.
    **/

  public:
    SearchRet pickValue(SimpleList<K*>& result, const uint64_t selector, M* const moveTo, const int depth) {

      // get a local copy of the container ID
      const uint32_t& moveToID = moveTo->getID();

      // figure out if there is any need to proceed further down the tree.
      // if this nodes all taken set has any values it in, no further search is required.
     if (allTakenSet.size() > 0) {
        return AllTaken;
      }

      // init the return codes to undefined.
      SearchRet ltRetCode = Undefined;
      SearchRet gtRetCode = Undefined;
      // and inside to false.
      bool inside = false;

      // set the preferred recursion direction for this level from the selector
      bool goGtThan = (selector & 0x1ULL) == 1;

      // If preferred is lt or the preferred is gt but that path is blocked, 
      // continue the search on the appropriate the lt child
      if (!goGtThan) {
        if (ltChild != nullptr) {
          ltRetCode = ltChild->pickValue(result, selector >> 1, moveTo, depth + 1);
        }
        else {
          ltRetCode = AllTakenByID;
        }
        if (result.size() == 0LLU) {
          if (gtChild != nullptr) {
            gtRetCode = gtChild->pickValue(result, selector >> 1, moveTo, depth + 1);
          }
          else {
            gtRetCode = AllTakenByID;
          }
        }
      }
      else {
        if (gtChild != nullptr) {
          gtRetCode = gtChild->pickValue(result, selector >> 1, moveTo, depth + 1);
        }
        else {
          gtRetCode = AllTakenByID;
        }
        if (result.size() == 0LLU) {
          if (ltChild != nullptr) {
            ltRetCode = ltChild->pickValue(result, selector >> 1, moveTo, depth + 1);
          }
          else {
            ltRetCode = AllTakenByID;
          }
        }
      }

      if (ltRetCode == Undefined) { // indicating the lt child was not called above.
        if (ltChild != nullptr)
          ltRetCode = ltChild->allTakenSet.size() == 0 ? NoResult : AllTaken;
        else
          ltRetCode = AllTakenByID;
      }
      if (gtRetCode == Undefined) {// indicating the gt child was not called above.
        if (gtChild != nullptr) 
          gtRetCode = gtChild->allTakenSet.size() == 0 ? NoResult : AllTaken;
        else
          gtRetCode = AllTakenByID;
      }


      // if no result was found on either of the children, see if this value is available.
      if (result.size() == 0LLU) {  // have not found a value yet
        if (movedTo == nullptr) {  // see if we should grab a lock
          // grab a lock
          LOCK_THIS_BLOCK;
          // check again to make sure this node was not taken
          if (movedTo == nullptr) {
            // capture the data from this node.
            result.push_back(tuple);
            moveTo->addToCluster(tuple, *values);
            // mark the node as taken by the searching container.
            movedTo = moveTo;  
            inside = true;
          }
          UNLOCK_THIS_BLOCK;
        }
      }

      // now figure out what to record in the allTakenSet and what code to return
      // this first condition covers the case where the current and all nodes below are have been taken by the same ID. 
      if (movedTo == moveTo && ltRetCode == AllTakenByID && gtRetCode == AllTakenByID) {
        LOCK_THIS_BLOCK;
        allTakenSet.insert(moveToID);
        UNLOCK_THIS_BLOCK;
        return AllTakenByID;
      }
      // this case handles the condition where all this node and the children have been taken by something but not all the same ID. 
      if (movedTo != nullptr && ltRetCode != NoResult && gtRetCode != NoResult) {
        LOCK_THIS_BLOCK;
        allTakenSet.insert(allTakenID);
        UNLOCK_THIS_BLOCK;
        return AllTaken;
      }
      // if none of the above return nothing to change.
      return NoResult;
    }

  }; // class KdNodePDB


private:
  size_t numPoints = 0;  // total number of points submitted to the KdTree.
  std::vector<KdNodePDB*> kdNodePtrs;  // vector of kdNodes resulting from submissions to the KdTree.
  SimpleAlloc simpleAlloc;  // Local memory allocator.  
  STLAdaptor<V>* adaptorV;  // pointer to an adaptor that wraps simpleAlloc for use with STL containers of type V
  KdNodePDB* root = nullptr;  // root of the KdTree after it's built
  std::vector<int64_t> permutation; // vector of pre-calculated levels
  signed_size_t numThreads = 1;     // number of threads used in the process
  signed_size_t maximumSubmitDepth = -1; // maximum depth in tree building or searching where new threads are submitted.

  // calculates the maximumSubmitDepth given the number of threads.
  signed_size_t calcMaximumSubmitDepth() {
    signed_size_t n = 0;
    if (numThreads > 0) {
      while (numThreads > 0) {
        n++;
        numThreads >>= 1L;
      }
      numThreads = 1LL << (n - 1);
    }
    else {
      numThreads = 0;
    }
    signed_size_t childThreads = numThreads - 1LL;
    maximumSubmitDepth = -1;
    if (numThreads < 2) {
      maximumSubmitDepth = -1; // The sentinel value -1 specifies no child threads.
    }
    else if (numThreads == 2) {
      maximumSubmitDepth = 0;
    }
    else {
      maximumSubmitDepth = (signed_size_t)floor(log((double)childThreads) / log(2.));
    }
    return maximumSubmitDepth;
  }

  std::vector<int64_t> getPermutations(size_t size, size_t dimension) {
    // Determine the maximum depth of the k-d tree, which is log2(size).
    unsigned int maxDepth = 1;
    int64_t lsize = size;
    while (lsize > 0) {
      maxDepth++;
      lsize >>= 1;
    }
    // The partition coordinate permutes n the order 0, 1, 2, 3, 0, 1, 2, 3, etc.
    // for e.g. 4-dimensional data.
    permutation = std::vector<int64_t>(maxDepth);
    for (size_t i = 0; i < permutation.size(); ++i) {
      permutation.at(i) = i % dimension;
    }
    return permutation;
  }

public:
  // main constructor
  KdTreePDB() {
    numThreads = std::thread::hardware_concurrency();
    calcMaximumSubmitDepth();
    adaptorV  = new STLAdaptor<V>(simpleAlloc);
  }

  // sets the number of threads to use but that number is forced to the next lower power of 2.
  void setNumThreads(signed_size_t thd) {
    numThreads = thd;
    calcMaximumSubmitDepth();
  }

  // add points and values to the KdTree.
  size_t add(std::vector<K>& tuple, V value) {
    if (tuple.size() != N) {
      return 0;
    }
    // allocate a KdNode from the local allocator
    // provide the allocator and adaptor 
    void* tp = simpleAlloc.allocate(sizeof(KdNodePDB), alignof(KdNodePDB));
    auto kp = new(tp) KdNodePDB(tuple, value, adaptorV, &simpleAlloc);
    // add a KdNode pointer to the pointer list.
    kdNodePtrs.push_back(kp);
    numPoints = kdNodePtrs.size();
    return numPoints;
  }

  /*
 * The searchRegionAndRemove function searches the k-d tree to find the KdNodes that
 * lie within a hyper-rectangle defined by the query lower and upper bounds.
 *
 * Calling parameters:
 *
 * result:      reference to a list of where to put pointers to the tuples of the taken points.
 * queryLower:  the query lower bound vector
 * queryUpper:  the query upper bound vector
 * moveTo:      pointer to the cluster into which the found items will be placed
 *
 * return:      a Boolean that is always true;
 */

  bool searchRegionAndRemove(SimpleList<K*>& result, std::vector<K>& queryLower, std::vector<K>& queryUpper, M* moveTo) {
    // if the tree is not built yet, build it
    if (root == nullptr) {
      buildTree();
    }

    // Ensure that each query lower bound <= the corresponding query upper bound.
    for (size_t i = 0; i < queryLower.size(); ++i) {
      if (queryLower[i] > queryUpper[i]) {
        K tmp = queryLower[i];
        queryLower[i] = queryUpper[i];
        queryUpper[i] = tmp;
      }
    }

    // Search the tree and return the resulting list of KdNodes.
    root->regionSearchAndRemove(result, queryLower, queryUpper, moveTo, 0LL, permutation, 1ULL);
    return true;
  }


  bool buildTree(signed_size_t numThreads = -1) {
    if (numThreads != -1) {
      setNumThreads(numThreads);
    }
    permutation = getPermutations(numPoints, N);
    root = KdNodePDB::createKdTree(kdNodePtrs, numThreads, maximumSubmitDepth);
    return true;
  }

  /*
  * The pickValue picks a value from the kdTree by descending to the lowest point in the tree that it finds
  * a node where a value has not been taken yet.
  * Parameters:
  *   result:     reference to a list of where to put a pointer to the tuple of the picked point.
  *   moveTo:     pointer to cluster where that will contain the value picked

  *   return - a Boolean indicating whether it found a node.
  **/
  bool pickValue(SimpleList<K*>& result, uint64_t selectionBias, M* moveTo) {
    // if the tree is not built yet, build it
    if (root == nullptr) {
      buildTree();
      // if root is still null; return a null
      if (root == nullptr) return false;
    }
    // search the tree to find the node to return and possibly delete.
    {
      SearchRet returnResult = root->pickValue(result, selectionBias, moveTo, 0);
    }
    // Check to see if no data was returned.  If so then return false.
    if (result.size() == 0) return false;

    return true;
  }

  ~KdTreePDB() {
    //Note: the KdNodes are deleted when the local simpleAlloc allocator is deleted
    //std::cout << nodesVisited << " Nodes visited\n";
    delete adaptorV;
  }

}; // KdTree

#endif // KDTREE_HPP
