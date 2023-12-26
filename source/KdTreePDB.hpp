/*
 * Copyright (c) 2023 John Robinson.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

//
//  KdTreePDB.h
//  KdTree
//
//  Created by John Robinson on 11/18/23.
//

#ifndef KDTREE_HPP
#define KDTREE_HPP
#include <random>
#include <limits>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <list>
#include <set>
#include <iostream>
#include <iomanip>
#include <exception>
#include <future>
#include "parallelSort.hpp"

/*
 * This type is the signed equivalent of size_t and might be equivalent to intmax_t
 */
typedef int64_t signed_size_t;

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

enum SearchRet { // the following are codes that the search routines will return fom a on the status of then node visited.
  DeadNode,                // indicates the values, ltChile and gtChild pointer are all null.
  NoResult,                 // indicates to result has been returned.
  ResultFound,              // indicates a reulst has been found.
  ResultFoundAndAllTaken,   // indicates reuls has been found and all nodes below that node have been taken
  AllTaken                  // indicates a result has not been found but all nodes below have been taken
};


template<typename K, typename V, typename M, size_t N>
class KdNodePDB {

  typedef std::pair<std::vector<K>, std::list<V>*> retPair_t;

public:
  K tuple[N];
private:
  KdNodePDB<K, V, M, N>* ltChild;
  KdNodePDB<K, V, M, N>* gtChild;
  std::list<V>* values = nullptr;
  std::mutex thisMutex;
  M* movedTo = nullptr;        // pointer to the container that is holding the values that were here.
  uint32_t  gtAllTaken = 0;    // ID of the container that has all values from the child nodes below
  uint32_t  ltAllTaken = 0;    // ID of the container that has all values from the child nodes below

public:
  KdNodePDB(V const value) { // Pass non-primitive types as 'V const&'
    for (size_t i = 0; i < N; i++) this->tuple[i] = (K)0;
    ltChild = gtChild = nullptr; // redundant
    values = new std::list<V>();
    values->push_back(value);
  }

  KdNodePDB(std::vector<K>& tuple, V const value) { // Pass non-primitive types as 'V const&'
    for (size_t i = 0; i < N; i++) this->tuple[i] = tuple[i];
    ltChild = gtChild = nullptr; // redundant
      values = new std::list<V>();
    values->push_back(value);
  }

public:
  ~KdNodePDB() {
    delete values;
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
    static size_t removeDuplicates(KdNodePDB<K, V, M, N>** kdNodes, size_t i, size_t size) {
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
        // append vales to the end of fte kpet kdNode and delete the now unused node
        kdNodes[end]->values->splice(kdNodes[end]->values->end(), *kdNodes[j]->values);
        delete kdNodes[j];
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
  static KdNodePDB<K, V, M, N>* buildKdTree(KdNodePDB<K, V, M, N>*** references,
    std::vector< std::vector<size_t> > const& permutation,
    size_t start, size_t end,
    signed_size_t maximumSubmitDepth, signed_size_t depth) {

    KdNodePDB<K, V, M, N>* node = nullptr;

    // The partition permutes as x, y, z, w... and specifies the most significant key.
    size_t p = permutation.at(depth).at(permutation.at(0).size() - 1);

    // Obtain the reference array that corresponds to the most significant key.
    KdNodePDB<K, V, M, N>** reference = references[permutation.at(depth).at(N)];

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

        // If depth < N-1, copy references[permut[N]] to references[permut[0]]
        // where permut is the permutation vector for this level of the tree.
        // Sort the two halves of references[permut[0]] with p+1 as the most
        // significant key of the super key. Use as the temporary array
        // references[permut[1]] because that array is not used for partitioning.
        // Partition a number of reference arrays equal to the tree depth because
        // those reference arrays are already sorted.
        if (depth < N - 1) {
          startIndex = N - depth;
          KdNodePDB<K, V, M, N>** dst = references[permutation.at(depth).at(0)];
           for (size_t i = start; i <= end; ++i) {
            dst[i] = reference[i];
          }
          // Sort the lower half of references[permut[0]] with the current thread.
          size_t numthreads = 1;
          size_t pp1 = p + 1;
          if ((maximumSubmitDepth - depth) > 0) numthreads = (1LL << (maximumSubmitDepth - depth));
          parallelSort(dst + start, dst + median, [&](KdNodePDB<K, V, M, N>* a, KdNodePDB<K, V, M, N>* b) {
            return 0 > superKeyCompare(a->tuple, b->tuple, pp1);
            }, numthreads);

          // Sort the upper half of references[permut[0]] with the current thread.
          parallelSort(dst + median + 1, dst + end + 1, [&](KdNodePDB<K, V, M, N>* a, KdNodePDB<K, V, M, N>* b) {
            return 0 > superKeyCompare(a->tuple, b->tuple, pp1);
            }, numthreads);
        }

        // Partition the reference arrays specified by 'startIndex' in
        // a priori sorted order by comparing super keys.  Store the
        // result from references[permut[i]]] in references[permut[i-1]]
        // where permut is the permutation vector for this level of the
        // tree, thus permuting the reference arrays. Skip the element
        // of references[permut[i]] that equals the tuple that is stored
        // in the new KdNode.
        K* tuple = node->tuple;
        for (size_t i = startIndex; i < N; ++i) {
          // Specify the source and destination reference arrays.
          KdNodePDB<K, V, M, N>** src = references[permutation.at(depth).at(i)];
          KdNodePDB<K, V, M, N>** dst = references[permutation.at(depth).at(i - 1)];

          // Fill the lower and upper halves of one reference array
          // in ascending order with the current thread.
          for (size_t j = start, lower = start - 1, upper = median; j <= end; ++j) {
            KdNodePDB<K, V, M, N>* src_j = src[j];
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

        // If depth < N-1, copy references[permut[N]] to references[permut[0]]
        // where permut is the permutation vector for this level of the tree.
        // Sort the two halves of references[permut[0]] with p+1 as the most
        // significant key of the super key. Use as the temporary array
        // references[permut[1]] because that array is not used for partitioning.
        // Partition a number of reference arrays equal to the tree depth because
        // those reference arrays are already sorted.
        if (depth < N - 1) {
          startIndex = N - depth;
          KdNodePDB<K, V, M, N>** dst = references[permutation.at(depth).at(0)];
          // Copy and sort the lower half of references[permut[0]] with a child thread.
          std::future<void> copyFuture =
            std::async(std::launch::async, [&] {
            for (size_t i = start; i <= median - 1L; ++i) {
              dst[i] = reference[i];
            }
            size_t numthreads = 1;
            size_t pp1 = p + 1;
            if (maximumSubmitDepth > -1) numthreads = (1LL << (maximumSubmitDepth - depth + 1));
            parallelSort(dst + start, dst + median, [&](KdNodePDB<K, V, M, N>* a, KdNodePDB<K, V, M, N>* b) {
              return 0 > superKeyCompare(a->tuple, b->tuple, pp1);
              }, numthreads / 2);
          });

          // Copy and sort the upper half of references[permut[0]] with the current thread.
          for (size_t i = median + 1; i <= end; ++i) {
            dst[i] = reference[i];
          }
          int64_t numthreads = 1;
          size_t pp1 = p + 1;
          if (maximumSubmitDepth > -1) numthreads = (1LL << (maximumSubmitDepth - depth + 1));
          parallelSort(dst + median + 1, dst + end + 1, [&](KdNodePDB<K, V, M, N>* a, KdNodePDB<K, V, M, N>* b) {
            return 0 > superKeyCompare(a->tuple, b->tuple, pp1);
            }, numthreads / 2);

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
        K* point = new K[N];
        for (size_t i = 0; i < N; ++i) {
          point[i] = tuple[i];
        }

        // Partition the reference arrays specified by 'startIndex' in
        // a priori sorted order by comparing super keys.  Store the
        // result from references[permut[i]]] in references[permut[i-1]]
        // where permut is the permutation vector for this level of the
        // tree, thus permuting the reference arrays. Skip the element
        // of references[permut[i]] that equals the tuple that is stored
        // in the new KdNode.
        for (size_t i = startIndex; i < N; ++i) {
          // Specify the source and destination reference arrays.
          KdNodePDB<K, V, M, N>** src = references[permutation.at(depth).at(i)];
          KdNodePDB<K, V, M, N>** dst = references[permutation.at(depth).at(i - 1)];

          // Two threads may be used to partition the reference arrays, analogous to
          // the use of two threads to merge the results for the merge sort algorithm.
          // Fill one reference array in ascending order with a child thread.
          std::future<void> partitionFuture =
            std::async(std::launch::async, [&] {
            for (size_t lower = start - 1, upper = median, j = start; j <= median; ++j) {
              KdNodePDB<K, V, M, N>* src_j = src[j];
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
            KdNodePDB<K, V, M, N>* src_k = src[k];
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

        // Delete the point array.
        delete[] point;

        // Recursively build the < branch of the tree with a child thread.
        // The recursive call to buildKdTree must be placed in a lambda
        // expression because buildKdTree is a template not a function.
        std::future<KdNodePDB<K, V, M, N>*> buildFuture =
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
     * The verifyKdTree function walks the k-d tree and checks that the
     * children of a node are in the correct branch of that node.
     *
     * Calling parameters:
     *
     * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
     * depth - the depth in the k-d tree
     *
     * returns: a count of the number of kdNodes in the k-d tree
     */
private:
  size_t verifyKdTree(signed_size_t maximumSubmitDepth, signed_size_t depth) {

    size_t count = 1;

    // The partition cycles as x, y, z, w...
    size_t p = depth % N;

    if (ltChild != nullptr) {
      if (ltChild->tuple[p] > tuple[p]) {
        std::cout << "child is > node!" << std::endl;
        exit(1);
      }
      if (superKeyCompare(ltChild->tuple, tuple, p) >= 0) {
        std::cout << "child is >= node!" << std::endl;
        exit(1);
      }
    }
    if (gtChild != nullptr) {
      if (gtChild->tuple[p] < tuple[p]) {
        std::cout << "child is < node!" << std::endl;
        exit(1);
      }
      if (superKeyCompare(gtChild->tuple, tuple, p) <= 0) {
        std::cout << "child is <= node" << std::endl;
        exit(1);
      }
    }

    // Verify the < branch with a child thread at as many levels of the tree as possible.
    // Create the child thread as high in the tree as possible for greater utilization.

    // Is a child thread available to build the < branch?
    if (maximumSubmitDepth < 0 || depth > maximumSubmitDepth) {

      // No, so verify the < branch with the current thread.
      if (ltChild != nullptr) {
        count += ltChild->verifyKdTree(maximumSubmitDepth, depth + 1);
      }

      // Then verify the > branch with the current thread.
      if (gtChild != nullptr) {
        count += gtChild->verifyKdTree(maximumSubmitDepth, depth + 1);
      }
    }
    else {

      // Yes, so verify the < branch with a child thread. Note that a
      // lambda is required to instantiate the verifyKdTree template.
      std::future<size_t> verifyFuture;
      if (ltChild != nullptr) {
        verifyFuture =
          std::async(std::launch::async, [&] {
          return ltChild->verifyKdTree(maximumSubmitDepth, depth + 1);
            });
      }

      // And simultaneously verify the > branch with the current thread.
      size_t gtCount = 0;
      if (gtChild != nullptr) {
        gtCount = gtChild->verifyKdTree(maximumSubmitDepth, depth + 1);
      }

      // Wait for the child thread to finish execution.
      size_t ltCount = 0;
      if (ltChild != nullptr) {
        try {
          ltCount = verifyFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }
      }
      count += ltCount + gtCount;
    }

    return count;
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
  static KdNodePDB<K, V, M, N>* createKdTree(std::vector<KdNodePDB<K, V, M, N>*>& kdNodes,
    size_t numThreads, signed_size_t maximumSubmitDepth) {

    // Create the references arrays including one additional array for use in building the k-d tree.
    KdNodePDB<K, V, M, N>*** references = new KdNodePDB<K, V, M, N>** [N + 1];

    // The first references array is the .data() array of the kdNodes vector.
    references[0] = kdNodes.data();

    // Allocate the remaining references arrays.
    for (size_t i = 1; i < N + 1; ++i) {
      references[i] = new KdNodePDB<K, V, M, N>* [kdNodes.size()];
    }

    // Sort the first reference array using multiple threads. Importantly,
    // for compatibility with the 'permutation' vector initialized below,
    // use the first dimension (0) as the leading key of the super key.
    int numthreads = 1;
    if (maximumSubmitDepth >= 0) numthreads = (1 << (maximumSubmitDepth + 1));
    parallelSort(references[0], references[0] + kdNodes.size(),
      [&](KdNodePDB<K, V, M, N>* a, KdNodePDB<K, V, M, N>* b) {
        return (0 > superKeyCompare(a->tuple, b->tuple, 0));
      }, numthreads);

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
    KdNodePDB<K, V, M, N>* root = buildKdTree(references, permutation, 0, end,
      maximumSubmitDepth, 0);
    // Verify the k-d tree and report the number of kdNodes.
    size_t numberOfNodes;
    numberOfNodes = root->verifyKdTree(maximumSubmitDepth, 0);
    // Delete all but the first of the references arrays.
    for (size_t i = 1; i < N + 1; ++i) {
      delete[] references[i];
    }
    delete[] references;

    // Return the pointer to the root of the k-d tree.
    return root;
  }

  /*
    * Walk the k-d tree to delete each KdNode.
    */
public:
  void deleteKdTree() {

    // Delete the < sub-tree.
    if (ltChild != nullptr) {
      ltChild->deleteKdTree();
    }
    // Delete the > sub-tree.
    if (gtChild != nullptr) {
      gtChild->deleteKdTree();
    }
    // Delete the current KdNode.
    delete this;
  }

  /*
    * The insideBounds function determines whether KdNode::tuple lies inside the
    * hyper-rectangle defined by the query lower and upper bound vectors.
    *
    * Calling parameters:
    *
    * queryLower - the query lower bound vector
    * queryUpper - the query upper bound vector
    * enable - a vector that specifies the dimensions on which to test for insidedness
    *
    * return true if inside, false if outside
    */
private:
  bool insideBounds(std::vector<K> const& queryLower, std::vector<K> const& queryUpper,
    std::vector<bool> const& enable) {
    bool inside = true;
    for (size_t i = 0; i < queryLower.size(); ++i) {
      if (enable[i] && (queryLower[i] > tuple[i] || queryUpper[i] < tuple[i])) {
        inside = false;
        break;
      }
    }
    return inside;
  }

    /*
     * The regionSearch function searches the k-d tree to find the KdNodes that
     * lie within a hyper-rectangle defined by the query lower and upper bounds.
     *
     * Calling parameters:
     *
     * queryLower - the query lower bound vector
     * queryUpper - the query upper bound vector
     * permutation - vector that specifies permutation of the partition coordinate
     * maximumSubmitDepth - the maximum tree depth at which a child task may be launched
     * depth - the depth in the k-d tree
     * enable - a vector that specifies the dimensions on which to prune the region search
     *
     * return a list that contains the KdNodes that lie within the cutoff distance of the query node
     */
#ifdef OLD_SEARCH
public:
  void regionSearch(std::list<retPair_t>& result,
    const std::vector<K>& queryLower, const std::vector<K>& queryUpper,
    const signed_size_t maximumSubmitDepth, const size_t depth, std::vector<int64_t>& permutation) {

    // Look up the primary coordinate index.
    unsigned int p = (unsigned int)permutation.at(depth);

    // the branchCode will be used later to select the actual branch configuration in the switch statement
    // below.  0 = no branch, 1 = < branch only, 2 = > branch only, 3 = both branches.
    int branchCode = 0;

    // Search the < branch of the k-d tree if the partition coordinate of the queryPlus is
    // <= the partition coordinate of the k-d node.  The < branch
    // must be searched when the cutoff distance equals the partition coordinate because the super
    // key may assign a point to either branch of the tree if the sorting or partition coordinate,
    // which forms the most significant portion of the super key, shows equality.
    if (queryLower[p] <= tuple[p]) {
      // but only search if the ltChild pointer is not null;
      if (ltChild != nullptr) branchCode = 1;
      // Search the > branch of the k-d tree if the partition coordinate of the queryPlus is
      // >= the partition coordinate of the k-d node.  The < branch
      // must be searched when the cutoff distance equals the partition coordinate because the super
      // key may assign a point to either branch of the tree if the sorting or partition coordinate,
      // which forms the most significant portion of the super key, shows equality.
      if (queryUpper[p] >= tuple[p]) {
        // but only if the gtChild pointer is not null;
        if (gtChild != nullptr) branchCode += 2;
        // while here check to see if the local tuple is inside the the hypercube.
        if (values != nullptr) {
          // If the distance from the query node to the k-d node is within the cutoff distance
          // in all k dimensions, add the k-d node to a list.
          bool inside = true;
          for (size_t i = 0; i < queryUpper.size(); i++) {
            if ((queryUpper[i] < tuple[i]) || (queryLower[i] >  tuple[i])) {
              inside = false;
              break;
            }
          }
          if (inside) {
            std::vector<K> tempTuple(tuple, tuple + queryLower.size());
            retPair_t tmpPair(tempTuple, *values);
            result.push_back(tmpPair);
          }
        }
      }
    }
    else { // will not decend the lt branch so lets check the gt.
      if (gtChild != nullptr && queryUpper[p] >= tuple[p]) branchCode = 2;
    }

    // now implenent the branching decided on earlier
    switch (branchCode) {
    case 0: // child pointer are both null so just return
      break;
    case 1: // only go down the less than branch
      ltChild->regionSearch(result, queryLower, queryUpper, maximumSubmitDepth, depth + 1, permutation);
      break;
    case 2: // only go down the greater than branch
      gtChild->regionSearch(result, queryLower, queryUpper, maximumSubmitDepth, depth + 1, permutation);
      break;
    case 3: // go down both branches
      if (depth <= maximumSubmitDepth) {
        // get a future and another list ready for a child thread
        std::future< void > searchFuture;
        std::list<retPair_t> ltResult;
        // Search the < branch asynchronously with a child thread.
        searchFuture = std::async(std::launch::async, [&] {
          return ltChild->regionSearch(ltResult, queryLower, queryUpper,
          maximumSubmitDepth, depth + 1, permutation); });
        // Search the > branch  with the master thread.
        gtChild->regionSearch(result, queryLower, queryUpper, maximumSubmitDepth, depth + 1, permutation);

        // Get the result of searching the < branch with the child thread.
        try {
          searchFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }
        result.splice(result.end(), ltResult);
      }
      else { // if below the maximum submit depth
        ltChild->regionSearch(result, queryLower, queryUpper, maximumSubmitDepth, depth + 1, permutation);
        gtChild->regionSearch(result, queryLower, queryUpper, maximumSubmitDepth, depth + 1, permutation);
      }
      break;
    }
    return;
  }

public:
  int64_t regionSearchAndRemove(std::list<retPair_t>& result,
    const std::vector<K>& queryLower, const std::vector<K>& queryUpper,
    const signed_size_t maximumSubmitDepth, const signed_size_t depth, std::vector<int64_t>& permutation) {

    // Look up the primary coordinate index.
    unsigned int p = (unsigned int)permutation[(unsigned int)depth];

    // the branchCode will be used later to select the actual branch configuration in the switch statement
    // below.  0 = no branch, 1 = < branch only, 2 = > branch only, 3 = both branches.
    int branchCode = 0;
    // the return codes indicate whether the status of the node that was just returned from
    // 0 = active node
    // 1 = the node and all nodes below are dead and it's can be removed from the tree.
    int64_t ltRetCode = 0;
    int64_t gtRetCode = 0;

    // Search the < branch of the k-d tree if the partition coordinate of the queryPlus is
    // <= the partition coordinate of the k-d node.  The < branch
    // must be searched when the cutoff distance equals the partition coordinate because the super
    // key may assign a point to either branch of the tree if the sorting or partition coordinate,
    // which forms the most significant portion of the super key, shows equality.
    if (queryLower[p] <= tuple[p]) {
      // but only search if the ltChild pointer is not null;
      if (ltChild != nullptr) branchCode = 1;
      // Search the > branch of the k-d tree if the partition coordinate of the queryPlus is
      // >= the partition coordinate of the k-d node.  The < branch
      // must be searched when the cutoff distance equals the partition coordinate because the super
      // key may assign a point to either branch of the tree if the sorting or partition coordinate,
      // which forms the most significant portion of the super key, shows equality.
      if (queryUpper[p] >= tuple[p]) {
        // but only if the gtChild pointer is not null;
        if (gtChild != nullptr) branchCode += 2;
        // while here check to see if the local tuple is inside the the hypercube.
        if (values != nullptr) {
          // If the distance from the query node to the k-d node is within the query rect
          // in all k dimensions, add the k-d node to a list.
          bool inside = true;
          for (size_t i = 0; i < queryUpper.size(); i++) {
            if ((queryUpper[i] < tuple[i]) || (queryLower[i] > tuple[i])) {
              inside = false;
              break;
            }
          }
          if (inside) {
            std::vector<K> tempTuple(tuple, tuple + queryLower.size());
            retPair_t tmpPair(tempTuple, *values);
            result.push_back(tmpPair);
            delete values;
            values = nullptr;   // mark the node dead by nulling the pointer
          }
        }
      }
    }
    else { // will not decend the lt branch so lets check the gt.
      if (gtChild != nullptr && queryUpper[p] >= tuple[p]) branchCode = 2;
    }

    // now implenent the branching decided on earlier
    switch (branchCode) {
    case 0: // child pointer are both null so just return
      break;
    case 1: // only go down the less than branch
      ltRetCode = ltChild->regionSearchAndRemove(result, queryLower, queryUpper, maximumSubmitDepth, depth + 1, permutation);
      break;
    case 2: // only go down the greater than branch
      gtRetCode = gtChild->regionSearchAndRemove(result, queryLower, queryUpper, maximumSubmitDepth, depth + 1, permutation);
      break;
    case 3: // go down both branches
      if (depth <= maximumSubmitDepth) {
        // get a future and another list ready for a child thread
        std::future< int64_t > searchFuture;
        std::list<retPair_t> ltResult;
        // Search the < branch asynchronously with a child thread.
        searchFuture = std::async(std::launch::async, [&] {
          return ltChild->regionSearchAndRemove(ltResult, queryLower, queryUpper,
          maximumSubmitDepth, depth + 1, permutation); });
        // Search the > branch  with the master thread.
        gtRetCode = gtChild->regionSearchAndRemove(result, queryLower, queryUpper, maximumSubmitDepth, depth + 1, permutation);

        // Get the result of searching the < branch with the child thread.
        try {
          ltRetCode = searchFuture.get();
        }
        catch (std::exception const& e) {
          std::cout << "caught exception " << e.what() << std::endl;
        }
        result.splice(result.end(), ltResult);
      }
      else { // if below the maximum submit depth
        ltRetCode = ltChild->regionSearchAndRemove(result, queryLower, queryUpper, maximumSubmitDepth, depth + 1, permutation);
        gtRetCode = gtChild->regionSearchAndRemove(result, queryLower, queryUpper, maximumSubmitDepth, depth + 1, permutation);
      }
      break;
    }
    if (ltRetCode == 1) {
      delete ltChild;
      ltChild = nullptr;
    }
    if (gtRetCode == 1) {
      delete gtChild;
      gtChild = nullptr;
    }
    // retun a dead node code if all the pointers in this node are null;
    if (values == nullptr && ltChild == nullptr && gtChild == nullptr) {
      return 1;
    }
    return 0;
  }
#endif
public:
  SearchRet regionSearchAndRemove(std::list<retPair_t*>& result,
    const std::vector<K>& queryLower, const std::vector<K>& queryUpper,
    M* const moveTo, const signed_size_t depth, const std::vector<int64_t>& permutation) {
    
    // get a local copy of fthe cluster ID;
    const uint32_t& moveToID = moveTo->getID();

    // Look up the primary coordinate index.
    unsigned int p = (unsigned int)permutation[(unsigned int)depth];

    // the branchCode will be used later to select the actual branch configuration in the switch statement
    // below.  0 = no branch, 1 = < branch only, 2 = > branch only, 3 = both branches.
    int branchCode = 0;
    // the return codes indicate whether the status of the node that was just returned from
    // init the return result to NoResut
    SearchRet ltRetCode = NoResult;
    SearchRet gtRetCode = NoResult;
    SearchRet returnResult = AllTaken;

    // Search the < branch of the k-d tree if the partition coordinate of the queryPlus is
    // <= the partition coordinate of the k-d node.  The < branch
    // must be searched when the cutoff distance equals the partition coordinate because the super
    // key may assign a point to either branch of the tree if the sorting or partition coordinate,
    // which forms the most significant portion of the super key, shows equality.
    if (queryLower[p] <= tuple[p]) {
      // but only search if the ltChild pointer is not null and all nodes below are not taken for this container
      if (ltChild != nullptr && (ltAllTaken != moveToID)) branchCode = 1;
      // Search the > branch of the k-d tree if the partition coordinate of the queryPlus is
      // >= the partition coordinate of the k-d node.  The < branch
      // must be searched when the cutoff distance equals the partition coordinate because the super
      // key may assign a point to either branch of the tree if the sorting or partition coordinate,
      // which forms the most significant portion of the super key, shows equality.
      if (queryUpper[p] >= tuple[p]) {
        // but only search if the gtChild pointer is not null and all nodes below are not taken for this container
        if (gtChild != nullptr && (gtAllTaken != moveToID)) branchCode += 2;
        // while here check to see if the local tuple is inside the the hypercube.
        // If the distance from the query node to the k-d node is within the query rect
        // in all k dimensions, add the k-d node to a list.
        bool inside = true;
        for (size_t i = 0; i < queryUpper.size(); i++) {
          if ((queryUpper[i] < tuple[i]) || (queryLower[i] > tuple[i])) {
            inside = false;
            break;
          }
        }
        if (inside) {
          // grab a lock on this node
          std::lock_guard<std::mutex> guard(thisMutex);
          // check if the node had not already been taken 
          if (values != nullptr) {
            // if it has not, add this nodes data to the return list and mark the node as taken
            retPair_t* tmpPair = new retPair_t(std::vector<K>(tuple, tuple + N), values);
            result.push_back(tmpPair);
            movedTo = moveTo;  // mark the value as taken by another container.
            values = nullptr;   // mark the node dead by nulling the pointer
          } else if (moveTo != movedTo) {
            // else indicate this cluster overlaps another cluster
            if (moveTo->overlaps == nullptr) {
              moveTo->overlaps = new std::set<M*>;
            }
            
            auto rslt = moveTo->overlaps->insert(movedTo);
#ifdef DEBUG_PRINT
            static std::mutex printMutex;
            if (rslt.second == true) {
              std::lock_guard<std::mutex> guard(printMutex);
              std::cout << "cluster " << moveToID << " overlaps with cluster " << movedTo->getID() << "\n";
            }
#endif
          }
        }
        returnResult = ResultFoundAndAllTaken;
        
      }
    }
    else { // will not decend the lt branch so lets check the gt.
      if (queryUpper[p] >= tuple[p]) {
        // but only search if the gtChild pointer is not null and all nodes below are not taken for this container
        if (gtChild != nullptr && (gtAllTaken != moveToID)) branchCode = 2;
      }
    }

    // now implenent the branching decided on earlier
    switch (branchCode) {
    case 0: // child pointer are both null so just return
      break;
    case 1: // only go down the less than branch
      ltRetCode = ltChild->regionSearchAndRemove(result, queryLower, queryUpper, moveTo, depth + 1, permutation);
      if (ltRetCode == AllTaken || ltRetCode == ResultFoundAndAllTaken) {
        ltAllTaken = moveToID;
      }
      break;
    case 2: // only go down the greater than branch
      gtRetCode = gtChild->regionSearchAndRemove(result, queryLower, queryUpper, moveTo, depth + 1, permutation);
      if (gtRetCode == AllTaken || gtRetCode == ResultFoundAndAllTaken) {
        gtAllTaken = moveToID;
      }
      break;
    case 3: // go down both branches
      ltRetCode = ltChild->regionSearchAndRemove(result, queryLower, queryUpper, moveTo, depth + 1, permutation);
      if (ltRetCode == AllTaken || ltRetCode == ResultFoundAndAllTaken) {
        ltAllTaken = moveToID;
      }
      gtRetCode = gtChild->regionSearchAndRemove(result, queryLower, queryUpper, moveTo, depth + 1, permutation);
      if (gtRetCode == AllTaken || gtRetCode == ResultFoundAndAllTaken) {
        gtAllTaken = moveToID;
      }
      break;
    }

    // Now figure out what to return.  If something was found either here or below return 1 or
    // if this node is still active as indicated by non null pointer set returnResult to 1.
    //Otherwise return a 0.
    if (returnResult == ResultFoundAndAllTaken && (movedTo != moveTo || !(ltAllTaken == moveToID || ltChild == nullptr) || !(gtAllTaken == moveToID || gtChild == nullptr))) {
      returnResult = ResultFound;
    }
    else if (returnResult == AllTaken && (movedTo != moveTo || !(ltAllTaken == moveToID || ltChild == nullptr) || !(gtAllTaken == moveToID || gtChild == nullptr))) {
      returnResult = NoResult;
    }
    return returnResult;
  }


  //TODO  Add search functons with the enable vector

  /**
  The pickValue picks a value from the kdTree by decending the tree until it finds
  a knNode where both child pointers are null and returns that value
  Paramaeters:
    returnKV:   refernece to a retPair where to place the value.
    selector:   unsigned in where the bits are used to guide which clikd node to decend
    removePick: bool where if true will cause the node found to be removed
    depth:      integer counting the number of levels decended.

    return:     integer indicating the state of the node

  **/
#ifdef OLD_SEARCH
public:
  SearchRet pickValue(retPair_t& returnKV,
    const uint64_t selector, bool removePick, const int depth) {

    // init the return result to NoResut
    SearchRet returnResult = NoResult;

    bool goGtThan = (selector & 0x1) == 1;

    // if greater than or less than, continue the search on  the appropriate child node
    if ((!goGtThan || gtChild == nullptr) && ltChild != nullptr) {
      returnResult = ltChild->pickValue(returnKV, selector >> 1, removePick, depth + 1);
      // if the node below declared itself dead, remove the link
      if (removePick && returnResult == DeadNode) {
        delete ltChild;
        ltChild = nullptr;
      }
    }
    else if ((goGtThan || ltChild == nullptr) && gtChild != nullptr) {
      returnResult = gtChild->pickValue(returnKV, selector >> 1, removePick, depth + 1);
      // if the node below declared itself dead, remove the link
      if (removePick && returnResult == DeadNode) {
        delete gtChild;
        gtChild = nullptr;
      }
    }
    else {// both child pointers must be null to get here so this is a leaf node
      if (values != nullptr) {
        std::vector<K> key(N);
        for (size_t i = 0; i < N; i++) {
          key[i] = tuple[i];
        }
        returnKV = retPair_t(key, *values);
        if (removePick) {
          delete values;
          values = nullptr;
          returnResult = DeadNode;  //flag for possible node removal
        }
        else {
          returnResult = ResultFound;
        }
      }
    }
    // Now figure out what to return.  If something was found either here or below return 1 or
    // if this node is still active as indicated by non null pointer set returnResult to 1.
    //Otherwise return a 0.
    if (returnResult == DeadNode && (values != nullptr || ltChild != nullptr || gtChild != nullptr)) {
      returnResult = ResultFound;
    }
    return returnResult;
  }


  /**
  The pickValue picks a value from the kdTree by decending the tree until it finds
  a knNode where both child pointers are null and returns that value
  Paramaeters:
    returnKV:   refernece to a retPair where to place the value.
    selector:   unsigned in where the bits are used to guide which clikd node to decend
    removePick: bool where if true will cause the value to be removed
    moveTo:     mointer to some object that will hold the value found
    depth:      integer counting the number of levels decended.

    return:     SearchRet indicating the state of the node

  **/
#endif
public:
  SearchRet pickValue(retPair_t& returnKV,
    const uint64_t selector, M* const moveTo, const bool removePick, const int depth) {

    const uint32_t& moveToID = moveTo->getID();
      
    // init the return result to 0
    SearchRet returnResult = AllTaken;

    bool goGtThan = (selector & 0x1) == 1;

    // if greater than or less than, continue the search on  the appropriate child node
    if ((!goGtThan || gtAllTaken != 0 || gtChild == nullptr ) && (ltAllTaken == 0 && ltChild != nullptr)) {
      returnResult = ltChild->pickValue(returnKV, selector >> 1, moveTo, removePick, depth + 1);
      // if the node below declared itself allTaken, set the gtAllTaken flag
      if (removePick && (returnResult == AllTaken || returnResult == ResultFoundAndAllTaken)) {  // i
        ltAllTaken = moveToID;
      }
      if ((returnResult == AllTaken || returnResult == NoResult) && gtChild != nullptr) {
        returnResult = gtChild->pickValue(returnKV, selector >> 1, moveTo, removePick, depth + 1);
        // if the node below declared itself allTaken, set the gtAllTaken flag
        if (removePick && (returnResult == AllTaken || returnResult == ResultFoundAndAllTaken)) {  // i
          gtAllTaken = moveToID;
        }
      }
    }
    else if ((goGtThan || ltAllTaken != 0 || ltChild == nullptr) && (gtAllTaken == 0 && gtChild != nullptr)) {
      returnResult = gtChild->pickValue(returnKV, selector >> 1, moveTo, removePick, depth + 1);
      // if the node below declared itself allTaken, set the gtAllTaken flag
      if (removePick && (returnResult == AllTaken || returnResult == ResultFoundAndAllTaken)) {  // i
        gtAllTaken = moveToID;
      }
      if ((returnResult == AllTaken || returnResult == NoResult) && ltChild != nullptr) {
        returnResult = ltChild->pickValue(returnKV, selector >> 1, moveTo, removePick, depth + 1);
        // if the node below declared itself allTaken, set the gtAllTaken flag
        if (removePick && (returnResult == AllTaken || returnResult == ResultFoundAndAllTaken)) {  // i
          ltAllTaken = moveToID;
        }
      }
    }
    if (returnResult == NoResult || returnResult == AllTaken) {
      if (values != nullptr) {
        // grab a lock
        std::lock_guard<std::mutex> guard(thisMutex);
        // check again to make sure this node was not taken
        if (values != nullptr) {
          returnKV = retPair_t(std::vector<K>(tuple, tuple + N), values);
          movedTo = moveTo;  // mark the node as taken by the searching container.
          if (removePick) {
            values = nullptr;
            returnResult = ResultFoundAndAllTaken;  //flag for possible node all taken
          } else {
            returnResult = ResultFound;
          }
        }  else {
          returnResult = AllTaken;  // possibly all taken, check later
        }
      } else {
        returnResult = AllTaken;  // possibly all taken, check later
      }
    }
    // Now figure out what to return.  
    // If the returnResult indcates ResultFoundAllTaken from a lower call but this not is not all taken, change returnResult to only ResultFfund
     if (returnResult == ResultFoundAndAllTaken && (movedTo != moveTo || !(ltAllTaken == moveToID || ltChild == nullptr) || !(gtAllTaken == moveToID || gtChild == nullptr))) {
      returnResult = ResultFound;
    } else if (returnResult == AllTaken && (movedTo != moveTo || !(ltAllTaken == moveToID || ltChild == nullptr) || !(gtAllTaken == moveToID || gtChild == nullptr))) {
      returnResult = NoResult;
    }
    return returnResult;
  }



 }; // class KdNodePDB

  
  

  template<typename K, typename V, typename M, size_t N>
class KdTreePDB {

  typedef std::pair<std::vector<K>, std::list<V>*> retPair_t;

private:
  size_t numPoints = 0;  // total number of points submitted to the KdTree.
  std::vector<KdNodePDB<K, V, M, N>*> kdNodes;  // vector of kdNodes reslting from submissions to the KdTree.
  KdNodePDB<K, V, M, N>* root = nullptr;  // root of the KdTree after it's built
  std::vector<int64_t> permutation; // vector of precalculated levels
  signed_size_t numThreads = 1;     // number of threads used in the process
  signed_size_t maximumSubmitDepth = -1; // maximum depth in tree building or searching where new threads are submitted.

  // calculats the maximumSubmitDepth given the number of threads.
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
  KdTreePDB<K, V, M, N>() {
    numThreads = std::thread::hardware_concurrency();
    calcMaximumSubmitDepth();
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
    auto kp = new KdNodePDB<K, V, M, N>(tuple, value);
    kdNodes.push_back(kp);
    numPoints = kdNodes.size();
    return numPoints;
  }

 /**
 * <p>
 * The {@code searchKdTree} method searches the k-d tree and finds the KdNodes
 * that lie within a cutoff distance from a query node in all k dimensions.
 * </p>
 *
 * @param result - List of tuple,value pairs that are within the search region.
 * @param queryLower - Array containing the lager search bound for each dimension
 * @param queryUpper - Array containing the smaller search bound for each dimension
 * @param numThreads - the maximum tree depth at which a thread may be launched
 * @return bool indicating success.
 */

  bool searchRegion(std::list<retPair_t>& retPair, std::vector<K>& queryLower, std::vector<K>& queryUpper) {
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
    root->regionSearch(retPair, queryLower, queryUpper, maximumSubmitDepth, 0, permutation);
    return true;
  }

  /**
  * <p>
  * The {@code searchKdTree} method searches the k-d tree and finds the KdNodes
  * that lie within a cutoff distance from a query node in all k dimensions.
  * </p>
  *
  * @param result - values where the associated tuples that are within the search region.
  * @param queryLower - Array containing the lager search bound for each dimension
  * @param queryUpper - Array containing the smaller search bound for each dimension
  * @param numThreads - the maximum tree depth at which a thread may be launched
  * @return bool indicating success.
  */

  bool searchRegion(std::list<V>& retVal, std::vector<K>& queryLower, std::vector<K>& queryUpper) {
      // if the tree is not built yet, build it
      if (root == nullptr) {
          buildTree();
      }
      std::list<retPair_t> retPair;
      bool retFlag = searchRegion(retPair, queryLower, queryUpper, 1);
      for (auto pp : retPair) {
        retVal.splice(retVal.begin(), *pp.second);
      }
      return retFlag;
  }

  /*
  Same as searchRegion but the nodes that contained the items added to the retPair list are removed from the KdTree
  */
  bool searchRegionAndRemove(std::list<retPair_t>& retPair, std::vector<K>& queryLower, std::vector<K>& queryUpper) {
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
    auto retCode = root->regionSearchAndRemove(retPair, queryLower, queryUpper, maximumSubmitDepth, 0l, permutation);
    return true;
  }

  bool searchRegionAndRemove(std::list<retPair_t*>& retPair, std::vector<K>& queryLower, std::vector<K>& queryUpper, M* moveTo) {
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
    root->regionSearchAndRemove(retPair, queryLower, queryUpper, moveTo, 0l, permutation);
    return true;
  }



  bool buildTree(signed_size_t numThreads = -1) {
    if (numThreads != -1) {
      setNumThreads(numThreads);
    }
    permutation = getPermutations(numPoints, N);
    root = KdNodePDB<K, V, M, N>::createKdTree(kdNodes, numThreads, maximumSubmitDepth);
    return true;
  }

  
  bool pickValue(retPair_t& returnKV, uint64_t selectionBias, M* moveTo, bool remove) {
      // if the tree is not built yet, build it
    if (root == nullptr) {
        buildTree();
        // if root is still null; return a null
        if (root == nullptr) return false;
      }
      // descent selector
    // search the tree to find the node to return and possibly delete.
    SearchRet returnResult = root->pickValue(returnKV, selectionBias, moveTo, remove, 0);
    // Check to see if no data was returned.  If so then return false.
    if (returnResult == NoResult || returnResult == AllTaken) return false;

    return true;
  }


  ~KdTreePDB() {
    root->deleteKdTree();
  }

}; // KdTree

#endif // KDTREE_HPP
