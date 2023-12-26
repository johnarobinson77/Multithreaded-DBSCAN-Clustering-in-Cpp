
/**
* parallelSot.hpp
* 
 * Copyright (c) 2023 John Robinson.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

#ifndef PARALLELSORT_HPP
#define PARALLELSORT_HPP

#include <cstring>
#include <stdint.h>
#include <algorithm>    // std::swap
#include <cmath>        // log2
#include "parallelFor.hpp"

#define sortFor false
#define sortRev true


// the next series of methods are utility functions that are used by the various sort routines.

#define minimum(a,b) ((a)<(b) ? (a) : (b))
#define maximum(a,b) ((a)>(b) ? (a) : (b))


  // this routine does an integer division rounded up.
  inline size_t iDivUp(size_t a, size_t b)
  {
    return ((a % b) == 0) ? (a / b) : (a / b + 1);
  }

  // mergeFF merges two sorted ranges in the src array into a single sorted range in the dst array
// aBeg and aEnd inclusive indicate one sorted range 
// bBeg and bEnd inclusive indicate the other sorted range
// dBeg indicates where in the dst array the merge list should start
// First, the elements at aEnd and bEnd are compared.  The end index with the smaller value is used.
// to cap the merge function.  After than, the reset of the other array is just copied to the dst array
// The two ranges being merged do not have to be adjacent in memory.
  template< class RandomItD, class RandomItS, class CF>
  inline void mergeFF(RandomItD dst, RandomItS src, size_t aBeg, size_t aEnd, size_t bBeg, size_t bEnd, size_t dBeg, CF compFunc) {
    if (compFunc(*(src+aEnd), *(src+bEnd))) { // determine which range will be completed first during a compare and copy loop
      while (aBeg <= aEnd) {  // the a range will be completely copied first so only compare up the end of a
        *(dst+dBeg++) = !compFunc(*(src+bBeg), *(src+aBeg)) ? *(src+aBeg++) : *(src+bBeg++);
      }
      while (bBeg <= bEnd) { // then copy the rest of b
        *(dst+dBeg++) = *(src+bBeg++);
      }
    }
    else {
      while (bBeg <= bEnd) {  // the b range will be completely copied first so only compare up the end of b
        *(dst+dBeg++) = compFunc(*(src+aBeg), *(src+bBeg)) ? *(src+aBeg++) : *(src+bBeg++);
      }
      while (aBeg <= aEnd) { // then copy the rest of a
        *(dst+dBeg++) = *(src+aBeg++);
      }
    }
  }


  // ParallelMerge and it's supporting functions getMergePaths and mergePath are a CPU implementation of
  // parallel merge function developed for GPUs described in "GPU Merge Path: A GPU Merging Algorithm" by Greenand, McColl, and Bader.
  // Proceedings of the 26th ACM International Conference on Supercomputing

  // For a particular output element of a merge, find all the elements in valA and ValB that will be below it.  
  template <class RandomIt, class CF>
  size_t mergePath(RandomIt valA, int64_t aCount, RandomIt valB, int64_t bCount, int64_t diag, CF compFunc, size_t threads) {

    size_t begin = maximum(0, diag - bCount);
    size_t end = minimum(diag, aCount);

    while (begin < end) {
      size_t mid = begin + ((end - begin) >> 1);
      //T  aVal = valA[mid];
      //T  bVal = valB[diag - 1 - mid];
      bool pred = compFunc( *(valA + mid), *(valB + (diag - 1 - mid)) );
      if (pred) begin = mid + 1;
      else end = mid;
    }
    return begin;
  }

  // Divide the output of the merge into #-of-threads equal segments and determine the elements of valA and valB below each segment boundary. 
  template <class RandomIt, class CF>
  inline void getMergePaths(size_t mpi[], RandomIt valA, int64_t aCount, RandomIt valB, int64_t bCount, int64_t spacing, CF compFunc, size_t threads) {

    int64_t partitions = threads;
    mpi[0] = 0;
    mpi[partitions] = aCount;

    for (int i = 1; i < partitions; i++) {
      //    parallelFor(1, partitions, [=](int64_t i, int64_t nothing) {
      int64_t diag = i * spacing;
      mpi[i] = mergePath(valA, aCount, valB, bCount, diag, compFunc, threads);
      //      }, threads);
    }
  }

  // parallelMerge() merges two sorted ranges in the src array into a single sorted range in the dst array
  // aBeg and aEnd inclusive indicate one sorted range in the src array 
  // bBeg and bEnd inclusive indicate the other sorted range in the src array
  // dBeg indicates where in the dst array the merge list should start
  // After determining the ranges of the src segments that will go nto each output segment, use the mergeFF() function to 
  // merge those segments in parallel.
  template< class RandomItS, class RandomItD, class CF>
  inline void parallelMerge(RandomItD dst, RandomItS src, size_t aBeg, size_t aEnd, size_t bBeg, size_t bEnd, size_t dBeg, CF compFunc, size_t threads) {

    size_t mpi[1024];  // an array that hold the intermediate ranges during a parallel merge.

    int64_t aCount = aEnd - aBeg + 1;
    int64_t bCount = bEnd - bBeg + 1;
    int64_t spacing = iDivUp(aCount + bCount, threads);
    getMergePaths(mpi, src+aBeg, aCount, src+bBeg, bCount, spacing, compFunc, threads);

    parallelFor((int64_t)0, (int64_t)threads, [=](int64_t thread) {
      size_t grid = thread * spacing;		// Calculate the relevant index ranges into the source array
      size_t a0 = mpi[thread] + aBeg;				// for this partition
      size_t a1 = mpi[thread + 1] + aBeg;
      size_t b0 = (grid - mpi[thread]) + bBeg;
      size_t b1 = (minimum(aCount + bCount, (int64_t)grid + spacing) - mpi[thread + 1]) + bBeg;
      size_t wtid = dBeg + thread * spacing;  //Place where this thread will start writing the data

      if (a0 == a1) {							// If no a data just copy b
        for (size_t b = b0; b < b1; b++) dst[wtid++] = src[b];
      }
      else if (b0 == b1) {					// If no b data just copy a
        for (size_t a = a0; a < a1; a++) dst[wtid++] = src[a];
      }
      else { // else do a merge using the forward-forward merge
        mergeFF(dst, src, a0, a1 - 1, b0, b1 - 1, wtid, compFunc);
      }
      }, threads);
  }

  // this function is only for debug purposes.
  template<typename T> // 
  void print(std::string t, T* in, int n) {
    std::cout << t << " ";
    for (size_t i = 0; i < n; i++) std::cout << in[i] << " ";
    std::cout << std::endl;
  }

//#define BALANCED_MULTITHREADING
#ifdef BALANCED_MULTITHREADING

#pragma message ("Compiling  BALANCED_MULTITHREADING mode")

  template< class RandomIt, class CF>
  void parallelSort(RandomIt begin, RandomIt end, CF compFunc , size_t threads = 0) {
    // default number of threads iw the hardware number of cores.
    if (threads == 0) threads = std::thread::hardware_concurrency();
    // Get the total size;
    const size_t len = end - begin;

    // The current code fails if the number of elements is small relative to the number of threads.
    // So limit the number of threads to 1 for small sort cases.
    if (len < 100) threads = 1;

    // calculate the fractional size of each segment to sort.
    // calculating the arry segments using doubles results in segment sizes where the max segment size
    // is only one bigger than the min.
    double delta = double(len) / double(threads);

    //sort threads segments of the input arry using the sort method provided in the function pointer
    parallelFor((int64_t)0, (int64_t)threads, [begin, delta, compFunc](int64_t i) {
      int64_t lb = llround(i * delta);
      int64_t le = llround((i + 1) * delta);
      std::sort(begin + lb, begin + le, compFunc);
      }, threads);

    if (threads <= 1) return;

    // create a local swap space

    typedef typename std::iterator_traits<RandomIt>::value_type T;

     T *swap;  // pointer to array that that the data will be swapped to during a merge function
    swap = new typename std::iterator_traits<RandomIt>::value_type[len];

    const int64_t depth = (int64_t)ceil(log2(threads)); // calculate the number of depth iterations

    double next;
    double start;
    for (int64_t d = depth; d > 0; d-=2) {
      start = 0.0;
      // merged the full sized pairs of of segments for this level.
      for (next = 2.0 * delta; llround(next) <= (int64_t)len; next += 2.0 * delta) {
        parallelMerge(swap, begin, llround(start), llround(start + delta) - 1, llround(start + delta), llround(next) - 1, llround(start), compFunc, threads);
        start = next;
      }
      // if there is an odd number segments to be merged, take care of the left overs.
      if (llround(next) > (int64_t)len) {
        int64_t mid = llround(start + delta);
        if (mid > (int64_t)len) mid = len;   // do a merge if there is a segment plus a partial segment
        parallelMerge(swap, begin, llround(start), mid - 1, mid, len - 1, llround(start), compFunc, threads);
      }
      delta *= 2.0;

      start = 0.0;
      // merged the full sized pairs of of segments for this level.
      for (next = 2.0 * delta; llround(next) <= (int64_t)len; next += 2.0 * delta) {
        parallelMerge(begin, swap, llround(start), llround(start + delta) - 1, llround(start + delta), llround(next) - 1, llround(start), compFunc, threads);
        start = next;
      }
      // if there is an odd number segments to be merged, take care of the left overs.
      if (llround(next) > (int64_t)len) {
        int64_t mid = llround(start + delta);
        if (mid > (int64_t)len) mid = (int64_t)len;   // do a merge if there is a segment plus a partial segment
        parallelMerge(begin, swap, llround(start), mid - 1, mid, len - 1, llround(start), compFunc, threads);
      }
      delta *= 2.0;

    }
    // clean up
    delete[] swap;

  }

#else 

#pragma message ("Compiling MINIMIZED_THREAD_LAUNCH mode")

  template< class RandomIt, class CF>
  void parallelSort(RandomIt begin, RandomIt end, CF compFunc, size_t threads = 0) {
    // default number of threads iw the hardware number of cores.
    if (threads == 0) threads = std::thread::hardware_concurrency();
    // Get the total size;
    const size_t len = end - begin;

    // The current code fails if the number of elements is small relative to the number of threads.
    // So limit the number of threads to 1 for small sort cases.
    if (len < 100) threads = 1;

    // calculate the fractional size of each segment to sort.
    // calculating the arry segments using doubles results in segment sizes where the max segment size
    // is only one bigger than the min.
    double delta = double(len) / double(threads);

    //sort threads segments of the input arry using the sort method provided in the function pointer
    parallelFor((int64_t)0, (int64_t)threads, [begin, delta, compFunc](int64_t i) {
      int64_t lb = llround(i * delta);
      int64_t le = llround((i + 1) * delta);
      std::sort(begin + lb, begin + le, compFunc);
      }, threads);

    if (threads <= 1) return;

    // create a local swap space

    typedef typename std::iterator_traits<RandomIt>::value_type T;

    T* swap;  // pointer to array that that the data will be swapped to during a merge function
    swap = new typename std::iterator_traits<RandomIt>::value_type[len];

    const int64_t depth = (int64_t)ceil(log2(threads)); // calculate the number of depth iterations

    int64_t loops; // number of full size (2*delta) segment merges that need to be run
    int64_t tasks; // number of parallel tasks tht need to be run (loops + leftover is any)
    int64_t pmThreads; // number of threads to apply to full size parallelMerges (inside the parallel loop)
    int64_t loThreads; // number of threads to apply to leftover parallel_merge.
    double loStart;  // the starting point if the left overs
    for (int64_t d = depth; d > 0; d -= 2) {
      // This section computes teh number of tasks that need to be performed and how many threads to 
      // give to each task.  After the above sort, there are <threads> number of segments that need 
      // to be iteratively merged.  Each segment is initially delta large.  The first set of tasks is 
      // to merge two adjacent segments that are exactly delta each.  The variable loops the number
      // of these that need to be processed. is  However there are cases where 
      // this does not cover all the segments, i.e. the case of an odd number of segments.  In such cases,
      // the "left over" segments need to be handled outside the main loop and adds one more task to be handled
      // pmThreads is the number of threads given to each parallelMerge in the full delta sized segments
      // loThreads is the number of threads assigned to the left over parallelMerge and will be 0
      // if there is no left over or = to or 1 less than pmThreads.
      loops = (int64_t)floor(double(len) / (2.0 * delta));
      loStart = double(loops) * 2.0 * delta;
      tasks = loops + ((llround(loStart) < (int64_t)len) ? 1 : 0);
      pmThreads = iDivUp(threads, tasks);
      loThreads = threads - (loops * pmThreads);
      // merged the full sized pairs of of segments for this level.
      auto futures = parallelForNoWait((int64_t)0, loops, [=](int64_t i) {
        int64_t lb = llround(double(2 * i) * delta);
        int64_t lm = llround(double(2 * i + 1) * delta);
        int64_t le = llround(double(2 * i + 2) * delta);
        parallelMerge(swap, begin, lb, lm - 1, lm, le - 1, lb, compFunc, pmThreads);
        }, loops);
      // if there is an odd number segments to be merged, take care of the left overs.
      if (llround(loStart) < (int64_t)len) {
        int64_t mid = llround(loStart + delta);
        if (mid > (int64_t)len) mid = len;   // do a merge if there is a segment plus a partial segment
        parallelMerge(swap, begin, llround(loStart), mid - 1, mid, len - 1, llround(loStart), compFunc, loThreads);
      }
      parallelForFinish(futures);
      delta *= 2.0;

      loops = (int64_t)floor(double(len) / (2.0 * delta));
      loStart = double(loops) * 2.0 * delta;
      tasks = loops + ((llround(loStart) < (int64_t)len) ? 1 : 0);
      pmThreads = iDivUp(threads, tasks);
      loThreads = threads - (loops * pmThreads);
      // merged the full sized pairs of of segments for this level.
      //for (int64_t i = 0; i < loops; i++) {
      futures = parallelForNoWait((int64_t)0, loops, [=](int64_t i) {
        int64_t lb = llround(double(2 * i) * delta);
        int64_t lm = llround(double(2 * i + 1) * delta);
        int64_t le = llround(double(2 * i + 2) * delta);
        parallelMerge(begin, swap, lb, lm - 1, lm, le - 1, lb, compFunc, pmThreads);
        }, loops);
      // if there is an odd number segments to be merged, take care of the left overs.
      if (llround(loStart) < (int64_t)len) {
        int64_t mid = llround(loStart + delta);
        if (mid > (int64_t)len) mid = len;   // do a merge if there is a segment plus a partial segment
        parallelMerge(begin, swap, llround(loStart), mid - 1, mid, len - 1, llround(loStart), compFunc, loThreads);
      }
      parallelForFinish(futures);
      delta *= 2.0;  // double the size of delta for the level of merges
    }
    // clean up
    delete[] swap;

  }
#endif

  template< class RandomIt>
  void parallelSort(RandomIt begin, RandomIt end, size_t threads = 0) {
    parallelSort(begin, end, std::less<typename std::iterator_traits<RandomIt>::value_type>(), threads);
  }

#endif // PARALLELSORT_HPP