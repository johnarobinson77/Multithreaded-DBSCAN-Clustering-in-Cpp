/**
 * Copyright (c) 2023 John Robinson.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

#ifndef THIS_PARALLEL_FOR
#define THIS_PARALLEL_FOR

#include <iostream>
#include <future>
#include <mutex>
#include <vector>
#include <numeric>

//#define DEBUG_PRINT

// parallel_for runs the equivalent to the statement "for (int i = begin; i < end; i++) fn(i);"
// where fn is a lambda.  if num_segs is greater than 1, it divides the the begin-end number range 
// into that many  segments and start starts separate threads to execute the first num_segs-1 segments
// The last segment is run on the current thread.  
// parallel_for returns after all treads segments have completed.
template< class RandomIt, class FN >
void parallelFor(const RandomIt begin, const RandomIt end, FN fn, int64_t num_segs = 1) {

  // compute the number iterations and if 0, then return
  const int64_t n = end - begin;
  if (n == 0) return;
  // if the number of elements is less than the number of segments or threads, reduce the number of segments
  if (n < num_segs) num_segs = n;
  // create the array of futures.
  std::vector<std::future<void>> futures;
  // Compute the number of iterations to preform on each segment using a double to preserve the fraction.
  // using the method of iterating through bounds applied to each thread using a double and then rounding to
  // the nearest integer ensures that the difference in the number of iterations that each thread will do is 
  // no greater than one.
  double seg_size = (double)n / (double)num_segs;
  int64_t sBeg = 0;
  double sEnd = seg_size;
  // start threads for one less than the number of segments
  for (int64_t seg = 0; seg < num_segs-1; seg++) { 
    futures.push_back(std::async(std::launch::async, [begin, fn](int64_t lb, int64_t le) {
#ifdef DEBUG_PRINT
      {
        static std::mutex lock;
        std::lock_guard<std::mutex> guard(lock);
        std::cout << "parallel_for: n == " << le - lb << " thread: " << std::this_thread::get_id() << '\n';
      }
#endif // DEBUG_PRINT
      for (int64_t i = lb; i < le; i++) fn(begin + i);
      }, sBeg, llround(sEnd)));
    sBeg = llround(sEnd);
    sEnd += seg_size;
  }
  // execute the last segment in this thread
  for (int64_t i = sBeg; i < llround(sEnd); i++) fn(begin + i);

  // wait for the segments to complete
  for (int seg = 0; seg < futures.size(); seg++)
    futures[seg].wait();
}

// parallelFor runs the equivalent to the statement "for (int i = begin; i < end; i++) fn(i);"
// where fn is a lambda.  if num_segs is greater than 1, it divides the the begin-end number range 
// into that many  segments and start starts separate threads to execute the first num_segs-1 segments
// The last segment is run on the current thread.  
// parallelForBoWait returns before checking that all thread are done.  parallelForFinish() must be
// called make sure all of the threads have finished.  
template< class RandomIt, class FN >
std::vector<std::future<void>>*  parallelForNoWait(const RandomIt begin, const RandomIt end, FN fn, int64_t num_segs = 1) {
  
  auto* futuresNW = new std::vector<std::future<void>>;

  // compute the number iterations and if 0, then return
  const int64_t n = end - begin;
  if (n == 0) return futuresNW;
  // if the number of elements is less than the number of segments or threads, reduce the number of segments
  if (n < num_segs) num_segs = n;
  // Compute the number of iterations to preform on each segment using a double to preserve the fraction.
  // using the method of iterating through bounds applied to each thread using a double and then rounding to
  // the nearest integer ensures that the difference in the number of iterations that each thread will do is 
  // no greater than one.
  double seg_size = (double)n / (double)num_segs;
  int64_t sBeg = 0;
  double sEnd = seg_size;
  // start threads for one less than the number of segments
  static std::mutex lock;
  for (int64_t seg = 0; seg < num_segs; seg++) {
    futuresNW->push_back(std::async(std::launch::async, [begin, fn](int64_t lb, int64_t le) {
      for (int64_t i = lb; i < le; i++) fn(begin + i);
      }, sBeg, llround(sEnd)));
    sBeg = llround(sEnd);
    sEnd += seg_size;
  }

  return futuresNW;
  }

void parallelForFinish(std::vector<std::future<void>>* futuresNW) {
  if (futuresNW == nullptr) return;
  // wait for the segments to complete
  for (size_t seg = 0; seg < futuresNW->size(); seg++)
    futuresNW->at(seg).wait();
  delete futuresNW;
}

#include <cmath>
void testParallelFor(size_t threads) {
  auto values = std::vector<double>(1000);

  parallelFor((size_t)0, values.size(),
    [&](int64_t i)
    {
      values[i] = std::sin(i * 0.001);
    }, threads);

  double total = 0;

  for (int i = 0;  i < values.size(); i++)
  {
    total += values[i];
    values[i] = 0;
  }

  std::cout << total << std::endl;

  parallelFor(values.begin(), values.end() ,
    [&](std::vector<double>::iterator vi)
    {
     *vi = std::sin((vi - values.begin()) * 0.001);
    }, threads);

  total = 0;

  for (int i = 0; i < values.size(); i++)
  {
    total += values[i];
    values[i] = 0;
  }

  std::cout << total << std::endl;

  double avalues[1000];

  parallelFor(avalues, avalues+1000,
    [&](double* vi)
    {
      *vi = std::sin((vi - avalues) * 0.001);
    }, threads);

  total = 0;

  for (int i = 0;  i < 1000; i++)
  {
    total += avalues[i];
    avalues[i] = 0;
  }

  std::cout << total << std::endl;

  auto f = parallelForNoWait((size_t)0, values.size(),
    [&](int64_t i)
    {
      values[i] = std::sin(i * 0.001);
    }, threads);

  parallelForFinish(f);

  total = 0;

  for (int i = 0; i < values.size(); i++)
  {
    total += values[i];
    values[i] = 0;
  }

  std::cout << total << std::endl;

}

#endif // THIS_PARALLEL_FOR
