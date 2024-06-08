#pragma once

/**
 * Copyright (c) 2024 John Robinson.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

#ifndef SIMPLESET_H
#define SIMPLESET_H

#include <vector>
#include <mutex>

// SempleSet is a set container into which a value of a uint32_t integer can be inserted 
// once and only once.  Also a test can be run to see if a value is already in
// In this way it is similar to std::set<uint32_t> except that it is faster for small sets on
// the order of 128 or less.   
class SimpleSet {
private:
  std::vector<uint32_t> container = std::vector<uint32_t>(1);
  size_t nextIn = 0;

  // find() does a binary search of the values in container to find the index of where the value sv should be inserted.
  // The index is returned in idxRef;  If the value was already in container, find() returns true, otherwise false; 
  bool find(const uint32_t sv, size_t& idxRef) {
    size_t ub = nextIn;
    if (sv > container[ub - 1]) {
      idxRef = ub;
      return false;
    }
    size_t idx = 0;
    size_t lb = 0;
    bool ret = false;
    size_t k = ub;
    for (size_t k = ub; k > 0; k >>= 1) {
      idx = lb + ((ub - lb) >> 1);
      if (container[idx] == sv) {
        ret = true;
        break;
      }
      else if (container[idx] < sv) lb = idx;
      else ub = idx;
    }
    idxRef = idx;
    return ret;
  }

public:
  // contains() returns true if the value sv is in container or false if not.
  // call with a mutex reference for multithreaded protection
  // It's done with a binary search so complexity is O(log n). 

  inline bool contains(const uint32_t sv, std::mutex& mx) {
    const std::lock_guard<std::mutex> lock(mx);
    return contains(sv);
  }

  bool contains(const uint32_t sv) {
    size_t ub = nextIn;
    if (ub == 0ULL || sv > container[ub - 1]) {
      return false;
    }
    size_t idx = 0;
    size_t lb = 0;
    bool ret = false;
    for (size_t k = ub; k > 0; k >>= 1) {
      idx = lb + ((ub - lb) >> 1);
      if (container[idx] == sv) {
        ret = true;
        break;
      }
      else if (container[idx] > sv) ub = idx;
      else lb = idx;
    }
    return ret;
  }

  // insert() attemps to add the value in to the container.  If the value is not in container, 
  // the value is added and insert() returns true;  If the value is in the container, nothing
  // is added and insert() returns false;
  // call with a mutex reference for multithreaded protection
  // Complexity is O(log n) to find the index of the insert point plus O(N/2) if an insert it performed.
  inline bool insert(uint32_t in, std::mutex& mx) {
    const std::lock_guard<std::mutex> lock(mx);
    return insert(in);
  }

  bool insert(uint32_t in) {
    if (nextIn == 0) {
      container[0] = in;
      nextIn++;
      return true;
    }
    size_t idx;
    bool ret = find(in, idx);
    if (ret) return false;  // container already has the value
    size_t cs = nextIn;
    if (idx >= cs) {
      container.push_back(in);
      nextIn++;
    }
    else {
      if (container[idx] < in) idx++;
      container.push_back(container[cs - 1]); // increase the size of the vector with the current last value
      nextIn++;
      for (size_t i = cs - 1; i > idx; i--) // copy the rest of the values up to the insertion point to make a space for the new value
        container[i] = container[i - 1];
      container[idx] = in; // copy the new value in.
    }
    return true;
  }

  // print() print all the values in the container.
  void print() {
    for (size_t i = 0; i < container.size(); i++) {
      std::cout << container[i] << ", ";
    }
    std::cout << std::endl;
  }

  // size() returns the number of items in the container.
  size_t const size() const {
    return nextIn;
  }

  // value returns the value of the idx'th entry in container.  If idx > size()-1, maximum uint32_t is returned.
  // call with a mutex reference for multithreaded protection

  inline uint32_t value(const size_t idx, std::mutex& mx) {
    if (idx >= nextIn) return 0xFFFFFFFF;
    const std::lock_guard<std::mutex> lock(mx);
    return container[idx];
  }

  inline uint32_t value(const size_t idx) {
    if (idx >= nextIn) return 0xFFFFFFFF;
    return container[idx];
  }

  // clear empties the set
  // call with a mutex reference for multithreaded protection
  void clear(uint32_t in, std::mutex& mx) {
    const std::lock_guard<std::mutex> lock(mx);
    clear();
  }

  void clear() {
    container.clear();
    container.push_back(0UL);
    nextIn = 0;
  }

};

#endif // SIMPLESET_H

