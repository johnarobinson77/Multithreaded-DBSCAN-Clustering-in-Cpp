#pragma once

/**
 * Copyright (c) 2024 John Robinson.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */


#ifndef SIMPLELIST_H
#define SIMPLELIST_H

#include <list>
#include <mutex>
#include <exception>



/* 
* SimpleList is a list that is expected to be used as a FIFO.  It only implements the following
* push_back(new_value) to add the FIFO
* pop_front(&value) to simultaneously get the fist value and remove it from the FIFO.
* front() to get the fist value but not pop it off
* reset() to put the FIFO in an empty state.  This should be called whenever the FIFO is known to be empty.
* size() to get the number of values in the FIFO
* empty() returns size = 0;
* In the use case for doing the cluster search, the FIFO is expected to be highly dynamic.
* To reduce the pressure on the memory system, Memory is allocated in blocks and reused 
* for its entire lifetime.  The memory is only freed when the list is deleted.
*/
template <typename _LT>
class SimpleList {
private:
  struct node {
    _LT value;
    node* nextNode;
  };

  // memory allocation fields;
  size_t allocationBlockSize;  // size of memory blocks allocated
  std::list<uint8_t*> allocationBlockPtrs;  // list of pointers to allocated blocks
  uint8_t* lastNodePtr = nullptr;      // pointer to the last node allocated.
  uint8_t* nextNodePtr = nullptr;      // pointer to the next node to be allocated.
  const size_t alignment = 8;          // byte alignment of the node allocations.
  const size_t paddedNodeSize = (sizeof(node) + (alignment - 1) & ~(alignment - 1));

  // Returns pointer to an internal allocation of a node.
  // if there is not enough memory in the current allocation block, malloc another block.
  node* nodeAlloc() {
    lastNodePtr = nextNodePtr;
    nextNodePtr += paddedNodeSize;
    if (nextNodePtr - allocationBlockPtrs.back() > static_cast<int64_t>(allocationBlockSize * paddedNodeSize)) {
      lastNodePtr = (uint8_t*)std::malloc(allocationBlockSize * paddedNodeSize);
      if (lastNodePtr == nullptr)
        throw std::runtime_error("Ran out of memory.");
      allocationBlockPtrs.push_back(lastNodePtr);
      nextNodePtr = lastNodePtr + paddedNodeSize;
    }
    return (node*)lastNodePtr;
  }
  // free all the allocation blocks.
  void freeMemory() {
    for (void* ptr : allocationBlockPtrs) {
      free(ptr);
    }
  }

  // node read/write pointers and counters
  size_t lastInCount = 0;
  node* lastInPtr = nullptr;
  size_t lastOutCount = 0;
  node* lastOutPtr = nullptr;

  // called only by constructors to set up the memory
  void listInit(size_t blockSize) {
    allocationBlockSize = blockSize;
    nextNodePtr = (uint8_t*)std::malloc(allocationBlockSize * paddedNodeSize);
    if (nextNodePtr == nullptr)
      throw std::runtime_error("Ran out of memory.");
    allocationBlockPtrs.push_back(nextNodePtr);
    reset();
  }

public:

  SimpleList() {
    listInit(1000);
  }

  SimpleList(size_t blockSize) {
    listInit(blockSize);
  }

  ~SimpleList() {
    freeMemory();
  }

  void reset() {
    lastInCount = 0;
    lastInPtr = nullptr;
    lastOutCount = 0;
    lastOutPtr = nullptr;
  }

  void push_back(_LT value) {
    // fist check to see if a new node is needed and allocate it.
    if (lastInPtr == nullptr) {
      if (lastNodePtr == nullptr)
        lastInPtr = nodeAlloc();
      else
        lastInPtr = (node*)allocationBlockPtrs.front();
      lastInPtr->value = value;
      lastInCount++;
      return;
    }
    else if (lastInPtr == (node*)lastNodePtr) {
      node* nextInPtr = nodeAlloc();
      lastInPtr->nextNode = nextInPtr;
      nextInPtr->nextNode = nullptr;
      lastInPtr = nextInPtr;
    }
    // else just move to the next node.
    else {
      lastInPtr = lastInPtr->nextNode;
    }
    // store the value and increment the counter.
    lastInPtr->value = value;
    lastInCount++;
  }

  bool pop_front(_LT& value) {
    if (lastOutCount >= lastInCount)
      return false;
    lastOutCount++;
    if (lastOutPtr == nullptr) {
      lastOutPtr = (node*)allocationBlockPtrs.front();
    }
    else {
      lastOutPtr = lastOutPtr->nextNode;
    }
    value = lastOutPtr->value;
    return true;
  }

  _LT front() {
    if (lastOutCount + 1 > lastInCount)
      return NULL;
    return lastOutPtr->nextNode->value;
  }

  size_t size() {
    return lastInCount - lastOutCount;
  }

  bool empty() {
    return (lastOutCount >= lastInCount);
  }

  size_t max_size() {
    return allocationBlockSize * allocationBlockPtrs.size();
  }

}; // class SimpleList

#endif // SIMPLELIST_H


#ifndef SIMPLEALLOC_H
#define SIMPLEALLOC_H

/*
* SimpleAlloc is a linear memory allocator that is useful for grouping large numbers
* of small allocations into larger system memory blocks, this relieving pressure on 
* system memory.  THe memory it allocates is not freed until the SimpleAlloc object is deleted
* releasing all of the memory at once.  In this use case it is used to store the KdTree
* which is deleted after the cluster search is complete.
*/

class SimpleAlloc {
private:
  // memory allocation fields;
  size_t allocationBlockSize;  // size of memory blocks allocated
  std::list<uint8_t*> allocationBlockPtrs;  // list of pointers to allocated blocks
  uint8_t* nextNodePtr = nullptr;      // pointer to the next node to be allocated.
  std::mutex oneAtATime;

  // Returns pointer to an allocation of size.
  // if there is not enough memory in the current allocation block, malloc another block.
public:
  void* allocate(size_t size, size_t alignment) {
    // lock the allocator
    oneAtATime.lock();
    //calculate size of actual allocation to include padding.
    uint8_t* lastNodePtr = (uint8_t *)(((size_t)nextNodePtr + (alignment - 1)) & ~(alignment - 1));
    nextNodePtr = lastNodePtr + size;
    // if this allocation request would extend past the current block, allocate another block.
    if (nextNodePtr - allocationBlockPtrs.back() > static_cast<int64_t>(allocationBlockSize) ) {
      lastNodePtr = (uint8_t*)std::malloc(allocationBlockSize);
      if (lastNodePtr == nullptr)
        throw std::runtime_error("Ran out of memory.");
      allocationBlockPtrs.push_back(lastNodePtr);
      nextNodePtr = lastNodePtr + size;
    }
    // unlock and return
    oneAtATime.unlock();
    return (void*)lastNodePtr;
  }

  // free all the allocation blocks.
public:  
  void freeMemory() {
    for (void* ptr : allocationBlockPtrs) {
      free(ptr);
    }
  }

  // called only by constructors to set up the memory
private:
  void allocInit(size_t blockSize) {
    allocationBlockSize = blockSize;
    nextNodePtr = (uint8_t*)std::malloc(allocationBlockSize);
    if (nextNodePtr == nullptr)
      throw std::runtime_error("Ran out of memory.");
    allocationBlockPtrs.push_back(nextNodePtr);
  }

public:

  SimpleAlloc() {
    allocInit(10000000);
  }

  SimpleAlloc(size_t blockSize) {
    allocInit(blockSize);
  }

  ~SimpleAlloc() noexcept {
    freeMemory();
  }

  void deallocate(void* p) {

  }

  SimpleAlloc(SimpleAlloc&& other) noexcept {
    allocationBlockSize = other.allocationBlockSize;
    allocationBlockPtrs = other.allocationBlockPtrs;
    nextNodePtr = other.nextNodePtr;
    other.allocationBlockPtrs.clear();
    other.nextNodePtr = nullptr;
  }

  SimpleAlloc& operator=(SimpleAlloc&& rhs) noexcept
  {
    allocationBlockSize = rhs.allocationBlockSize;
    allocationBlockPtrs = rhs.allocationBlockPtrs;
    nextNodePtr = rhs.nextNodePtr;
    rhs.allocationBlockPtrs.clear();
    rhs.nextNodePtr = nullptr;
    return *this;
  }

  void* getBase() noexcept {
    return (void*)allocationBlockPtrs.front();
  }

  size_t getMaxSize() noexcept {
    return allocationBlockPtrs.size() * allocationBlockSize;
  }

}; // SimpleAlloc


/*
* STLAdaptor is an class the wraps a SimpleAllocator for use in STL containers.
* In this use case it is provided so that the containers in the KdTree object
* use the same memory and the nodes in the KdTree.
*/
template<typename T>
class STLAdaptor
{
public:

  typedef T value_type;

  SimpleAlloc& simpleAlloc;

  STLAdaptor() = delete;

  STLAdaptor(SimpleAlloc& alloc) : 
    simpleAlloc(alloc)
  { }

  template<typename U>
  STLAdaptor(const STLAdaptor < U >& other) : simpleAlloc(other.simpleAlloc) {}

  [[nodiscard]] constexpr T* allocate(std::size_t n) {
    return reinterpret_cast<T*>
      (simpleAlloc.allocate(n * sizeof(T), alignof(T)));
  }

  constexpr void deallocate(T* p, [[maybe_unused]] std::size_t n) noexcept {
    simpleAlloc.deallocate(p);
  }

  std::size_t MaxAllocationSize() const noexcept {
    return simpleAlloc.getMaxSize();
  }

  bool operator==(const STLAdaptor<T>& rhs) const noexcept {
    return simpleAlloc.getBase() == rhs.simpleAlloc.getBase();
  }

  bool operator!=(const STLAdaptor<T>& rhs) const noexcept {
    return !(*this == rhs);
  }

}; // STLAdaptor

#endif // SIMPLEALLOC_H


#ifndef SIMPLESET_H
#define SIMPLESET_H

/*
* SimpleSet is a set container into which a value of a uint32_t integer can be inserted 
* once and only once.  Also a test can be run to see if a value is already in
* In this way it is similar to std::set<uint32_t> except that it is faster for small sets on
* the order of 128 or less.
* This version of SimpleSet uses the SimpleAlloc for memory allocation.  In the use case of
* the node of a KdTree it uses the same allocator as the rest of the node. 
*/

class SimpleSet {
private:
  SimpleAlloc* alloc = nullptr;  // pointer to the memory allocator
  uint32_t* container = nullptr;  // pointer to the array that holds the values
  size_t containerSize = 0;  // the current size of the container.
  size_t nextIn = 0;  // the number of values in the set and where the next entry will be placed.

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

  // push_back adds to values to the container and grows the
  // container size when needed.
  void push_back(uint32_t newValue) {
    // first check to see if there is enough memory
    if (nextIn == containerSize) {
      uint32_t *oldContainer = container;
      size_t oldContainerSize = containerSize;
      containerSize *= 4;
      container = (uint32_t*)alloc->allocate(sizeof(uint32_t) * containerSize, alignof(uint32_t));
      for (size_t i = 0; i < nextIn; i++) {
        container[i] = oldContainer[i];
      }
    }
    container[nextIn] = newValue;
  }

public:
  // Required
  // setup() initializes the memory for a SimpleSet object.  It must be 
  // called before inserting any values to provide a pointer to the SimpleAlloc allocator.
  void setup(SimpleAlloc* newAlloc) {
    alloc = newAlloc;
    nextIn = 0;
    containerSize = 1;
    container = (uint32_t*)alloc->allocate(sizeof(uint32_t) * containerSize, alignof(uint32_t));
  }

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
      push_back(in);
      nextIn++;
    }
    else {
      if (container[idx] < in) idx++;
      push_back(container[cs - 1]); // increase the size of the vector with the current last value
      nextIn++;
      for (size_t i = cs - 1; i > idx; i--) // copy the rest of the values up to the insertion point to make a space for the new value
        container[i] = container[i - 1];
      container[idx] = in; // copy the new value in.
    }
    return true;
  }

  // print() print all the values in the container.
  void print() {
    for (size_t i = 0; i < nextIn; i++) {
      std::cout << container[i] << ", ";
    }
    std::cout << std::endl;
  }

  // size() returns the number of items in the container.
  size_t const size() const {
    return nextIn;
  }

  // value() returns the value of the idx'th entry in container.  If idx > size()-1, maximum uint32_t is returned.
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
    nextIn = 0;
    containerSize = 1;
    container = (uint32_t*)alloc->allocate(sizeof(uint32_t) * containerSize, alignof(uint32_t));
  }
};

#endif // SIMPLESET_H

