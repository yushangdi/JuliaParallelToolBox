#include <iostream>
#include "gettime.h"

#define CILKP

#if defined(CILKP)
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <sstream>
#include <iostream>
#include <cstdlib>
#define parallel_for cilk_for

static int getWorkers() {
  return __cilkrts_get_nworkers();
}

static int getWorkerId() {
  return __cilkrts_get_worker_number();
}

#else
#define parallel_for for

static int getWorkers() {
  return 1;
}

static int getWorkerId() {
  return 0;
}


#endif

#define newA(__E,__n) (__E*) malloc((__n)*sizeof(__E))
#define nblocks(_n,_bsize) (1 + ((_n)-1)/(_bsize))

#define blocked_for(_i, _s, _e, _bsize, _body)  { \
    int _ss = _s;          \
    int _ee = _e;          \
    int _n = _ee-_ss;          \
    int _l = nblocks(_n,_bsize);     \
    parallel_for (int _i = 0; _i < _l; _i++) {   \
      int _s = _ss + _i * (_bsize);      \
      int _e = std::min(_s + (_bsize), _ee);      \
      _body           \
  }           \
  }


#define _SCAN_LOG_BSIZE 10
#define _SCAN_BSIZE (1 << _SCAN_LOG_BSIZE)

template <class T>
T prefixSumSerialSimple(T* data, int s, int e) {

  for (int i = s+1; i < e; ++i) {
    data[i] = data[i]+data[i-1];
  }
  return data[e-1];
}



template <class T>
T prefixSumSerial(T* data, int s, int e) {
  T res = 0;
  for (int i = s; i < e; ++i) {
    res += data[i];
    data[i] = res - data[i];
  }
  return res;
}

template <class T>
void addSerial(T* data, int s, int e, T val) {
  for (int i = s; i < e; ++i)
    data[i] += val;
}

template <class T>
T prefixSum(T* data, int s, int e) {
  int l = nblocks(e-s, _SCAN_BSIZE);
  if (l <= 1) return prefixSumSerial(data, s, e);
  T* sums = newA(T, l);
  blocked_for (i, s, e, _SCAN_BSIZE,
      sums[i] = prefixSumSerial<T>(data, s, e););
  T res = prefixSumSerial(sums, 0, l);
  blocked_for (i, s, e, _SCAN_BSIZE,
      addSerial(data, s, e, sums[i]););
  free(sums);
  return res;
}

void run_test(int N){
  std::cout << N/100000 << std::endl;
//
  int *arr = newA(int, N);
  for (int i = 0; i < N; ++ i){
    arr[i] = i;
  }
  timer t1;

  t1.start();
  prefixSumSerialSimple(arr, 0, N);
  std::cout << t1.next() << std::endl;

  for (int i = 0; i < N; ++ i){
    arr[i] = i;
  }

  t1.next();
  prefixSumSerial(arr, 0, N);
  std::cout << t1.next() << std::endl;


  for (int i = 0; i < N; ++ i){
    arr[i] = i;
  }

  t1.next();
  prefixSum(arr, 0, N);
  std::cout << t1.next() << std::endl;

  free(arr);
}

int main()
{
  std::cout << getWorkers() << std::endl;
  int N = 100000000;
  run_test(N);
}

//
// g++ -O3 -fcilkplus -ldl prefixsum.cpp -o a.out
// CILK_NWORKERS=4 ./a.out
