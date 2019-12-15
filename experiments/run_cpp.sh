#!/bin/bash

g++ -O3 -fcilkplus -ldl prefixsum.cpp -o a.out
CILK_NWORKERS=$1 ./a.out

