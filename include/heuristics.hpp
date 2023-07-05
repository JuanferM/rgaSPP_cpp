#ifndef HEURISTICS_H
#define HEURISTICS_H

#include "movements.hpp"

#include <omp.h>
#include <cmath>
#include <cstdlib>
#include <math.h>
#include <string>
#include <cstring>

struct Data {
  Data() : iter(-1), maxIter(-1), P(nullptr) {}
  Data(int i, int m, float* p) : iter(i), maxIter(m), P(p) {}
  int iter, maxIter;
  float* P;
};

// Greedy randomized construction of a feasible solution
std::tuple<char*, int, char*> GreedyRandomized(
    int m,
    int n,
    const int* C,
    const char* A,
    const float* U,
    const float alpha,
    float* phi = nullptr);

// Elaborate solution for ACO
std::tuple<char*, int, char*> elaborateSolution(
    int m,
    int n,
    const int* C,
    const char* A,
    const float* U,
    float* phi = nullptr,
    Data selection = Data());

// Greedy improvement of a feasible solution through (deep) local search
void GreedyImprovement(
    int m,
    int n,
    const int* C,
    const char* A,
    char* x,
    int* z,
    bool deep = true,
    char* column = nullptr,
    float* phi = nullptr);

// ReactiveGRASP for the Set Packing Problem
void ReactiveGRASP(
    const int m,
    const int n,
    const int* C,
    const char* A,
    const float* U,
    std::vector<int>& zInits,
    std::vector<int>& zAmels,
    std::vector<int>& zBests,
    float* phi,
    float phiOffset,
    const std::vector<double>& alpha,
    std::vector<double>& proba,
    const int probaUpdate,
    const double delta,
    int nbIter = 100,
    bool deep = true,
    bool parallel = true);

// ACO for the Set Packing Problem
int ACO(
    const int m,
    const int n,
    const int* C,
    const char* A,
    const float* U,
    const int zGRASP,
    std::vector<int>& zInits,
    std::vector<int>& zAmels,
    std::vector<int>& zBests,
    std::vector<float>& pExploit,
    float* phi,
    float** phi_bef,
    float** phi_aft,
    const float phiNul,
    const float rhoE,
    const float rhoD,
    const int iterStagnant,
    int maxAnts = 15,
    int maxIter = 100,
    int numIter = 100,
    bool deep = true,
    bool parallel = true,
    bool restartStop = true,
    int maxRestart = 2,
    bool capturePhi = false);

#endif /* end of include guard: HEURISTICS_H */
