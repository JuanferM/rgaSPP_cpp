#include "heuristics.hpp"
#include "librarySPP.hpp"
#include <climits>
#include <cmath>
#include <cstdlib>

struct Data;

std::tuple<char*, int, char*> GreedyRandomized(
    int m,
    int n,
    const int* C,
    const char* A,
    const float* U,
    const float alpha,
    float* phi) {
  bool valid;
  int i(0), j(0), k(0), s(0), e(0), min_u(n-1), max_u(0);
  float limit(0.0f);
  std::vector<int> RCL;
  char *x = new char[n], *column = new char[m];
  for(i = 0; i < n; i++) x[i] = 0;
  for(j = 0; j < m; j++) column[j] = 0;

  std::vector<int> u_order = argsort(n, U); // DON'T FORGET TO DELETE

  k = 0;
  while(s != m && k < n) {
    // Indices of max and min utilities
    for(j = max_u, max_u = 0; j < n && u_order[max_u] == -1; j++)
      max_u = j;
    for(j = min_u, min_u = n-1; j >= 0 && u_order[min_u] == -1; j--)
      min_u = j;
    limit = U[u_order[min_u]] + alpha * (U[u_order[max_u]] - U[u_order[min_u]]);
    for(j = 0; j < n; j++)
      // If variable's index is candidate and the utility
      // greater than limit then we add the index in RCL
      if(u_order[j] != -1 && U[u_order[j]] >= limit)
         RCL.push_back(j);

    // Select an element e from RCL at random
    e = (RCL.size()) ? *select_randomly(RCL.begin(), RCL.end()) : max_u;
    for(j = 0, valid = true; j < m && valid; j++)
      valid = !(column[j] & A[INDEX(u_order[e], j)]);
    for(j = 0, s = 0; valid && j < m; s += column[j], j++)
      column[j] += A[INDEX(u_order[e], j)];
    x[u_order[e]] = valid, phi[u_order[e]] += valid, u_order[e] = -1;
    k += 1; RCL.clear();
  }

  return std::make_tuple(x, dot(n, x, C), column);
}

std::tuple<char*, int, char*> elaborateSolution(
    int m,
    int n,
    const int* C,
    const char* A,
    const float* U,
    float* phi,
    Data selection) {
  bool valid(true), valid2(true);
  int i(0), j(0), k(0), s(0);
  float sum_phi(0);
  std::vector<float> probs;
  char *x = new char[n], *column = new char[m], *phi_util = new char[n];
  for(i = 0; i < n; i++) {
    x[i] = 0;
    if(phi) phi_util[i] = /* U[i] * */ phi[i];
  }

  // Indices of utilities in utilities decreasing order
  std::vector<int> order = phi ? argsort(n, phi_util) : argsort(n, U);
  // We set the variable with the greatest utility to 1
  x[order[0]] = 1;
  // Selecting that variable means that we must select the
  // corresponding column in the matrix A and check if the
  // constraints are still verified
  for(j = 0; j < m; j++) {
    column[j] = A[INDEX(order[0], j)];
    s += column[j];
  }

  // Si on a bien l'adresse de P, on modifie sa valeur (mode sélection)
  if(selection.P) {
    *selection.P = !selection.iter ? selection.iter
           : log10(selection.iter) / log10(selection.maxIter);
    // Init roulette wheel probabilities
    sum_phi = 0, probs = std::vector<float>(n);
    for(k = 0; k < n; k++) {
      probs[k] = sum_phi + phi[k];
      sum_phi += phi[k];
    }

    probs[order[0]] = -1;
  }

  // Repeat the same process with each utility until constraints
  // are eventually violated
  i = 1;
  while(s != m && i < n) {
    if(selection.P && ((float)rand() / (float)RAND_MAX) > *selection.P) {
      float r = (float)rand() / (float)RAND_MAX;

      // Selon la roulette, on vérifie si une variable candidate existe
      for(k = 0, valid = false; k < n && !(valid && valid2); k++) {
        for(j = 0, valid = true; j < m && valid && valid2; j++) {
          valid = !(column[j] & A[INDEX(k, j)]);
          valid2 = (r < probs[k]/sum_phi);
        }

        // si la variable candidate existe, la prendre en compte
        for(j = 0, s = 0; j < m && valid && valid2; s += column[j], j++)
          column[j] += A[INDEX(k, j)];
        if(valid2) x[k] = valid, probs[k] = -1;
      }

      i++;
    } else {
      for(j = 0, valid = true; j < m && valid; j++)
        valid = !(column[j] & A[INDEX(order[i], j)]);
      for(j = 0, s = 0; j < m && valid; s += column[j], j++)
        column[j] += A[INDEX(order[i], j)];
      if(selection.P) probs[order[i]] = -1;
      x[order[i++]] = valid;
    }
  }

  return std::make_tuple(x, dot(n, x, C), column);
}

void GreedyImprovement(
    int m,
    int n,
    const int* C,
    const char* A,
    char* x,
    int* z,
    bool deep,
    char* column,
    float* phi) {
  int i(2);
  bool (*f[3])(int, int, const int*, const char*, char*, int*, bool, char*, float*) = {
      zero_oneExchange,
      one_oneExchange,
      two_oneExchange
    };

  // We modify x and z directly (no copy)
  while(i >= 0){
    if(!f[i](m, n, C, A, x, z, deep, column, phi)) i--;
  }
}

void managePheromones(
    const int n,
    const float rhoE,
    const float rhoD,
    const float phiNul,
    const int iter,
    const int maxIter,
    const int itStag,
    float* phi,
    float** phi_bef,
    float** phi_aft,
    char* xbest_iter,
    int* nbRestart,
    bool capturePhi) {
  int i(0);
  bool existPhiNul(false);

  for(i = 0; i < n; i++) {
    // Pheromone evaporation
    phi[i] = phi[i] * rhoE;
    // Pheromone deposition
    if(xbest_iter[i]) phi[i] = phi[i] + rhoD;
    if(phi[i] <= phiNul) existPhiNul = true;
  }

  // Territory disturbance
  if(itStag == 0 && existPhiNul) {
    (*nbRestart)++;
    // If phi_bef not initialized, copy phi into phi_bef
    if(!(*phi_bef) && capturePhi) {
      *phi_bef = new float[n];
      std::copy(phi, phi+n, *phi_bef);
    }
    int pn = rand() % (int)ceil(0.1 * n);

    // Disturb the pheromones
    for(i = 0; i < n; i++) {
      phi[i] = phi[i] * 0.95 * (!iter ? iter : log10(iter)/log10(maxIter));
      if(i < pn) {
        float r = (float)rand() / (float)RAND_MAX,
            d = (1.0 - (float)iter/(float)maxIter)*0.5 - 0.05;
        phi[rand() % n] = 0.05 + r * d;
      }
      // Offset on the pheromones with low level
      if(phi[i] < 0.1) {
        float r = (float)rand() / (float)RAND_MAX,
            d = (1.0 - (float)iter/(float)maxIter)*0.5 - 0.05;
        phi[i] = phi[i] + 0.05 + r * d;
      }
    }

    // If phi_aft not initialized, copy phi into phi_aft
    if(!(*phi_aft) && capturePhi) {
      *phi_aft = new float[n];
      std::copy(phi, phi+n, *phi_aft);
    }
  }
}

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
    int nbIter,
    bool deep,
    bool parallel) {
  int iter(0), zBest(-1), chunkLeft(probaUpdate), upd(0);
  double mean(0.0), diff(0.0), frac(0.0), sum(0.0), zmax(0.0), zmin(0.0);
  std::vector<std::vector<int>> pool(alpha.size(), std::vector<int>(probaUpdate));
  std::vector<double> valuation(pool.size(), 0.0);
  std::vector<int> poolData_i(probaUpdate, 0);
  std::vector<int> poolData_z(probaUpdate, 0);

  for(iter = 0; iter < nbIter; iter += chunkLeft) {
    if(iter + chunkLeft > nbIter) chunkLeft = nbIter-iter;

    #pragma omp parallel for if(parallel)
    for(upd = iter; upd < iter+chunkLeft; upd++) {
      char *x(nullptr), *column(nullptr);
      int i(0);
      float sel_alpha(-1.0), idx((double)rand() / RAND_MAX), s(0);
      for(i = 0; i < (int)proba.size() && sel_alpha == -1.0; i++) {
        s += proba[i];
        if(idx < s) sel_alpha = (float)alpha[i];
      }
      if(i == (int)proba.size()) i--;
      if(sel_alpha == -1.0) { // Make sure sel_alpha is well defined in any
        i = rand() % alpha.size(); // case
        sel_alpha = alpha[i];
      }
      std::tie(x, zInits[upd], column) = GreedyRandomized(m, n, C, A,
                        U, sel_alpha, phi);
      zAmels[upd] = zInits[upd];
      GreedyImprovement(m, n, C, A,
          x, &zAmels[upd], deep, column, phi);
      // Pool data (will help to reconstruct the pool after the parallel for)
      poolData_i[upd-iter] = i;
      poolData_z[upd-iter] = zAmels[upd];

      /* MOST IMPORTANT SECTION */
      if(x) delete[] x, x = nullptr;
      if(column) delete[] column, column = nullptr;
    }

    // Reconstruct pool
    for(upd = 0; upd < chunkLeft; upd++)
      pool[poolData_i[upd]].push_back(poolData_z[upd]);

    // Section de code difficilement parallélisable
    zmax = (double)*std::max_element(zAmels.begin(), zAmels.begin()+chunkLeft),
    zmin = (double)*std::min_element(zAmels.begin(), zAmels.begin()+chunkLeft);
    for(upd = 0, sum = 0.0; upd < (int)pool.size(); upd++, mean = 0.0) {
      for(double e : pool[upd]) mean += e;
      mean = pool[upd].size() ? mean/pool[upd].size() : zmin;
      diff = zmax - zmin;
      frac = diff ? (mean - zmin)/diff : diff;
      valuation[upd] = std::pow(std::abs(frac), (double)delta);
      sum += valuation[upd];
    }

    for(upd = 0; upd < (int)proba.size(); upd++)
      proba[upd] = sum ? valuation[upd]/sum : proba[upd];

    for(auto e : pool) e.clear(); // Clear each pool
  }

  zmax = (float)INT_MIN; // Use to compute max(phi)
  // Compute zBests using zAmels
  for(iter = 0; iter < nbIter; iter++) {
    zBest = std::max(zBest, zAmels[iter]);
    zBests[iter] = zBest;
    // Find max(phi) (part 1)
    if(iter < n && phi[iter] > zmax) zmax = phi[iter];
  }

  // Finish finding max(phi) (part 2)
  for(; iter < n; iter++)
    if(phi[iter] > zmax) zmax = phi[iter];

  // Compute phi
  for(iter = 0; iter < n; iter++) {
    float r = rand()/(float)RAND_MAX,
        d = phiOffset - phiOffset/4;
    if(phi[iter]/zmax < phiOffset)
      phi[iter] = phiOffset/4 + r*d + phi[iter]/zmax;
    else phi[iter] = phi[iter]/zmax;
    if(phi[iter] > 1) phi[iter] = 1;
  }
}

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
    int maxAnts,
    int maxIter,
    int numIter,
    bool deep,
    bool parallel,
    bool restartStop,
    int maxRestart,
    bool capturePhi) {
  float P(0); // ant's curiosity level
  bool keep_going(true);
  int zbest(zGRASP), zbest_iter;
  int nbRestart(0), iter(0), itStag(iterStagnant);
  char *xbest(nullptr), *xbest_iter(nullptr);

  for(iter = 0; iter < maxIter && keep_going; iter++) {
    xbest_iter = new char[n], zbest_iter = -1;

    #pragma omp parallel for if(parallel)
    for(int ant = 0; ant < maxAnts; ant++) {
      char *x(nullptr), *column_ant(nullptr);
      bool condition = ((float)rand() / (float)RAND_MAX) >= P;
      Data selection = Data(iter, maxIter, &P);

      std::tie(x, zInits[numIter+iter*maxAnts+ ant], column_ant) =
        condition ? elaborateSolution(m, n, C, A, U, phi, selection)
              : elaborateSolution(m, n, C, A, U, phi);
      pExploit[numIter+iter*maxAnts+ant] = P;
      zAmels[numIter+iter*maxAnts+ant] = zInits[numIter+iter*maxAnts+ant];
      GreedyImprovement(m, n, C, A, x, &zAmels[numIter+iter*maxAnts+ant],
          deep, column_ant);
      if(column_ant) delete[] column_ant, column_ant = nullptr;


      #pragma omp critical
      if(zAmels[iter * maxAnts + ant] > zbest_iter) {
        std::copy(x, x+n, xbest_iter);
        zbest_iter = zAmels[numIter+iter*maxAnts+ant];
        if(zbest_iter > zbest) {
          itStag = iterStagnant;
          zbest = zbest_iter;
          if(!xbest) xbest = new char[n];
          std::copy(xbest_iter, xbest_iter+n, xbest);
        }
      }

      if(x) delete[] x, x = nullptr;
    }

    managePheromones(n, rhoE, rhoD, phiNul, iter, maxIter, itStag--,
            phi, phi_bef, phi_aft, xbest_iter, &nbRestart,
            capturePhi);
    if(itStag <= 0) itStag = 0;
    if(restartStop && nbRestart == maxRestart) keep_going = false;
    if(xbest_iter) delete[] xbest_iter, xbest_iter = nullptr;
  }

  // Compute zBests using zAmels
  zbest_iter = std::max(zGRASP, zAmels[numIter]);
  zBests[numIter] = zbest_iter;
  for(itStag = 0; itStag < maxIter*maxAnts; itStag++) {
    zbest_iter = std::max(zbest_iter, zAmels[numIter+itStag]);
    zBests[numIter+itStag] = zbest_iter;
  }

  return numIter+iter*maxAnts;
}
