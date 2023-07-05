#ifndef PLOTS_H
#define PLOTS_H

#include "librarySPP.hpp"

#include <cmath>
// Modified version of https://github.com/alandefreitas/matplotplusplus
#include <matplot/matplot.h>
#include <matplot/util/common.h>

// Plot l'examen d'un run de RGA sur
// une instance
void plotRunRGA(
    const std::string instance,
    const std::vector<int>& ozInits,
    const std::vector<int>& pzAmels,
    const std::vector<int>& ozBests,
    const std::vector<float>& oprobas,
    const int zinit,
    const int zbest,
    const int num_iter,
    const int done_iter,
    std::string save_path = "",
    bool silent_mode = false);

// Plot le bilan de tous les runs de RGA
// sur une instance (plot exactement NUM_DIVISION
// points avec NUM_DIVISION <= NUM_ITER)
void plotAnalyseRGA(
    const std::string instance,
    const std::vector<double>& divs,
    const std::vector<int>& zMin,
    const std::vector<double>& zMoy,
    const std::vector<int>& zMax,
    const int allrunzmin,
    const float allrunzmoy,
    const int allrunzmax,
    std::string save_path = "",
    bool silent_mode = false);

// Plot le bilan CPUt pour chaque instance (le
// temps d'exécution moyen d'un run)
void plotCPUt(
    std::vector<std::string>& fnames,
    std::vector<float>& tMoy,
    std::string save_path = "",
    bool silent_mode = false);

// Plot des niveaux de phéromones avant
// et après le premier restart
void plotPhiRunACO(
    const std::string instance,
    const int n,
    const float* phi,
    bool before = true,
    std::string save_path = "",
    bool silent_mode = false);

// Plot des probabilités pour chaque valeur
// de alpha
void plotProbaRunGRASP(
    const std::string instance,
    const std::vector<double>& alpha,
    const std::vector<double>& proba,
    std::string save_path = "",
    bool silent_mode = false);

#endif /* end of include guard: PLOTS_H */
