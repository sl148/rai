/*  ------------------------------------------------------------------
    Copyright (c) 2011-2024 Marc Toussaint
    email: toussaint@tu-berlin.de

    This code is distributed under the MIT License.
    Please see <root-path>/LICENSE for details.
    --------------------------------------------------------------  */

#pragma once

#include "ConfigurationProblem.h"
#include "../Optim/NLP.h"
#include "../Algo/ann.h"

namespace rai {

struct PRM_PathFinder_Options {
  RAI_PARAM("prm/", int, verbose, 0)
  RAI_PARAM("prm/", int, nSamples, 500)
  RAI_PARAM("prm/", int, kNeighbors, 10)
  RAI_PARAM("prm/", double, collisionTolerance, 1e-4)
  RAI_PARAM("prm/", bool, useBroadCollisions, true)
  RAI_PARAM("prm/", int, subsamples, 4)
};

struct PRM_PathFinder : NonCopyable {
  PRM_PathFinder_Options opt;

  shared_ptr<ConfigurationProblem> P;
  ANN ann;
  rai::Array<uintA> adjacency; 
  
  //query info
  arr start, goal;

  //output
  arr path;
  shared_ptr<SolverReturn> ret;

  //setup
  void setProblem(const Configuration& C);
  void setStartGoal(const arr& _starts, const arr& _goals);
  void setExplicitCollisionPairs(const StringA& collisionPairs);

  //solve
  shared_ptr<SolverReturn> solve();

  //helpers
  void buildRoadmap();
  uint addNode(const arr& q, bool connect=true);
  bool checkEdge(const arr& q1, const arr& q2);
  void view(bool pause, const char* txt=0, bool play=false);

private:
  rai::Configuration DISP;
  void ensure_DISP();
};

} //namespace
