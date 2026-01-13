/*  ------------------------------------------------------------------
    Copyright (c) 2011-2024 Marc Toussaint
    email: toussaint@tu-berlin.de

    This code is distributed under the MIT License.
    Please see <root-path>/LICENSE for details.
    --------------------------------------------------------------  */

#include "PRM_PathFinder.h"

#include "../Gui/opengl.h"
#include "../Kin/viewer.h"
#include "../KOMO/pathTools.h"

#include <queue>
#include <map>

namespace rai {

// Copied from RRT_PathFinder.cpp to ensure independence
// Renamed to avoid collision in unity build
double prm_corput(uint n, uint base) {
  double q = 0.;
  double bk = 1./double(base);

  while(n > 0) {
    q += (n % base)*bk;
    n /= base;
    bk /= double(base);
  }
  return q;
}

bool prm_checkConnection(ConfigurationProblem& P,
                     const arr& start,
                     const arr& end,
                     const uint num,
                     const bool binary) {
  if(binary) {
    for(uint i=1; i<num; ++i) {
      double ind = prm_corput(i, 2);
      arr p = start + ind * (end-start);
      if(!P.query(p)->isFeasible) return false;
    }
  } else {
    for(uint i=1; i<num-1; ++i) {
      arr p = start + 1.0 * i / (num-1) * (end-start);
      if(!P.query(p)->isFeasible) return false;
    }
  }
  return true;
}

// Helpers for spherical coordinates (Copied from RRT_PathFinder.cpp)
// Renamed to avoid collision in unity build
void prm_normalizeSphericalCoordinates(arr& x, const uintA& idx){
  arr xsub = x({idx(0), idx(0)+idx(1)-1});
  op_normalize(xsub);
}

void prm_randomSphericalCoordinates(arr& x, const uintA& idx){
  arr xsub = x({idx(0), idx(0)+idx(1)-1});
  xsub = randn(xsub.N);
  op_normalize(xsub);
}

void prm_flipSphericalCoordinates(arr& x, const uintA& idx){
  arr xsub = x({idx(0), idx(0)+idx(1)-1});
  xsub *= -1.;
}

void PRM_PathFinder::setProblem(const Configuration& C){
  // Manually sync options from global params just in case
  opt.collisionTolerance = rai::getParameter<double>("prm/collisionTolerance", opt.collisionTolerance);
  opt.nSamples = rai::getParameter<int>("prm/nSamples", opt.nSamples);
  opt.kNeighbors = rai::getParameter<int>("prm/kNeighbors", opt.kNeighbors);

  P = make_shared<ConfigurationProblem>(C, opt.useBroadCollisions, opt.collisionTolerance, 1);
  P->verbose=0;
  LOG(0) << "PRM Problem set. Tolerance: " << opt.collisionTolerance;
}

void PRM_PathFinder::setStartGoal(const arr& _starts, const arr& _goals){
  start = _starts;
  goal = _goals;
  
  LOG(0) << "PRM setStartGoal. P->collisionTolerance: " << P->collisionTolerance;

  // Verify feasibility
  auto q0ret = P->query(start);
  auto qTret = P->query(goal);
  if(!q0ret->isFeasible) { LOG(0) <<"initializing with infeasible start: " << q0ret->totalCollision; q0ret->writeDetails(std::cout, *P); }
  if(!qTret->isFeasible) { LOG(0) <<"initializing with infeasible goal: " << qTret->totalCollision; qTret->writeDetails(std::cout, *P); }
}

void PRM_PathFinder::setExplicitCollisionPairs(const StringA& collisionPairs) {
  CHECK(P, "need to set problem first");
  P->setExplicitCollisionPairs(collisionPairs);
}

uint PRM_PathFinder::addNode(const arr& q, bool connect){
  auto ret = P->query(q);
  if(!ret->isFeasible){
      LOG(0) << "Failed to add node. Collision: " << ret->totalCollision << " Tolerance: " << P->collisionTolerance;
      return -1;
  }

  uint idx = ann.X.d0;
  ann.append(q);
  LOG(0) << "Added node " << idx << ". Total nodes: " << ann.X.d0;

  // If spherical, add flipped versions
  if(P->sphericalCoordinates.N){
      arr q_flipped = q;
      for(uint i=0; i<P->sphericalCoordinates.d0; i++) {
          prm_flipSphericalCoordinates(q_flipped, P->sphericalCoordinates[i]);
          ann.append(q_flipped);
      }
  }

  if(connect){
      adjacency.resize(ann.X.d0);
      uint k = opt.kNeighbors;
      if(k >= ann.X.d0) k = ann.X.d0 - 1;
      
      LOG(0) << "Connecting node " << idx << " to " << k << " neighbors.";

      // Connect newly added nodes (idx to ann.X.d0-1)
      for(uint i=idx; i<ann.X.d0; i++){
          uintA idxs;
          arr dists;
          ann.getkNN(dists, idxs, ann.X[i], k); // Finds kNN in the WHOLE tree
          
          LOG(0) << "Node " << i << " found " << idxs.N << " neighbors.";
          
          uint connectedCount = 0;
          for(uint j : idxs){
              if(i == j) continue;
              
              // Check existing?
              bool exists = false;
              for(uint neighbor : adjacency(i)){ if(neighbor == j) { exists=true; break; } }
              if(exists) continue;
              
              if(checkEdge(ann.X[i], ann.X[j])){
                  adjacency(i).append(j);
                  adjacency(j).append(i);
                  connectedCount++;
              }
          }
          LOG(0) << "Node " << i << " connected to " << connectedCount << " new neighbors.";
      }
  }
  
  return idx;
}

void PRM_PathFinder::buildRoadmap(){
  ann.clear();
  adjacency.clear();

  // 1. (Formerly) Add Start and Goal first -- NOW SKIPPED in buildRoadmap

  uint dim = P->C.getJointStateDimension();

  // 2. Sample random nodes
  uint n = opt.nSamples;
  for(uint i=0; i<n; ++i){
    arr q(dim);
    // Random sampling based on limits
    for(uint j=0; j<q.N; j++) {
      double lo=P->limits(0, j), up=P->limits(1, j);
      if(up-lo < 1e-6) { // Fixed joint or small range
          q(j) = lo;
      } else {
          q(j) = lo + rnd.uni()*(up-lo);
      }
    }
    
    // Handle spherical coordinates
    for(uint j=0; j<P->sphericalCoordinates.d0; j++) {
        prm_randomSphericalCoordinates(q, P->sphericalCoordinates[j]);
    }
    
    addNode(q, false);
  }
  
  // Resize adjacency
  adjacency.resize(ann.X.d0);
  LOG(0) << "Roadmap nodes: " << ann.X.d0;
  
  if(ann.X.d0 > 2){
      LOG(0) << "First few roadmap configurations:";
      for(uint i=0; i<std::min(ann.X.d0, (uint)5); i++){
          LOG(0) << i << ": " << ann.X[i];
      }
  }

  // 3. Connect nodes (Batch)
  uint k = opt.kNeighbors;
  if(ann.X.d0 > 0 && k >= ann.X.d0) k = ann.X.d0 - 1;
  
  std::cout << "DEBUG: Connecting roadmap nodes (k=" << k << ")..." << std::endl;

  uint totalEdges = 0;
  for(uint i=0; i<ann.X.d0; ++i){
    uintA idxs;
    arr dists;
    ann.getkNN(dists, idxs, ann.X[i], k);
    
    for(uint j : idxs){
      if(i == j) continue;
      
      // Check if edge already exists (undirected)
      bool exists = false;
      for(uint neighbor : adjacency(i)){
        if(neighbor == j) { exists = true; break; }
      }
      if(exists) continue;

      if(checkEdge(ann.X[i], ann.X[j])){
        adjacency(i).append(j);
        adjacency(j).append(i);
        totalEdges++;
      }
    }
  }
  LOG(0) << "Roadmap construction complete. Total edges: " << totalEdges;
}

bool PRM_PathFinder::checkEdge(const arr& q1, const arr& q2){
    return prm_checkConnection(*P, q1, q2, opt.subsamples, true);
}

shared_ptr<SolverReturn> PRM_PathFinder::solve(){
  if(!ret) ret = make_shared<SolverReturn>();
  P->useBroadCollisions = opt.useBroadCollisions;

  ret->time -= rai::cpuTime();

  // Build roadmap if not already built
  if(ann.X.d0 == 0) buildRoadmap();

  // Ensure start and goal are in the map
  uint startNode = (uint)-1;
  uint goalNode = (uint)-1;

  // Try to find existing start
  if(ann.X.d0 > 0) {
      uint nn = ann.getNN(start);
      if(length(ann.X[nn] - start) < 1e-6) startNode = nn;
  }
  if(startNode == (uint)-1) {
      LOG(0) << "Adding Start Node...";
      startNode = addNode(start, true);
  } else {
      LOG(0) << "Start Node found at " << startNode;
  }

  // Try to find existing goal
  if(ann.X.d0 > 0) {
      uint nn = ann.getNN(goal);
      if(length(ann.X[nn] - goal) < 1e-6) goalNode = nn;
  }
  if(goalNode == (uint)-1) {
       LOG(0) << "Adding Goal Node...";
       goalNode = addNode(goal, true);
  } else {
      LOG(0) << "Goal Node found at " << goalNode;
  }

  // Check feasibility of start/goal
  if(startNode == (uint)-1){
       LOG(0) << "Start configuration is infeasible!";
       ret->feasible = false;
       ret->time += rai::cpuTime();
       return ret;
  }
  if(goalNode == (uint)-1){
       LOG(0) << "Goal configuration is infeasible!";
       ret->feasible = false;
       ret->time += rai::cpuTime();
       return ret;
  }
  
  // Debug Adjacency
  LOG(0) << "Start Node " << startNode << " degree: " << adjacency(startNode).N;
  LOG(0) << "Goal Node " << goalNode << " degree: " << adjacency(goalNode).N;


  // Dijkstra
  // Verify distances are zero (or close)
  if(length(ann.X[startNode] - start) > 1e-6) {
      LOG(-1) << "PRM: Start node not found in roadmap?";
  }
  if(length(ann.X[goalNode] - goal) > 1e-6) {
      LOG(-1) << "PRM: Goal node not found in roadmap?";
  }
  
  // Dijkstra Search
  uint N = ann.X.d0;
  arr dist = zeros(N);
  for(uint i=0; i<N; i++) dist(i) = -1.0; // -1 as infinity
  
  uintA parent(N);
  for(uint i=0; i<N; i++) parent(i) = i;

  // Priority queue: <distance, node_index>
  std::priority_queue<std::pair<double, uint>> pq;
  
  dist(startNode) = 0.0;
  pq.push({0.0, startNode});

  bool success = false;
  uint reachedGoalNode = 0;
  
  while(!pq.empty()){
    double d = -pq.top().first;
    uint u = pq.top().second;
    pq.pop();

    if(d > dist(u) && dist(u) != -1.0) continue;
    
    if(u == goalNode){
      success = true;
      reachedGoalNode = u;
      break;
    }

    for(uint v : adjacency(u)){
      double weight = length(ann.X[u] - ann.X[v]);
      if(dist(v) == -1.0 || dist(u) + weight < dist(v)){
        dist(v) = dist(u) + weight;
        parent(v) = u;
        pq.push({-dist(v), v});
      }
    }
  }

  path.clear();
  if(success){
    uint curr = reachedGoalNode;
    while(curr != startNode){
      path.append(ann.X[curr]);
      curr = parent(curr);
    }
    path.append(ann.X[startNode]);
    
    // Reshape to ensure it is 2D [T, dim]
    path.reshape(-1, ann.X.d1);

    // Reverse path
    uint nPath = path.d0;
    arr temp;
    for(uint i=0; i<nPath/2; i++){
      temp = path[i];
      path[i] = path[nPath-1-i];
      path[nPath-1-i] = temp;
    }
  } else {
      LOG(0) << "Dijkstra failed to find path between node " << startNode << " and " << goalNode;
  }

  ret->time += rai::cpuTime();
  ret->done = true;
  ret->feasible = success;
  ret->x = path;
  ret->evals = N; // Number of nodes in graph

  return ret;
}

void PRM_PathFinder::view(bool pause, const char* txt, bool play){
  ensure_DISP();
  DISP.get_viewer() -> updateConfiguration(DISP);
  if(path.N) DISP.get_viewer() -> setMotion(DISP, path);
  DISP.get_viewer() -> view(pause, txt);
}

void PRM_PathFinder::ensure_DISP(){
  if(DISP.getJointStateDimension() != P->C.getJointStateDimension()){
    DISP.copy(P->C);
  }
}

} //namespace
