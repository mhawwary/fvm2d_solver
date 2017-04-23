#ifndef SIMCASE_H
#define SIMCASE_H

#include"Mesh.hpp"
#include"general_tools.h"
#include"SimData.hpp"
#include"SpaceSolver.hpp"
#include"Euler2DSolver.hpp"
#include"NS2DSolver.hpp"
#include"ExplicitTimeSolver.hpp"

using namespace std;

class SimCase {

public:
//  Construction Functions :
  SimCase(void);
  ~SimCase(void);
  void setup(const std::string &input_fname_);
  void InitSim();
  void RunSim();
  void PostProcess();

protected:
  void logo();

protected:
  //std::string input_fname;  // input file name
  SimData simdata ;
  GasProb gasdata ;
  Mesh     *grid=nullptr;
  MeshData *grid_data=nullptr;
  SpaceSolver *fvm_space_solver=nullptr;
  ExplicitTimeSolver *time_solver=nullptr;

};

#endif
