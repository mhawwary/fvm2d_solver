#ifndef SIMCASE_H
#define SIMCASE_H

#include"Mesh.hpp"
//#include"MeshData.h"
#include"general_tools.h"
#include"SimData.hpp"

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
  std::string input_fname;  // input file name
  SimData simdata_ ;
  Mesh     *grid_=NULL;
  MeshData *grid_data_=NULL;
  //int simIter = 0;
  //double sim_phy_time=0.0;
  //double sim_dt=0.0;


};

#endif
