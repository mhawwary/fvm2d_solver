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
  void PostProcess(const int& iter_,  double& Resid_sum_);

protected:
  void logo();
  void dump_resid_L2norm(const int& iter_, double& resid_sum_
                         , double& cont_resid,double& Xmom_resid
                         , double& Ymom_resid,double& Energy_resid);
//  void dump_resid_L1norm(const int& iter_, double& resid_sum_
//                         , double& cont_resid,double& Xmom_resid
//                         , double& Ymom_resid,double& Energy_resid);
  void dump_wall_data(const int& oiter, double** Qv
                      ,double* p_wall_, double* tau_xx_wall_
                      , double* tau_xy_wall_, double* tau_yy_wall_
                      , double* Cf_);

  void dump_field_data(double **Qv, const int& oiter);
  void copy_problem_inputdata();
  void copyFile(const std::string& fileNameFrom, const std::string& fileNameTo);


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
