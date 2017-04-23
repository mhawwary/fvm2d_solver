#ifndef SIMDATA_H
#define SIMDATA_H

#include"../include/getpot.h"
#include"global_var.h"

struct GasProb{

    // Gas Properties:
    //------------------

    double gama=1.4;
    double R_gas = 287;
    double Prndtl = 0.72;
    double Ptot = 101325;
    double Rho_tot = 1.225;
    double T_tot = 288.2030861;

    // Flow Properties:
    //--------------------
    double alpha_deg_=0.0;
    double Mach_ = 0.5;
    double Re_no = 1.0e6;

    void Parse(const std::string &fname);

    void CalculateFlowProperties(double& rho_,
                                 double& u_, double& v_,
                                 double& p_, double& T_,
                                 double& E_);

    void CalculateFlowProperties(double& rho_,
                                 double& u_, double& v_,
                                 double& E_);

};

struct SimData {

    std::string case_title;
    std::string mesh_fname;   // mesh file name
    std::string case_postproc_dir;  // postprocessing directory

    int scheme_order=0;   // FVM Scheme order
    std::string ReimannSolver_type_;
    std::string eqn_set;

    int restart_flag=0;  //0: start a new simulation; 1: restart a previous simulation
    int restart_iter=0;
    int forces_print_freq=10;
    int fields_print_freq=10;

    double dt_ = 1e-3;  // dt time step
    double t_init_ = 0;  // initial time
    double t_end_ =1e20;  // end time
    double maxIter_ = 1e10; // maximum number of iterations
    double CFL_=0.5;
    int RK_order=0;       // Runge-Kutta type (0: euler FT, 2: SSPRK22, 3: SSPRK33)

    void Parse(const std::string &fname);

};

#endif
