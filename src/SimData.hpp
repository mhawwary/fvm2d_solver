#ifndef SIMDATA_H
#define SIMDATA_H

#include"../include/getpot.h"
#include"global_var.h"
#include"general_tools.h"

struct GasProb{

    // Gas Properties:
    //------------------

    // Total Conditions:
    double gama=1.4;
    double R_gas = 287.0;
    double Prndtl = 0.72;
    double Ptot = 101325.0;
    double Rho_tot = 1.225;
    double T_tot = 288.2030861;

    // Static Conditions:
    double Ps=Ptot;
    double Ts=T_tot;
    double Rho_s= Rho_tot;

    // Flow Properties:
    //--------------------
    double alpha_deg_=0.0;
    double Mach_ = 0.5;
    double Re_no = 1.0e6;

    void Parse(const std::string &fname);
    void compute_static_conditions();

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
    std::string FarFieldBC;
    std::string WallBC;

    int restart_flag=0;  //0: start a new simulation; 1: restart a previous simulation
    int restart_iter=0;
    int forces_print_freq=10;
    int fields_print_freq=10;
    int conv_hist_pfreq=10;
    double conv_threshold=1e-15;

    double dt_ = 1e-3;  // dt time step
    double t_init_ = 0;  // initial time
    double t_end_ =1e20;  // end time
    double maxIter_ = 1e10; // maximum number of iterations
    double CFL_=0.5;
    int RK_order=0;       // Runge-Kutta type (0: euler FT, 2: SSPRK22, 3: SSPRK33)
    int use_local_timeStep=0;

    void Parse(const std::string &fname);

};

#endif
