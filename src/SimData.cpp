#include"SimData.hpp"

void SimData::Parse(const std::string &fname){


    GetPot gp_input(fname.c_str());

    case_title = gp_input("Case/title","default");

    case_postproc_dir = gp_input("Case/postprocess_dir","default");

    mesh_fname = gp_input("Case/mesh_file_name","default");

    restart_flag=gp_input("Simulation/restart_flag",0);
    restart_iter=gp_input("Simulation/restart_iter",0);
    forces_print_freq=gp_input("Simulation/forces_print_freq",0);
    fields_print_freq=gp_input("Simulation/fields_print_freq",0);

    scheme_order=gp_input("space_solver/scheme_order",1);
    ReimannSolver_type_ = gp_input("space_solver/Riemann_solver","Rusanov");
    eqn_set = gp_input("space_solver/eqn_set","Euler");
    FarFieldBC = gp_input("space_solver/FarField_Boundary_condition","Characteristics");
    WallBC = gp_input("space_solver/Wall_Boundary_condition","Slip");

    dt_ = gp_input("time_solver/dt",1e-9);
    t_init_ = gp_input("time_solver/initial_time",1e-9);
    t_end_ = gp_input("time_solver/final_time",1e-9);
    maxIter_ = gp_input("time_solver/maximum_iteration",1e-9);
    CFL_ = gp_input("time_solver/CFL",0.5);
    RK_order=gp_input("time_solver/explicit/RK_order",0);

    return;
}

void GasProb::Parse(const std::string &fname){

    GetPot gp_input(fname.c_str());

    gama = gp_input("Fluid/gamma",1.4);
    R_gas = gp_input("Fluid/R_gas_constant",287);
    Prndtl = gp_input("Fluid/Prandtl_no",0.72);
    Ptot = gp_input("Fluid/P_infinity",100000.0);
    Rho_tot = gp_input("Fluid/rho_infinity",1.225);
    T_tot = gp_input("Fluid/T_inifinity",300);


    alpha_deg_ = gp_input("Flow/alpha",0.);
    Mach_ = gp_input("Flow/Mach",0.5);
    Re_no = gp_input("Flow/Reynolds_no",0.);

    return;
}

void GasProb::CalculateFlowProperties(double& rho_,
                                      double& u_, double& v_,
                                      double& p_, double& T_,
                                      double& E_){

    double c_=0.0,g_gm1=0.0, gm1_2=0.0;

    g_gm1 = gama/(gama-1.);

    gm1_2 = (gama-1.)/2.;

    p_ = Ptot / pow((1.+ (gm1_2 * Mach_* Mach_)),g_gm1);
    T_ = T_tot/(1.+ (gm1_2 * Mach_* Mach_));

    rho_ = p_/(R_gas * T_);

    c_ = sqrt(gama*p_/rho_);

    u_ = Mach_*c_*cos(alpha_deg_ * PI /180.0);
    v_ = Mach_*c_*sin(alpha_deg_ * PI /180.0);

    E_ = (p_/(gama-1.0)) + 0.5 * rho_ * (pow(u_,2)+pow(v_,2)) ;

    return;
}


void GasProb::CalculateFlowProperties(double& rho_,
                                      double& u_, double& v_,
                                      double& E_){

    double c_=0.0,g_gm1=0.0, gm1_2=0.0, p_=0.0,T_=0.0;

    g_gm1 = gama/(gama-1.);

    gm1_2 = (gama-1.)/2.;

    p_ = Ptot / pow((1.+ (gm1_2 * Mach_* Mach_)),g_gm1);
    T_ = T_tot/(1.+ (gm1_2 * Mach_* Mach_));

    rho_ = p_/(R_gas * T_);

    c_ = sqrt(gama*p_/rho_);

    u_ = Mach_*c_*cos(alpha_deg_ * PI /180.0);
    v_ = Mach_*c_*sin(alpha_deg_ * PI /180.0);

    E_ = (p_/(gama-1.0)) + 0.5 * rho_ * (pow(u_,2)+pow(v_,2)) ;

    return;
}
