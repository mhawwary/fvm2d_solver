#include"SimData.hpp"

void SimData::Parse(const std::string &fname){


    GetPot gp_input(fname.c_str());

    case_title = gp_input("Case/title","default");

    case_postproc_dir = gp_input("Case/postprocess_dir","default");

    mesh_fname = gp_input("Case/mesh_file_name","default");

    restart_flag=gp_input("Simulation/restart_flag",0);

    print_freq=gp_input("Simulation/print_freq",0);

    scheme_order=gp_input("space_solver/order",1);

    ReimannSolver_type_ = gp_input("space_solver/Riemann_solver","upwind");

    RK_order=gp_input("time_solver/explicit/RK_order",0);

    dt_ = gp_input("time_solver/dt",1e-9);

    t_init_ = gp_input("time_solver/initial_time",1e-9);

    t_end_ = gp_input("time_solver/final_time",1e-9);

    maxIter_ = gp_input("time_solver/maximum_iteration",1e-9);

}
