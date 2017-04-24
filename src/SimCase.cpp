#include"SimCase.hpp"
//#include"general_tools.h"


SimCase::SimCase(void){
    emptypointer(grid);
    emptypointer(grid_data);
}

SimCase::~SimCase(void){
    //FatalError("Inside SimCase Destructor");

    emptypointer(fvm_space_solver);

    //FatalError("After empty fvm_space solver");

    emptypointer(grid);

    //FatalError("After empty grid");

    //printf("\n*grid_data: %d\n",grid_data);

    emptypointer(grid_data);

    //printf("\n*grid_data: %d\n",grid_data);

    emptypointer(time_solver);
}

void SimCase::setup(const std::string &input_fname_){

    simdata.Parse(input_fname_);
    gasdata.Parse(input_fname_);

    fvm_space_solver = new Euler2DSolver;

    time_solver = new ExplicitTimeSolver;

    logo();

    cout <<"\n-----------------------------------------------------------\n";
    cout <<"                   Checking Input File Parsing                 \n";
    cout <<"CaseTitle:   "<<simdata.case_title<<endl;
    cout <<"Casemesh: "<<simdata.mesh_fname<<endl;
    cout <<"Casepostprocessdir: "<<simdata.case_postproc_dir<<endl;
    cout <<"CaseRestart: "<<simdata.restart_flag<<endl;
    cout <<"CaseFvmOrder: "<<simdata.scheme_order<<endl;
    cout <<"CaseRK: "<<simdata.RK_order<<endl;
    cout <<"CaseReimann: "<<simdata.ReimannSolver_type_<<endl;
    cout <<"Casedt_: "<<simdata.dt_<<endl;
    cout <<"Caset_init_: "<<simdata.t_init_<<endl;
    cout <<"Caset_end_: "<<simdata.t_end_<<endl;
    cout <<"CasemaxIter: "<<simdata.maxIter_<<endl;
    cout <<"\n-----------------------------------------------------------\n";


    allocator<char> allchar; // default allocator for char

    mkdir(simdata.case_postproc_dir.c_str(),0777);

    char *current_working_dir=allchar.allocate(1000);
    getcwd(current_working_dir,1000);

    chdir(simdata.case_postproc_dir.c_str());

    mkdir("./BINARY",0777);
    mkdir("./feild_output",0777);
    mkdir("./var_output",0777);
    mkdir("./conv_hist",0777);
    mkdir("./mesh_metrics",0777);

    chdir(current_working_dir);
    //chdir("..");

    cout<<"\n--> Currnet working directory: "<<current_working_dir<<endl;
    cout<<"--> Post processing directory: "<<simdata.case_postproc_dir<<endl;

    return;
}

void SimCase::InitSim(){

    // Preparing MeshData:
    //------------------------

    grid = new Mesh;

    grid->Read(simdata.mesh_fname);

    grid->generate_meshData();

    std::string write_fname_;

    write_fname_="./post_process/gridtest.out";

    grid->WriteMesh(write_fname_);

    grid_data = grid->Release_meshData();

    emptypointer(grid);

    // Setupping Space and Time Solvers and initializing the solution:
    //-----------------------------------------------------------------

    fvm_space_solver->setup_solver(grid_data,simdata,gasdata);

    fvm_space_solver->InitSol();

    time_solver->setupTimeSolver(fvm_space_solver,&simdata);

    return;
}


void SimCase::RunSim(){

    //double gtime=fvm_space_solver->GetPhyTime();

    //double dt_= fvm_space_solver->GetTimeStep();

    double Resid_norm_=10.0;

    int n=0;

    while (Resid_norm_>1e-15){

            time_solver->SolveOneStep(fvm_space_solver->GetNumSolution());

            //time_solver->space_solver_->UpdatePhyTime(dt_);

            //gtime=fvm_space_solver->GetPhyTime();

            Resid_norm_= time_solver->GetResNorm();

             if(n%10000==1) _compare(n,Resid_norm_);

            n++;

            if(n>1e7) break;
        }

    _compare(n,Resid_norm_);

    return;
}

void SimCase::PostProcess(){

    return;
}

void SimCase::logo() {

    cout<<"                                                                                         "<<endl;
    cout<<"_________________________________________________________________________________________"<<endl;
    cout<<"                                                                                         "<<endl;
    cout<<"                                                                                         "<<endl;
    cout<<"         "<<"Welcome to the High-Order Finite Volume Fluid Simulation Platform"<<"       "<<endl;
    cout<<"                                                                                         "<<endl;
    cout<<"                                                                                         "<<endl;
    cout<<"         Author:               Mohammad Alhawwary, PhD. Student                          "<<endl;
    cout<<"    Affiliation:   Aerospace Engineering Department, University of Kansas, USA           "<< endl;
    cout<<"                                                                                         "<<endl;
    cout<<"                                                                                         "<<endl;
    cout<<"_________________________________________________________________________________________"<<endl;

    return ;
}


