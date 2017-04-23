#include"SimCase.hpp"
//#include"general_tools.h"


SimCase::SimCase(void){
    emptypointer(grid_);
    emptypointer(grid_data_);
}

SimCase::~SimCase(void){
    emptypointer(grid_);
    emptypointer(grid_data_);
}

void SimCase::setup(const std::string &input_fname_){

    simdata_.Parse(input_fname_);

    logo();

    cout <<"\n-----------------------------------------------------------\n";
    cout <<"                   Checking Input File Parsing                 \n";
    cout <<"CaseTitle:   "<<simdata_.case_title<<endl;
    cout <<"Casemesh: "<<simdata_.mesh_fname<<endl;
    cout <<"Casepostprocessdir: "<<simdata_.case_postproc_dir<<endl;
    cout <<"CaseRestart: "<<simdata_.restart_flag<<endl;
    cout <<"CaseFvmOrder: "<<simdata_.scheme_order<<endl;
    cout <<"CaseRK: "<<simdata_.RK_order<<endl;
    cout <<"CaseReimann: "<<simdata_.ReimannSolver_type_<<endl;
    cout <<"Casedt_: "<<simdata_.dt_<<endl;
    cout <<"Caset_init_: "<<simdata_.t_init_<<endl;
    cout <<"Caset_end_: "<<simdata_.t_end_<<endl;
    cout <<"CasemaxIter: "<<simdata_.maxIter_<<endl;
    cout <<"\n-----------------------------------------------------------\n";


    allocator<char> allchar; // default allocator for char

    mkdir(simdata_.case_postproc_dir.c_str(),0777);

    char *current_working_dir=allchar.allocate(1000);
    getcwd(current_working_dir,1000);

    chdir(simdata_.case_postproc_dir.c_str());

    mkdir("./BINARY",0777);
    mkdir("./feild_output",0777);
    mkdir("./var_output",0777);
    mkdir("./conv_hist",0777);
    mkdir("./mesh_metrics",0777);

    chdir(current_working_dir);
    //chdir("..");

    cout<<"\n--> Currnet working directory: "<<current_working_dir<<endl;
    cout<<"--> Post processing directory: "<<simdata_.case_postproc_dir<<endl;


    return;
}

void SimCase::InitSim(){

    grid_ = new Mesh;

    grid_->Read(simdata_.mesh_fname);

    grid_->generate_meshData();

    std::string write_fname__;

    write_fname__="./post_process/gridtest.out";

    grid_->WriteMesh(write_fname__);

    grid_data_ = grid_->Release_meshData();

    emptypointer(grid_);

    return;
}


void SimCase::RunSim(){



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


