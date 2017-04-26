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
    cout <<"MeshFile: "<<simdata.mesh_fname<<endl;
    cout <<"Postprocessdir: "<<simdata.case_postproc_dir<<endl;
    cout <<"SchemeOrder: "<<simdata.scheme_order<<endl;
    cout <<"RK order: "<<simdata.RK_order<<endl;
    cout <<"Reimann Solver: "<<simdata.ReimannSolver_type_<<endl;
    cout <<"Fixed_dt_: "<<simdata.dt_<<endl;
    cout <<"CFL: "<<simdata.CFL_<<endl;
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

    std::string write_fname_tec;

    write_fname_="./post_process/gridtest.out";
    write_fname_tec="./post_process/gridtest.tp";

    grid->WriteMesh(write_fname_);

    grid->WriteMeshTecplot(write_fname_tec);

    grid_data = grid->Release_meshData();

    emptypointer(grid);

    // Setupping Space and Time Solvers and initializing the solution:
    //-----------------------------------------------------------------

    fvm_space_solver->setup_solver(grid_data,simdata,gasdata);

    fvm_space_solver->InitSol();

    time_solver->setupTimeSolver(fvm_space_solver,&simdata, grid_data);

    return;
}

void SimCase::RunSim(){

    //double gtime=fvm_space_solver->GetPhyTime();

    double dt_min,dt_max;

    double Resid_norm_=10.0, cont_Resid_norm,Xmom_Resid_norm,Ymom_Resid_norm,Energy_Resid_norm;
    double **Qv=nullptr;

    int n=0;

    // Computing and dumping Initial Residual:
    //----------------------------------------------
    time_solver->ComputeInitialResid(fvm_space_solver->GetNumSolution());
    time_solver->SolveOneStep(fvm_space_solver->GetNumSolution());
    n=time_solver->GetIter();

    Resid_norm_= time_solver->GetResNorm();
    cont_Resid_norm = time_solver->GetContinuityResNorm();
    Xmom_Resid_norm = time_solver->GetMomentumXResNorm();
    Ymom_Resid_norm = time_solver->GetMomentumYResNorm();
    Energy_Resid_norm = time_solver->GetEnergyResNorm();

    printf("\nIter: %d,  Resid_sum:%e , rho:%e, rhoU:%e, rhoV:%e, E:%e"
           ,n,Resid_norm_,cont_Resid_norm
           , Xmom_Resid_norm,Ymom_Resid_norm
           ,Energy_Resid_norm);

    dump_resid_norm(n,Resid_norm_,cont_Resid_norm
                    ,Xmom_Resid_norm,Ymom_Resid_norm
                    ,Energy_Resid_norm);

    // Main Solution Loop

    while (n<1e7){

            time_solver->SolveOneStep(fvm_space_solver->GetNumSolution());

            n=time_solver->GetIter();

            if(n%simdata.conv_hist_pfreq==0){

                Resid_norm_= time_solver->GetResNorm();
                cont_Resid_norm = time_solver->GetContinuityResNorm();
                Xmom_Resid_norm = time_solver->GetMomentumXResNorm();
                Ymom_Resid_norm = time_solver->GetMomentumYResNorm();
                Energy_Resid_norm = time_solver->GetEnergyResNorm();

                printf("\nIter: %d,  Resid_sum:%e , rho:%e, rhoU:%e, rhoV:%e, E:%e"
                       ,n,Resid_norm_,cont_Resid_norm
                       , Xmom_Resid_norm,Ymom_Resid_norm
                       ,Energy_Resid_norm);

                dump_resid_norm(n,Resid_norm_,cont_Resid_norm
                                ,Xmom_Resid_norm,Ymom_Resid_norm
                                ,Energy_Resid_norm);

                if(n%simdata.forces_print_freq==0){
                    fvm_space_solver->Compute_vertex_sol();
                    Qv = fvm_space_solver->GetVertexSolution();
                    dump_wall_data(Qv);
                }

                if(simdata.use_local_timeStep==1 && n<=5000){
                    time_solver->update_local_timestep();
                    dt_min = time_solver->getdt_min();
                    dt_max = time_solver->getdt_max();
                    printf("\ndt_min: %e,  dt_max: %e",dt_min, dt_max);
                }
            }

//            if(n%100==1) {
//                printf("\nIter: %d,  Resid:%e , ResidC:%e, ResidX:%e, ResidY:%e, ResidEnergy:%e"
//                       ,n-1,Resid_norm_,cont_Resid_norm, Xmom_Resid_norm,Ymom_Resid_norm,Energy_Resid_norm);
//                if(simdata.use_local_timeStep==1){
//                dt_min = time_solver->getdt_min();
//                dt_max = time_solver->getdt_max();
//                printf("\nMinimum dt: %e,  Maximum dt: %e",dt_min, dt_max);
//                }
//            }

            if(n>5e7) break;
        }

    // Printing the final Residual:
    //------------------------------------
    Resid_norm_= time_solver->GetResNorm();
    cont_Resid_norm = time_solver->GetContinuityResNorm();
    Xmom_Resid_norm = time_solver->GetMomentumXResNorm();
    Ymom_Resid_norm = time_solver->GetMomentumYResNorm();
    Energy_Resid_norm = time_solver->GetEnergyResNorm();

    printf("\nIter: %d,  Resid_sum:%e , rho:%e, rhoU:%e, rhoV:%e, E:%e"
           ,n,Resid_norm_,cont_Resid_norm
           , Xmom_Resid_norm,Ymom_Resid_norm
           ,Energy_Resid_norm);

    dump_resid_norm(n,Resid_norm_,cont_Resid_norm
                    ,Xmom_Resid_norm,Ymom_Resid_norm
                    ,Energy_Resid_norm);

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

void SimCase::dump_resid_norm(const int& iter_, double& resid_norm, double& cont_resid
                              ,double& Xmom_resid, double& Ymom_resid, double& Energy_resid){

    char* fname=nullptr;
    fname =new char[100];

    sprintf(fname,"%sresidual_hist.dat",simdata.case_postproc_dir.c_str());

    FILE*  outfile=fopen(fname,"at+");

    fprintf(outfile,"%d %e %e %e %e %e\n",iter_,resid_norm,cont_resid
            ,Xmom_resid,Ymom_resid,Energy_resid);

    fclose(outfile);

    emptyarray(fname);

    return;
}

void SimCase::dump_wall_data(double **Qv){

    register int i; int nID;
    double c,M;

    //const char* fname = write_fname_.c_str();

    char *fname=nullptr; fname=new char[100];

    //fname = write_fname_.c_str();

    sprintf(fname,"%supper_wall_data.dat",simdata.case_postproc_dir.c_str());

    FILE*  outfile1=fopen(fname,"wt");

    fprintf(outfile1, "VARIABLES = \"X\",\"Y\",\"RHO\",\"u\",\"v\",\"P\",\"M\"");

    for(i=0; i<grid_data->NupperWallnodes; i++){
        nID = grid_data->upper_wall_nodelist[i];

        c= sqrt(gasdata.gama*Qv[nID][3]/Qv[nID][0]);
        M = sqrt( pow(Qv[nID][1],2)+pow(Qv[nID][2],2) ) / c;

        fprintf(outfile1,"\n%e %e %e %e %e %e %e",grid_data->Xn[nID] ,grid_data->Yn[nID]
                , Qv[nID][0], Qv[nID][1], Qv[nID][2], Qv[nID][3], M);
    }

    fclose(outfile1);
    emptyarray(fname);

    fname=new char[100];

    sprintf(fname,"%slower_wall_data.dat",simdata.case_postproc_dir.c_str());

    FILE*  outfile=fopen(fname,"wt");

    fprintf(outfile, "VARIABLES = \"X\",\"Y\",\"RHO\",\"u\",\"v\",\"P\",\"M\"");

    for(i=0; i<grid_data->NlowerWallnodes; i++){
        nID = grid_data->lower_wall_nodelist[i];

        c= sqrt(gasdata.gama*Qv[nID][3]/Qv[nID][0]);
        M = sqrt( pow(Qv[nID][1],2)+pow(Qv[nID][2],2) ) / c;

        fprintf(outfile,"\n%e %e %e %e %e %e %e",grid_data->Xn[nID] ,grid_data->Yn[nID]
                , Qv[nID][0], Qv[nID][1], Qv[nID][2], Qv[nID][3], M);
    }

    fclose(outfile);
    emptyarray(fname);

    return;
}
