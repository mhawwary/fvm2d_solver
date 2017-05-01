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

    if(simdata.eqn_set=="Euler")
        fvm_space_solver = new Euler2DSolver;
    else if(simdata.eqn_set=="NavierStokes")
        fvm_space_solver = new NS2DSolver;
    else
        _notImplemented("An unidentified eqn_set");

    time_solver = new ExplicitTimeSolver;

    logo();

    allocator<char> allchar; // default allocator for char

    char eqn_type[100], case_dir[100];

    if(simdata.eqn_set=="Euler"){
        sprintf(eqn_type,"inviscid");
        sprintf(case_dir,"%s_%s_%1.1f_%1.1f",simdata.case_title.c_str(),eqn_type
                ,gasdata.alpha_deg_,gasdata.Mach_);
    }else{
        sprintf(eqn_type,"viscous");
        sprintf(case_dir,"%s_%s_%1.1f_%1.1f_%3.0f",simdata.case_title.c_str(),eqn_type
                ,gasdata.alpha_deg_,gasdata.Mach_,gasdata.Re_no);
    }

    mkdir(simdata.case_postproc_dir.c_str(),0777);  // parent directory if not already exists

    char *current_working_dir=allchar.allocate(1000);
    getcwd(current_working_dir,1000);

    chdir(simdata.case_postproc_dir.c_str());

    mkdir(case_dir,0777);

    simdata.case_postproc_dir+=case_dir;

    chdir(case_dir);

    mkdir("./BINARY",0777);
    mkdir("./field_output",0777);
    mkdir("./surface_output",0777);
    mkdir("./converg_hist",0777);
    mkdir("./mesh_metrics",0777);

    chdir(current_working_dir);

    copy_problem_inputdata();

    cout <<"\n-----------------------------------------------------------\n";
    cout <<"                   Checking Input File Parsing                 \n";
    cout <<"CaseTitle:   "<<case_dir<<endl;
    cout <<"Postprocessdir: "<<simdata.case_postproc_dir<<endl;

    cout <<"\nEquationSet: "<<simdata.eqn_set<<endl;
    cout <<"SchemeOrder: "<<simdata.scheme_order<<endl;
    cout <<"RK order: "<<simdata.RK_order<<endl;
    cout <<"Reimann Solver: "<<simdata.ReimannSolver_type_<<endl;
    if(simdata.use_local_timeStep==0){
        cout <<"fixed global time step dt: "<<simdata.dt_<<endl;
    }else{
        cout <<"local time stepping with global CFL: "<<simdata.CFL_<<endl;
    }
    cout <<"\n-----------------------------------------------------------\n";

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

    write_fname_=simdata.case_postproc_dir+"/mesh_metrics/grid_info.dat";
    write_fname_tec=simdata.case_postproc_dir+"/mesh_metrics/gridtest.plt";

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

    double Resid_sum_,cont_Resid_norm
            ,Xmom_Resid_norm,Ymom_Resid_norm
            ,Energy_Resid_norm;
    double **Qv=nullptr;

    int n=0;

    // Computing and dumping Initial Residual:
    //----------------------------------------------
    time_solver->ComputeInitialResid(fvm_space_solver->GetNumSolution());
    time_solver->SolveOneStep(fvm_space_solver->GetNumSolution());
    n=time_solver->GetIter();  // update iteration number

    cont_Resid_norm = time_solver->GetContinuityResNorm(2);
    Xmom_Resid_norm = time_solver->GetMomentumXResNorm(2);
    Ymom_Resid_norm = time_solver->GetMomentumYResNorm(2);
    Energy_Resid_norm = time_solver->GetEnergyResNorm(2);
    Resid_sum_ = time_solver->GetResNorm(2);

    printf("\nL2, Iter: %d, Rsum:%e, rho:%e, rhoU:%e, rhoV:%e, E:%e"
           ,n,Resid_sum_,cont_Resid_norm
           , Xmom_Resid_norm,Ymom_Resid_norm
           ,Energy_Resid_norm);

    dump_resid_L2norm(n,Resid_sum_,cont_Resid_norm
                    ,Xmom_Resid_norm,Ymom_Resid_norm
                    ,Energy_Resid_norm);

    // Main Solution Loop

    while (cont_Resid_norm>simdata.conv_threshold){

            time_solver->SolveOneStep(fvm_space_solver->GetNumSolution());

            n=time_solver->GetIter();

            if(n%simdata.conv_hist_pfreq==0){

                //Resid_norm_= time_solver->GetResNorm();
                cont_Resid_norm = time_solver->GetContinuityResNorm(2);
                Xmom_Resid_norm = time_solver->GetMomentumXResNorm(2);
                Ymom_Resid_norm = time_solver->GetMomentumYResNorm(2);
                Energy_Resid_norm = time_solver->GetEnergyResNorm(2);
                Resid_sum_ = time_solver->GetResNorm(2);

                printf("\nL2, Iter: %d, Rsum:%e, rho:%e, rhoU:%e, rhoV:%e, E:%e"
                       ,n,Resid_sum_,cont_Resid_norm
                       , Xmom_Resid_norm,Ymom_Resid_norm
                       ,Energy_Resid_norm);

                dump_resid_L2norm(n,Resid_sum_,cont_Resid_norm
                                ,Xmom_Resid_norm,Ymom_Resid_norm
                                ,Energy_Resid_norm);

                if(n%simdata.forces_print_freq==0){
                    fvm_space_solver->Compute_vertex_sol(n);
                    Qv = fvm_space_solver->GetVertexSolution();
                    dump_wall_data(Qv,n);
                }

                if(n%simdata.fields_print_freq==0){
                    fvm_space_solver->Compute_vertex_sol(n);
                    Qv = fvm_space_solver->GetVertexSolution();
                    dump_field_data(Qv,n);
                }

                if(simdata.use_local_timeStep==1 && n<=2000){
                    time_solver->update_local_timestep();
                    dt_min = time_solver->getdt_min();
                    dt_max = time_solver->getdt_max();
                    printf("\ndt_min: %e,  dt_max: %e",dt_min, dt_max);
                }
            }

            if(n>5e7) break;
        }

    // Printing the final Residual:
    //------------------------------------
    //Resid_norm_= time_solver->GetResNorm();
    cont_Resid_norm = time_solver->GetContinuityResNorm(2);
    Xmom_Resid_norm = time_solver->GetMomentumXResNorm(2);
    Ymom_Resid_norm = time_solver->GetMomentumYResNorm(2);
    Energy_Resid_norm = time_solver->GetEnergyResNorm(2);
    Resid_sum_ = time_solver->GetResNorm(2);

    printf("\nIter: %d, Rsum:%e, rho:%e, rhoU:%e, rhoV:%e, E:%e\n\n"
           ,n,Resid_sum_,cont_Resid_norm
           , Xmom_Resid_norm,Ymom_Resid_norm
           ,Energy_Resid_norm);

    dump_resid_L2norm(n,Resid_sum_,cont_Resid_norm
                    ,Xmom_Resid_norm,Ymom_Resid_norm
                    ,Energy_Resid_norm);

    fvm_space_solver->Compute_vertex_sol(n);
    Qv = fvm_space_solver->GetVertexSolution();
    dump_wall_data(Qv,n);
    dump_field_data(Qv,n);

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

void SimCase::dump_resid_L2norm(const int& iter_, double& resid_sum_
                              ,double& cont_resid,double& Xmom_resid
                              , double& Ymom_resid,double& Energy_resid){
    char* fname=nullptr;
    fname =new char[100];

    sprintf(fname,"%s/converg_hist/residual_hist.dat"
            ,simdata.case_postproc_dir.c_str());

    static int check=0;

    if(check==0){

        FILE*  outfile1=fopen(fname,"w");

        fprintf(outfile1,"%d %e %e %e %e %e\n",iter_,cont_resid
                    ,Xmom_resid,Ymom_resid,Energy_resid,resid_sum_);

        fclose(outfile1); check++;

    }else{

        FILE*  outfile=fopen(fname,"at+");

        fprintf(outfile,"%d %e %e %e %e %e\n",iter_,cont_resid
                ,Xmom_resid,Ymom_resid,Energy_resid,resid_sum_);

        fclose(outfile);
    }

    emptyarray(fname);

    return;
}

void SimCase::dump_wall_data(double **Qv, const int& oiter){

    register int i; int nID;

    char *fname=nullptr; fname=new char[100];

    sprintf(fname,"%s/surface_output/upper_wall_data_%d.dat"
            ,simdata.case_postproc_dir.c_str(),oiter);

    FILE*  outfile1=fopen(fname,"wt");

    for(i=0; i<grid_data->NupperWallnodes; i++){
        nID = grid_data->upper_wall_nodelist[i];

        fprintf(outfile1,"%e %e %e %e %e %e %e\n",grid_data->Xn[nID] ,grid_data->Yn[nID]
                , Qv[nID][0], Qv[nID][1], Qv[nID][2], Qv[nID][3], Qv[nID][4]);
    }

    fclose(outfile1);
    emptyarray(fname);

    fname=new char[100];

    sprintf(fname,"%s/surface_output/lower_wall_data_%d.dat"
            ,simdata.case_postproc_dir.c_str(),oiter);

    FILE*  outfile=fopen(fname,"wt");

    for(i=0; i<grid_data->NlowerWallnodes; i++){
        nID = grid_data->lower_wall_nodelist[i];

        fprintf(outfile,"%e %e %e %e %e %e %e\n",grid_data->Xn[nID] ,grid_data->Yn[nID]
                , Qv[nID][0], Qv[nID][1], Qv[nID][2], Qv[nID][3],  Qv[nID][4]);
    }

    fclose(outfile);
    emptyarray(fname);

    return;
}

void SimCase::dump_field_data(double **Qv, const int& oiter){

    register int i; int j;

    char *fname=nullptr; fname=new char[100];

    sprintf(fname,"%s/field_output/contour_data_%d.plt"
            ,simdata.case_postproc_dir.c_str(),oiter);

    FILE*  outfile=fopen(fname,"wt");

    // Printing the header first:
    //-------------------------------------
    fprintf(outfile, "VARIABLES = \"X\",\"Y\",\"RHO\",\"u\",\"v\",\"Cp\",\"M\"\n");

    fprintf(outfile, "ZONE N=%d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n"
            , grid_data->Nnodes_postproc
            , grid_data->NpostProc);

    // Printing the solution variables at the nodes:
    //------------------------------------------------
    for(i=0; i<grid_data->Nnodes_postproc; i++){

        fprintf(outfile,"%e %e %e %e %e %e %e\n"
                ,grid_data->post_Xn[i] ,grid_data->post_Yn[i]
                , Qv[i][0], Qv[i][1], Qv[i][2], Qv[i][3], Qv[i][4]);
    }

    fprintf(outfile,"\n");

    // printing connectivity information:
    //-------------------------------------
    int node_id;

    for(i=0; i<grid_data->NpostProc; i++)
    {
        for(j=0; j<grid_data->post_proc_elemlist[i].n_local_nodes; j++){
            node_id = grid_data->post_proc_elemlist[i].to_node[j];
            fprintf(outfile, "%d ",node_id+1);
        }

        if(grid_data->post_proc_elemlist[i].n_local_nodes==3)
            fprintf(outfile, "%d",node_id+1);
        fprintf(outfile,"\n");
    }

    fclose(outfile);
    emptyarray(fname);

    return;
}

void SimCase::copy_problem_inputdata(){

    char dirr[100];

    sprintf(dirr,"%s/case_input.in",simdata.case_postproc_dir.c_str());

    copyFile("./pre_process/case_input.in",dirr);

    return;
}

void SimCase::copyFile(const std::string& fileNameFrom, const std::string& fileNameTo)
{
     std::ifstream in (fileNameFrom.c_str());
     std::ofstream out (fileNameTo.c_str());
     out << in.rdbuf();
     out.close();
     in.close();
}
