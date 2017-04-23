
#include"Euler2DSolver.hpp"

// Constructor/Destructor/ Setup functions:
//------------------------------------------------
//Euler2DSolver::Euler2DSolver(void){

//}

Euler2DSolver::~Euler2DSolver(void){
    Reset_solver();
}

void Euler2DSolver::setup_solver(MeshData*& meshdata_, SimData& osimdata_
                                 , GasProb& ogasdata_){
    grid_ = meshdata_;
    simdata_ = &osimdata_;
    gasdata_ = &ogasdata_;

    Nelem = grid_->Nelem;
    Nfaces = grid_->Nfaces;
    Nnodes = grid_->Nnodes;
    NbndElem = grid_->NbndElem;
    Nbndfaces = grid_->Nbndfaces;

    register int i;

    Qn    =  new double* [Nelem];

    for(i=0; i<grid_->Nelem; i++)
        Qn[i]    = new double[Ndof];

    Qv = new double*[Nnodes];

    for(i=0; i<Nnodes; i++)
        Qv[i]    = new double[Ndof];

    flux_com = new double*[Nfaces];

    for(i=0; i<Nfaces; i++)
        flux_com[i]    = new double[Ndof];

    SetPhyTime(simdata_->t_init_);  // initial physical time

    dt_ = simdata_->dt_;

    CalcTimeStep();

    // Screen Output of input and simulation parameters:
//    cout <<"\n===============================================\n";
//    cout << "CFL no.: "<<CFL<<endl;
//    cout << "dx     : "<<grid_->dx<<endl;
//    cout << "dt     : "<<dt_<<endl;
//    cout << "MaxIter: "<<simdata_->t_end_/time_step<<endl;
//    cout << "t_end  : "<<simdata_->t_end_<<endl;
//    cout << "\nNumber of Elements: "<< grid_->Nelem<<endl;
//    cout << "Polynomial  order : "<< simdata_->poly_order_  << endl;
//    cout << "Runge-Kutta order : "<< simdata_->RK_order_    << endl;
//    cout << "Upwind parameter  : "<< simdata_->upwind_param_<< endl <<"\n";

    return;
}

void Euler2DSolver::CalcTimeStep(){

//    if(simdata_->calc_dt_flag==1){

//        dt_ = (grid_->dx * simdata_->CFL_ )/ simdata_->a_wave_;
//        CFL = simdata_->CFL_;

//        T_period = (grid_->xf - grid_->x0) / simdata_->a_wave_;

//    }else if(simdata_->calc_dt_flag==0){

//        time_step = simdata_->dt_;

//        CFL = simdata_->a_wave_ * time_step / grid_->dx ;

//    }else {

//        FatalError("Wrong Calc_dt_flag");
//    }

    return;
}

void Euler2DSolver::Reset_solver(){

    //FatalError("Inside Euler2DSolver Destructor");

    emptyarray(Nelem,Qn);

    //FatalError("After empty Qn");

    emptyarray(Nnodes,Qv);

    //FatalError("After empty Qv");

    emptyarray(Nfaces,flux_com);

    return;
}

void Euler2DSolver::InitSol(){

    register int i;

    double rho,u,v,E;

    gasdata_->CalculateFlowProperties(rho,u,v,E);

    for(i=0; i<Nelem; i++){

        Qn[i][0] = rho;
        Qn[i][1] = rho*u;
        Qn[i][2] = rho*v;
        Qn[i][3] = E;
    }


    return;
}
