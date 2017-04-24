
#include"Euler2DSolver.hpp"

// Constructor/Destructor/ Setup functions:
//------------------------------------------------

Euler2DSolver::~Euler2DSolver(void){
    Reset_solver();
}

void Euler2DSolver::setup_solver(MeshData*& meshdata_, SimData& osimdata_
                                 , GasProb& ogasdata_){
    grid_ = meshdata_;
    simdata_ = &osimdata_;
    gasdata_ = &ogasdata_;

    Nelem = grid_->Nelem;
    Nelem_extend = grid_->Nelem_extend;
    Nfaces = grid_->Nfaces;
    Nnodes = grid_->Nnodes;
    NbndElem = grid_->NbndElem;
    Nbndfaces = grid_->Nbndfaces;

    scheme_order = simdata_->scheme_order;

    register int i;

    Qc  =  new double* [Nelem_extend];

    for(i=0; i<Nelem_extend; i++)
        Qc[i]    = new double[Ndof];

    if(scheme_order==2){
        dQdx = new double*[Nelem];
        dQdy = new double*[Nelem];

        for(i=0; i<Nelem; i++){
            dQdx[i]  = new double[Ndof];
            dQdy[i]  = new double[Ndof];
        }
    }else if(scheme_order!=1 && scheme_order!=2){
        _notImplemented("Error in parsing scheme order, scheme order is not implemented yet");
    }

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

    emptyarray(Nelem_extend,Qc);
    emptyarray(Nelem,dQdx);
    emptyarray(Nelem,dQdy);

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

        Qc[i][0] = rho;
        Qc[i][1] = rho*u;
        Qc[i][2] = rho*v;
        Qc[i][3] = E;
    }

    return;
}

void Euler2DSolver::UpdateResid(double **Resid_, double **Qc_){

    //FatalError("Inside Update Resid in Euler Solver");

    SetGhostVariables();

    //FatalError("Finished setting Ghost Variables");

    if(scheme_order==2) Reconstruct_sol();

    //FatalError("Finished Reconstruction");

    register int i,j;

    // Boundary Face loop to calculate the common interface fluxes:
    //----------------------------------------------------
    double Ql[4]={0.,0.,0.,0.}, Qr[4]={0.,0.,0.,0.};
    double nx=0.0,ny=0.0;

    // Boundary faces:

    for(j=0; j<Nbndfaces; j++){

        Compute_left_right_boundfacesol(j,&Ql[0],&Qr[0]);

        Compute_common_inviscidflux(&Ql[0],&Qr[0],nx,ny, &flux_com[j][0]);
    }

    // Interior faces
    for(j=Nbndfaces; j<Nfaces; j++){

        Compute_left_right_facesol(j,&Ql[0],&Qr[0]);

        Compute_common_inviscidflux(&Ql[0],&Qr[0],nx,ny, &flux_com[j][0]);
    }

    /* Face loop to update the Residuals in each respective cell
     * Two cells are updated at a time
     */

    // Set Residuals to zero before computing:
    Residulas_setZero(Resid_);

    int iL,iR;
    double Vol_l,Vol_r;

    double resid_temp[4]={0.,0.,0.,0.};

    // Note that there should be a negative sign added to all Residuals;

    // Bound Faces
    for(i=0; i<Nbndfaces; i++){

        iL = grid_->facelist[i].Lcell;
        Vol_l = grid_->elemlist[iL].Vc;

        compute_resid_OneFace(i, &resid_temp[0]);

        for(j=0; j<Ndof; j++){
            resid_temp[j] = (-resid_temp[j]/Vol_l); // to move it to the RHS of the equation

            Resid_[iL][j] += resid_temp[j];
        }
    }

    // Interior Faces:

    for(i=Nbndfaces; i<Nfaces; i++){

        iL = grid_->facelist[i].Lcell;
        iR = grid_->facelist[i].Rcell;

        Vol_l = grid_->elemlist[iL].Vc;
        Vol_r = grid_->elemlist[iR].Vc;

        if(Vol_l==0 || Vol_r==0)
            _compare(Vol_l,Vol_r);

        compute_resid_OneFace(i, &resid_temp[0]);

        for(j=0; j<Ndof; j++){
            resid_temp[j] = -resid_temp[j]; // to move it to the RHS of the equation

            Resid_[iL][j] +=   (resid_temp[j]/Vol_l);
            Resid_[iR][j] +=  (-resid_temp[j]/Vol_r);
        }
    }

    return;
}

void Euler2DSolver::SetGhostVariables(){

    register int i;

    int iL,iR;

    double nx=0.0,ny=0.0;

    for(i=Nelem; i<Nelem_extend; i++){

        iL = grid_->facelist[i-Nelem].Lcell;
        iR = grid_->facelist[i-Nelem].Rcell;
        nx = grid_->facelist[i-Nelem].nx;
        ny = grid_->facelist[i-Nelem].ny;

        if(simdata_->WallBC=="Extrapolation"){
            int j=0;
            for(j=0; j<Ndof; j++) Qc[iR][j] = Qc[iL][j];

        }else{

            if(grid_->elemlist[i].bnd_type==-1){ // wall B.C

                /*if(iR!=i)
            _compare(iR,i);

            if(grid_->facelist[i-Nelem].bnd_type!=-1)
                _print("face bndtype is wrong for inviscid wall"); */

                compute_GhostSol_inviscidWallBC(nx,ny,&Qc[iL][0],&Qc[iR][0]);

            }else if(grid_->elemlist[i].bnd_type==-2){ // FarField B.C

                /*if(iR!=i)
            _compare(iR,i);

            if(grid_->facelist[i-Nelem].bnd_type!=-2)
                _print("face bndtype is wrong for farfield");*/

                compute_GhostSol_farfieldBC(nx,ny,&Qc[iL][0],&Qc[iR][0]);

            }else{
                FatalError_exit("This ghost cell boundary type is wrong");
            }
        }
    }

    return;
}

void Euler2DSolver::compute_GhostSol_inviscidWallBC(const double& nx, const double& ny
                                                    ,const double *Ql, double *Qr){

    double ul=0.0,vl=0.0,Vnl,ur=0.0,vr=0.0;

    ul = Ql[1]/Ql[0];
    vl = Ql[2]/Ql[0];

    Vnl = ul*nx + vl*ny; // normal component of left cell velocity vector

    ur = ul - 2*Vnl * nx;
    vr = vl - 2*Vnl * ny;

    Qr[0] = Ql[0];
    Qr[1] = Qr[0] * ur;
    Qr[2] = Qr[0] * vr;
    Qr[3] = Ql[3];   // since they are of same velocity magnitude and Pr=Pl

    return;
}

void Euler2DSolver::compute_GhostSol_farfieldBC(const double& nx, const double& ny
                                                ,const double *Ql, double *Qr){

    if(simdata_->FarFieldBC=="Fixed"){

        double rho,u,v,E;
        gasdata_->CalculateFlowProperties(rho,u,v,E);

        Qr[0] = rho;
        Qr[1] = rho*u;
        Qr[2] = rho*v;
        Qr[3]=E;

    }else if(simdata_->FarFieldBC=="Characteristics"){

    }else{
        _notImplemented("This FarField Boundary Condition is not implemented");
    }

    return;
}

void Euler2DSolver::Compute_left_right_facesol(const int& fID
                                               , double *Ql_, double *Qr_){

    int j=0,iL,iR;

    double Xf,Yf,Xcl,Ycl,Xcr,Ycr;

    iL = grid_->facelist[fID].Lcell;
    iR = grid_->facelist[fID].Rcell;
    Xf = grid_->facelist[fID].Xf;
    Yf = grid_->facelist[fID].Yf;

    Xcl = grid_->elemlist[iL].Xc;
    Ycl = grid_->elemlist[iL].Yc;

    Xcr = grid_->elemlist[iR].Xc;
    Ycr = grid_->elemlist[iR].Yc;

    if(scheme_order==1){
        for(j=0; j<Ndof; j++){
            Ql_[j] = Qc[iL][j];
            Qr_[j] = Qc[iR][j];
        }

    }else if(scheme_order==2){
        for(j=0; j<Ndof; j++){
            Ql_[j] = Qc[iL][j] + dQdx[iL][j] * (Xf-Xcl) + dQdy[iL][j] * (Yf-Ycl);
            Qr_[j] = Qc[iR][j] + dQdx[iR][j] * (Xf-Xcr) + dQdy[iR][j] * (Yf-Ycr);
        }
    }

    return;
}

void Euler2DSolver::Compute_left_right_boundfacesol(const int& fID
                                               , double *Ql_, double *Qr_){

    int j=0,iL=0;

    double Xf,Yf,Xcl,Ycl,nx,ny;

    iL = grid_->facelist[fID].Lcell;
    Xf = grid_->facelist[fID].Xf;
    Yf = grid_->facelist[fID].Yf;
    nx = grid_->facelist[fID].nx;
    nx = grid_->facelist[fID].ny;

    Xcl = grid_->elemlist[iL].Xc;
    Ycl = grid_->elemlist[iL].Yc;

    if(scheme_order==1){
        for(j=0; j<Ndof; j++)
            Ql_[j] = Qc[iL][j];

    }else if(scheme_order==2){
        for(j=0; j<Ndof; j++)
            Ql_[j] = Qc[iL][j] + dQdx[iL][j] * (Xf-Xcl) + dQdy[iL][j] * (Yf-Ycl);
    }

    if(simdata_->WallBC=="Extrapolation"){
        int j=0;
        for(j=0; j<Ndof; j++) Qr_[j] = Ql_[j];

    }else{
        if(grid_->facelist[fID].bnd_type==-1){ // wall BC

            compute_GhostSol_inviscidWallBC(nx,ny,&Ql_[0],&Qr_[0]);

        }else if(grid_->facelist[fID].bnd_type==-2){ // FarField BC

            compute_GhostSol_farfieldBC(nx,ny,&Ql_[0],&Qr_[0]);
        }
    }

    return;
}

void Euler2DSolver::Reconstruct_sol(){

    register int i; int j;

    for(i=0; i<Nelem; i++){

        for(j=0; j<Ndof; j++){
            dQdx[i][j] = 0.0;
            dQdy[i][j] = 0.0;
        }
    }

    return;
}

void Euler2DSolver::Residulas_setZero(double** Residuals){
    register int i,j;
    for (i=0; i<Nelem; i++)
        for(j=0; j<Ndof; j++)
            Residuals[i][j]=0.0;
    return;
}

void Euler2DSolver::compute_resid_OneFace(const int &fID, double* resid_){

    int j=0;

    for(j=0; j<Ndof; j++)
        resid_[j] = flux_com[fID][j] * grid_->facelist[fID].Af;

    return ;
}

void Euler2DSolver::Compute_common_inviscidflux(const double *Ql, const double *Qr
                                                   , const double& nx, const double& ny
                                                   , double* flux_){

    if(simdata_->ReimannSolver_type_=="Rusanov"){

        Rusanov_flux(&Ql[0], &Qr[0], nx, ny, &flux_[0]);

    }else if(simdata_->ReimannSolver_type_=="Roe"){

        _notImplemented("( ROE ) Reimann Solver is not implemented yet");

    }else{
        _notImplemented("This Reimann Solver type is not implemented yet");
    }

    return;
}

void Euler2DSolver::Rusanov_flux(const double *Ql, const double *Qr
             , const double& nx, const double& ny, double* flux_){

    int j=0;
    double FL[4]={0.,0.,0.,0.}, FR[4]={0.,0.,0.,0.};
    double Vn_abs,Vnl,Vnr,cl,cr,c,dQ=0.0;

//    Vnl = (Ql[1]/Ql[0])*nx + (Ql[2]/Ql[0])*ny;
//    Vnr = (Qr[1]/Qr[0])*nx + (Qr[2]/Qr[0])*ny;

    compute_normal_inViscidFlux(nx,ny,&Ql[0],&FL[0], Vnl,cl);
    compute_normal_inViscidFlux(nx,ny,&Qr[0],&FR[0], Vnr,cr);

    Vn_abs = abs(Vnl+Vnr)/2.;

    c = (cl+cr)/2.;

    for(j=0; j<Ndof; j++){
        dQ = Qr[j] - Ql[j];
        flux_[j] = 0.5 * (FL[j] + FR[j] - (Vn_abs + c) * dQ ) ;
    }

    return;
}

void Euler2DSolver::compute_normal_inViscidFlux(const double& nx, const double& ny
                                                , const double *Q_, double *normInvflux_
                                                , double& Vn, double& c){

    double rho,u,v,p,E;
    double gama_=gasdata_->gama;

    rho = Q_[0];
    u = Q_[1]/Q_[0];
    v = Q_[2]/Q_[0];
    E = Q_[3];

    p = (gama_-1.) * ( E - 0.5 * rho * ( pow(u,2) + pow(v,2) ) );

    c = sqrt(gama_*p/rho);

    Vn = u*nx + v*ny;

    normInvflux_[0] = rho * Vn;
    normInvflux_[1] = rho * u * Vn + p*nx;
    normInvflux_[1] = rho * v * Vn + p*ny;
    normInvflux_[3] = Vn*(E+p);

    return;
}


















