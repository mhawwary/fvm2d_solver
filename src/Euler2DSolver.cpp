
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
    Ndof = 4;

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
        Qv[i]    = new double[Ndof+1];

    flux_com = new double*[Nfaces];

    for(i=0; i<Nfaces; i++)
        flux_com[i]    = new double[Ndof];

    SetPhyTime(simdata_->t_init_);  // initial physical time

    dt_ = simdata_->dt_;

    CFL_ = simdata_->CFL_;

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

    register int i,j;

    SetGhostVariables();

    if(scheme_order==2) Reconstruct_sol();

    // Boundary Face loop to calculate the common interface fluxes:
    //----------------------------------------------------
    double Ql[4]={0.,0.,0.,0.}, Qr[4]={0.,0.,0.,0.};
    double nx=0.0,ny=0.0;

    // Boundary faces:

    for(j=0; j<Nbndfaces; j++){

        nx = grid_->facelist[j].nx;
        ny = grid_->facelist[j].ny;

        Compute_left_right_boundfacesol(j,&Ql[0],&Qr[0]);

        Compute_common_inviscidflux(&Ql[0],&Qr[0],nx,ny, &flux_com[j][0]);
    }

    // Interior faces
    for(j=Nbndfaces; j<Nfaces; j++){

        nx = grid_->facelist[j].nx;
        ny = grid_->facelist[j].ny;

        Compute_left_right_facesol(j,&Ql[0],&Qr[0]);

        Compute_common_inviscidflux(&Ql[0],&Qr[0],nx,ny, &flux_com[j][0]);
    }

    /* Face loop to update the Residuals in each respective cell
     * Two cells are updated at a time
     */

    // Set Residuals to zero before computing:
    Residulas_setZero(Resid_);

    int iL=0,iR=0;

    double resid_temp[4]={0.,0.,0.,0.};

    // Note that there should be a negative sign added to all Residuals;

/*    // Bound Faces
//    for(i=0; i<Nbndfaces; i++){

//        iL = grid_->facelist[i].Lcell;
//        Vol_l = grid_->elemlist[iL].Vc;

//        compute_resid_OneFace(grid_->facelist[i].Af, &flux_com[i][0], &resid_temp[0]);

//        for(j=0; j<Ndof; j++)
//            Resid_[iL][j] += (resid_temp[j]/Vol_l);
//    }

//    // Interior Faces:

//    for(i=Nbndfaces; i<Nfaces; i++){

//        iL = grid_->facelist[i].Lcell;
//        iR = grid_->facelist[i].Rcell;

//        Vol_l = grid_->elemlist[iL].Vc;
//        Vol_r = grid_->elemlist[iR].Vc;

//        compute_resid_OneFace(grid_->facelist[i].Af, &flux_com[i][0], &resid_temp[0]);

//        for(j=0; j<Ndof; j++){
//            Resid_[iL][j] +=  ( resid_temp[j]/Vol_l);
//            Resid_[iR][j] +=  (-resid_temp[j]/Vol_r);
//        }
//    } */

    // Bound Faces
    for(i=0; i<Nbndfaces; i++){

        iL = grid_->facelist[i].Lcell;

        compute_resid_OneFace(grid_->facelist[i].Af, &flux_com[i][0], &resid_temp[0]);

        for(j=0; j<Ndof; j++)
            Resid_[iL][j] += (resid_temp[j]);
    }

    // Interior Faces:

    for(i=Nbndfaces; i<Nfaces; i++){

        iL = grid_->facelist[i].Lcell;
        iR = grid_->facelist[i].Rcell;

        compute_resid_OneFace(grid_->facelist[i].Af, &flux_com[i][0], &resid_temp[0]);

        for(j=0; j<Ndof; j++){
            Resid_[iL][j] +=  ( resid_temp[j]);
            Resid_[iR][j] +=  (-resid_temp[j]);
        }
    }

    return;
}

void Euler2DSolver::SetGhostVariables(){

    register int i; int j=0;

    int iL,iR;

    double nx=0.0,ny=0.0;

    for(i=0; i<Nbndfaces; i++){

        iL = grid_->facelist[i].Lcell;
        iR = grid_->facelist[i].Rcell;
        nx = grid_->facelist[i].nx;
        ny = grid_->facelist[i].ny;

        if(simdata_->WallBC=="Extrapolation"){

            for(j=0; j<Ndof; j++) Qc[iR][j] = Qc[iL][j];

            /*if(iR!=i+Nelem) { _compare(iR,i+Nelem); std::cin.get(); }

            if(grid_->facelist[i].bnd_type==0)
                FatalErrorST("face bndtype is wrong");*/

        }else{

            if(grid_->facelist[i].bnd_type==-1){ // wall B.C

                compute_GhostSol_WallBC(nx,ny,&Qc[iL][0],&Qc[iR][0]);

            }else if(grid_->facelist[i].bnd_type==-2){ // FarField B.C

                compute_GhostSol_farfieldBC(nx,ny,&Qc[iL][0],&Qc[iR][0]);

            }else{
                FatalError_exit("This ghost cell boundary type is wrong");
            }
        }
    }

    return;
}

void Euler2DSolver::compute_GhostSol_WallBC(const double& nx, const double& ny
                                                    ,double *Ql, double *Qr){
    // InViscid Slip Boundary Condition:
    //-----------------------------------------
    double ul=0.0,vl=0.0,Vnl,ur=0.0,vr=0.0;

    ul = Ql[1]/Ql[0];
    vl = Ql[2]/Ql[0];

    Vnl = (ul*nx + vl*ny); // normal component of left cell velocity vector

    ur = ul - (2.*Vnl * nx);
    vr = vl - (2.*Vnl * ny);

    Qr[0] = Ql[0];
    Qr[1] = Qr[0] * ur;
    Qr[2] = Qr[0] * vr;
    Qr[3] = Ql[3];   // since they are of same velocity magnitude and Pr=Pl

    return;
}

void Euler2DSolver::compute_GhostSol_SymmetryBC(const double& nx, const double& ny
                                                    ,double *Ql, double *Qr){
    // InViscid Slip (Symmetry) Boundary Condition:
    //----------------------------------------------
    double ul=0.0,vl=0.0,Vnl,ur=0.0,vr=0.0;

    ul = Ql[1]/Ql[0];
    vl = Ql[2]/Ql[0];

    Vnl = (ul*nx + vl*ny); // normal component of left cell velocity vector

    ur = ul - (2.*Vnl * nx);
    vr = vl - (2.*Vnl * ny);

    Qr[0] = Ql[0];
    Qr[1] = Qr[0] * ur;
    Qr[2] = Qr[0] * vr;
    Qr[3] = Ql[3];   // since they are of same velocity magnitude and Pr=Pl

    return;
}

void Euler2DSolver::compute_GhostSol_farfieldBC(const double& nx, const double& ny
                                                ,double *Ql, double *Qr){

    if(simdata_->FarFieldBC=="Fixed"){

        double rho,u,v,E;
        gasdata_->CalculateFlowProperties(rho,u,v,E);

        Qr[0] = rho;
        Qr[1] = rho*u;
        Qr[2] = rho*v;
        Qr[3] = E;

    }else if(simdata_->FarFieldBC=="Characteristics"){

        // First test If Inlet or Exit :
        // -----------------------------
        double u,v,Vn;
        u = Ql[1]/Ql[0]; v = Ql[2]/Ql[0];
        Vn = u*nx + v*ny;

        if(Vn < 0){ // Inlet Boundary condition:
            Compute_Inlet_charBC(nx,ny,Ql,Qr);
        } else if(Vn>0) { // Exit Boundary Condition
            Compute_Exit_charBC(nx,ny,Ql,Qr);
        } else{ // Vn=0, Symmetry boundary condition
            //FatalError("Very Strange Error Vn=0 at Farfeild BC");
            compute_GhostSol_SymmetryBC(nx,ny,Ql,Qr);
        }
    }else{
        _notImplemented("This FarField Boundary Condition is not implemented");
    }

    return;
}

void Euler2DSolver::Compute_Inlet_charBC(const double &nx, const double &ny
                                         , double *Ql, double *Qr){
    double Vn_l, u, v , c_l, Pl;

    u = Ql[1] / Ql[0];
    v = Ql[2] / Ql[0];

    Pl = (gasdata_->gama-1.0) * ( Ql[3] - ( 0.5 * ( pow(Ql[1],2) + pow(Ql[2],2) ) / Ql[0] ) );
    Vn_l = u*nx+v*ny;
    c_l  = sqrt ( gasdata_->gama * Pl / Ql[0] );

    double Mn = fabs(Vn_l) / c_l;

    if(Mn>=1) { // Supersonic Inflow
        //printf("\n Supersonic Inflow ");
        double rho_r, u_r, v_r, E_r;
        gasdata_->CalculateFlowProperties(rho_r,u_r,v_r,E_r);
        Qr[0] = rho_r;
        Qr[1] = rho_r * u_r;
        Qr[1] = rho_r * v_r;
        Qr[3] = E_r;

    }else {     // Subsonic Inflow
        double gama = gasdata_->gama;
        double Po = gasdata_->Ptot;
        double To = gasdata_->T_tot;

        double alpha = gasdata_->alpha_deg_ * PI / 180.0;

        double Pr,Tr,rho_r,Mr,c_r;

        Pr = Pl; // Pr = Pl

        double gm1_g = (gama-1.0)/gama;

        Mr = sqrt( 2.0 * ( pow((Po/Pr),gm1_g) - 1.0 ) / (gama-1.0) ) ;

        Tr = To / ( 1.0 + (0.5*(gama-1.0) * pow(Mr,2)) ) ;

        rho_r = Pr / ( gasdata_->R_gas * Tr );

        c_r = sqrt(gama*Pr/rho_r);

        Qr[0] = rho_r;
        Qr[1] = Qr[0] * Mr * c_r * cos(alpha);
        Qr[2] = Qr[0] * Mr * c_r * sin(alpha);
        Qr[3] = (Pr/(gama-1.0)) + ( 0.5 * ( pow(Qr[1],2) + pow(Qr[2],2) ) / Qr[0] ) ;
    }

    return;
}

void Euler2DSolver::Compute_Exit_charBC(const double &nx, const double &ny
                                        , double *Ql, double *Qr){
    double Vn_l, u, v , c_l, Pl;

    u = Ql[1] / Ql[0];
    v = Ql[2] / Ql[0];

    Pl = (gasdata_->gama-1.) * ( Ql[3] - ( 0.5 * ( pow(Ql[1],2) + pow(Ql[2],2) ) / Ql[0] ) );
    Vn_l = u*nx+v*ny;
    c_l  = sqrt ( gasdata_->gama * Pl / Ql[0] );

    double Mn = fabs(Vn_l) / c_l;

    if(Mn>=1) { // Supersonic Outflow
        //printf("\n Supersonic Outflow ");
        Qr[0] = Ql[0];
        Qr[1] = Ql[1];
        Qr[2] = Ql[2];
        Qr[3] = Ql[3];

    } else { // subsonic Outflow

        double gama = gasdata_->gama;
        double u_r,v_r,Pr,Vn_r,c_r,rho_r,Ro,Rp;
        Pr = gasdata_->Ps;

        Ro = Pl / pow( Ql[0],gama ) ;
        Rp = Vn_l + ( 2 * c_l / (gama-1.0) ) ;

        rho_r = pow( (Pr/Ro) , (1.0/gama) );
        c_r = sqrt( gama * Pr / rho_r );
        Vn_r = Rp - ( 2 * c_r / (gama-1.0) );
        u_r = u  + ( Vn_r - Vn_l) * nx;
        v_r = v  + ( Vn_r - Vn_l) * ny;

        Qr[0] = rho_r;
        Qr[1] = rho_r * u_r;
        Qr[2] = rho_r * v_r;
        Qr[3] = (Pr/(gama-1.0)) + 0.5 * rho_r * (pow(u_r,2)+pow(v_r,2));
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
    ny = grid_->facelist[fID].ny;

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

        for(j=0; j<Ndof; j++) Qr_[j] = Ql_[j];

    }else{

        if(grid_->facelist[fID].bnd_type==-1){ // wall BC

            compute_GhostSol_WallBC(nx,ny,&Ql_[0],&Qr_[0]);

        }else if(grid_->facelist[fID].bnd_type==-2){ // FarField BC

            compute_GhostSol_farfieldBC(nx,ny,&Ql_[0],&Qr_[0]);
        }
    }

    return;
}

void Euler2DSolver::Reconstruct_sol(){

    register int i; int j,k,fID,iL,iR,ii;
    double dUx=0.0,dUy=0.0,dU=0.0;

    for(i=0; i<Nelem; i++){

        for(j=0; j<Ndof; j++){

            dUx = 0.0; dUy=0.0;

            for(k=0; k<grid_->elemlist[i].n_local_faces; k++){

                fID = grid_->elemlist[i].to_face[k];
                iL = grid_->facelist[fID].Lcell;
                iR = grid_->facelist[fID].Rcell;

                if(iL!=i) ii=iL;
                else ii=iR;

                dU = ( Qc[ii][j] - Qc[i][j] );
                dUx += dU * grid_->elemlist[i].DX[k];
                dUy += dU * grid_->elemlist[i].DY[k];
            }

            dQdx[i][j] = grid_->Iyy[i] * dUx - grid_->Ixy[i] * dUy;
            dQdy[i][j] = grid_->Ixy[i] * dUx - grid_->Ixx[i] * dUy;
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

void Euler2DSolver::compute_resid_OneFace(const double &face_area
                                          ,double* face_flux_, double* resid_){

    int j=0;

    for(j=0; j<Ndof; j++)
        resid_[j] = (-1.0)* face_flux_[j] * face_area ;
    // the minus sign is required to move it to the RHS of the equation

    return ;
}

void Euler2DSolver::Compute_common_inviscidflux(double *Ql, double *Qr
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

void Euler2DSolver::Rusanov_flux(double *Ql, double *Qr
             , const double& nx, const double& ny, double* flux_){

    int j=0;
    double FL[4]={0.,0.,0.,0.}, FR[4]={0.,0.,0.,0.};
    double Vn_abs,Vnl,Vnr,pl,pr,c,dQ=0.0;

    compute_normal_inViscidFlux(nx,ny,&Ql[0],&FL[0], Vnl,pl);
    compute_normal_inViscidFlux(nx,ny,&Qr[0],&FR[0], Vnr,pr);

    Vn_abs = fabs(Vnl+Vnr)/2.;

    c = sqrt(gasdata_->gama *(pl+pr)/(Ql[0]+Qr[0]));

    for(j=0; j<Ndof; j++){
        dQ = Qr[j] - Ql[j];
        flux_[j] = 0.5 * ( FL[j] + FR[j] - (Vn_abs + c) * dQ ) ;
    }

    return;
}

void Euler2DSolver::compute_normal_inViscidFlux(const double& nx, const double& ny
                                                , double *Q_, double *normInvflux_
                                                , double& Vn, double& p){

    double rho,u,v,E;
    double gama_=gasdata_->gama;

    rho = Q_[0];
    u = Q_[1]/Q_[0];
    v = Q_[2]/Q_[0];
    E = Q_[3];

    p = (gama_-1.) * ( E - (0.5 * rho * ( pow(u,2) + pow(v,2) )) );

    //c = sqrt(gama_*p/rho);

    Vn = u*nx + v*ny;

    normInvflux_[0] = rho * Vn;
    normInvflux_[1] = (rho * u * Vn) + (p*nx);
    normInvflux_[2] = (rho * v * Vn) + (p*ny);
    normInvflux_[3] = Vn*(E+p);

    return;
}

void Euler2DSolver::evaluate_sol(double &Xp, double &Yp, const int& eID, double *qq_){

    int j;

    double Xc=grid_->elemlist[eID].Xc;
    double Yc=grid_->elemlist[eID].Yc;

    if(scheme_order==1){
        for(j=0; j<Ndof; j++)  qq_[j] = Qc[eID][j];

    }else if(scheme_order==2){

        for(j=0; j<Ndof; j++){
            qq_[j] = Qc[eID][j] + dQdx[eID][j] * (Xp-Xc) + dQdy[eID][j] * (Yp-Yc);
        }
    }

    return;
}

void Euler2DSolver::Compute_vertex_sol(){

    register int i; int j,eID,iL,fID;

    double rho=0.0,u=0.0,v=0.0,p=0.0,p_inf=0.0,rhoV_inf=0.,V_inf=0.,rho_inf,T,E;
    double qq_[4] ={0.,0.,0.,0.};
    double nx,ny;

    for(i=0; i<Nnodes; i++){

        Qv[i][0] =0.;
        Qv[i][1] =0.;
        Qv[i][2] =0.;
        Qv[i][3] =0.;
        Qv[i][4] =0.;

        for(j=0; j<grid_->Nnode_neighElem[i]; j++){

            eID = grid_->node_to_elemlist[i][j];

            gasdata_->CalculateFlowProperties(rho_inf,u,v,p_inf,T,E);

            V_inf = sqrt(pow(u,2)+pow(v,2));

            rhoV_inf= rho_inf * V_inf * V_inf ;

            if(eID>=Nelem){

                fID = eID-Nelem;
                iL = grid_->facelist[fID].Lcell;
                nx = grid_->facelist[fID].nx;
                ny = grid_->facelist[fID].ny;

                evaluate_sol(grid_->Xn[i],grid_->Yn[i],iL, &qq_[0]);

                if(grid_->facelist[fID].bnd_type==-1)
                    compute_GhostSol_WallBC(nx,ny,&qq_[0],&qq_[0]); // becareful of debendency between ql,qr
                else if(grid_->facelist[fID].bnd_type==-2)
                    compute_GhostSol_farfieldBC(nx,ny,&qq_[0],&qq_[0]);
                else
                    FatalError("Wrong face boundary condition for vertex solution");
            }else{
                evaluate_sol(grid_->Xn[i],grid_->Yn[i],eID, &qq_[0]);
            }

            rho = qq_[0];
            u = qq_[1]/qq_[0];
            v = qq_[2]/qq_[0];

            p = (gasdata_->gama-1.) * ( qq_[3] - 0.5 * rho * ( pow(u,2) + pow(v,2) ) );

            Qv[i][0] += rho/rho_inf;
            Qv[i][1] += u/V_inf;
            Qv[i][2] += v/V_inf;
            Qv[i][3] += ( 2.*(p-p_inf)/ rhoV_inf );
            Qv[i][4] += sqrt((u*u)+(v*v)) / sqrt(gasdata_->gama * p /rho);
        }

        for(j=0; j<Ndof; j++)
            Qv[i][j] = Qv[i][j] / grid_->Nnode_neighElem[i];
    }

    return;
}

void Euler2DSolver::Compute_local_TimeStep(double *dt_cell_){

    double *cell_radii=nullptr; // sum(|Vn|+c * Sf)

    cell_radii = new double[Nelem];

    register int i; int iL,iR;

    double nx,ny,Sf,Vn;
    double rho,u,v,p,E,c;
    double gama_=gasdata_->gama;

    double Ql_[4]={0.,0.,0.,0.},Qr_[4]={0.,0.,0.,0.};

    for(i=0; i<Nelem; i++) cell_radii[i]=0.0;

    for(i=0; i<Nbndfaces; i++){

        iL = grid_->facelist[i].Lcell;
        nx = grid_->facelist[i].nx;
        ny = grid_->facelist[i].ny;
        Sf = grid_->facelist[i].Af;

        Compute_left_right_boundfacesol(i,Ql_,Qr_);

        rho = Ql_[0];
        u = Ql_[1]/Ql_[0];
        v = Ql_[2]/Ql_[0];
        E = Ql_[3];

        p = (gama_-1.) * ( E - 0.5 * rho * ( pow(u,2) + pow(v,2) ) );

        c = sqrt(gama_*p/rho);

        Vn = u*nx + v*ny;

        rho = Qr_[0];
        u = Qr_[1]/Qr_[0];
        v = Qr_[2]/Qr_[0];
        E = Qr_[3];

        p = (gama_-1.) * ( E - 0.5 * rho * ( pow(u,2) + pow(v,2) ) );

        c += sqrt(gama_*p/rho);

        Vn += (u*nx + v*ny);

        Vn = fabs(Vn/2.);
        c = c/2.;

        cell_radii[iL] += ((Vn+c)*Sf);
    }

    for(i=Nbndfaces; i<Nfaces; i++){

        iL = grid_->facelist[i].Lcell;
        iR = grid_->facelist[i].Rcell;
        nx = grid_->facelist[i].nx;
        ny = grid_->facelist[i].ny;
        Sf = grid_->facelist[i].Af;

        Compute_left_right_facesol(i,Ql_,Qr_);

        rho = Ql_[0];
        u = Ql_[1]/Ql_[0];
        v = Ql_[2]/Ql_[0];
        E = Ql_[3];

        p = (gama_-1.) * ( E - 0.5 * rho * ( pow(u,2) + pow(v,2) ) );

        c = sqrt(gama_*p/rho);

        Vn = u*nx + v*ny;

        rho = Qr_[0];
        u = Qr_[1]/Qr_[0];
        v = Qr_[2]/Qr_[0];
        E = Qr_[3];

        p = (gama_-1.) * ( E - 0.5 * rho * ( pow(u,2) + pow(v,2) ) );

        c += sqrt(gama_*p/rho);

        Vn += (u*nx + v*ny);

        Vn = fabs(Vn/2.);
        c = c/2.;

        cell_radii[iL] += ((Vn+c)*Sf);
        cell_radii[iR] += ((Vn+c)*Sf);
    }

    for(i=0; i<Nelem; i++) {

        cell_radii[i] = cell_radii[i] / (2.0 * grid_->elemlist[i].Vc);

        dt_cell_[i] = CFL_ / cell_radii[i];
    }

    emptyarray(cell_radii);

    return;
}
















