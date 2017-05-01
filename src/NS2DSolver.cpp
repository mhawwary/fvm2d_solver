
#include"NS2DSolver.hpp"

// Constructor/Destructor/ Setup functions:
//------------------------------------------------

NS2DSolver::~NS2DSolver(void){
    Reset_solver();
}

void NS2DSolver::setup_solver(MeshData*& meshdata_, SimData& osimdata_
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

    Qv = new double*[grid_->Nnodes_postproc];

    for(i=0; i<grid_->Nnodes_postproc; i++)
        Qv[i]    = new double[Ndof+1];

    flux_com = new double*[Nfaces];
    visc_flux_com = new double*[Nfaces];

    for(i=0; i<Nfaces; i++){
        flux_com[i] = new double[Ndof];
        visc_flux_com[i] = new double[Ndof];
    }

    SetPhyTime(simdata_->t_init_);  // initial physical time

    dt_ = simdata_->dt_;

    CFL_ = simdata_->CFL_;

    return;
}

void NS2DSolver::Reset_solver(){

    //FatalError("Inside NS2DSolver Destructor");

    emptyarray(Nelem_extend,Qc);
    emptyarray(Nelem,dQdx);
    emptyarray(Nelem,dQdy);

    //FatalError("After empty Qn");

    emptyarray(grid_->Nnodes_postproc,Qv);

    //FatalError("After empty Qv");

    emptyarray(Nfaces,flux_com);
    emptyarray(Nfaces,visc_flux_com);

    return;
}

void NS2DSolver::InitSol(){

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

void NS2DSolver::UpdateResid(double **Resid_, double **Qc_){

    register int i,j; int k,iL,iR;

    SetGhostVariables();

    if(scheme_order==2) Reconstruct_sol();

    // Boundary Face loop to calculate the common interface fluxes:
    //----------------------------------------------------
    double Ql[4]={0.,0.,0.,0.}, Qr[4]={0.,0.,0.,0.};
    double nx=0.0,ny=0.0;

    // Boundary faces:

    double **dQl_=nullptr,**dQr_=nullptr;

    dQl_ = new double*[2]; dQr_ = new double*[2];
    dQl_[0] = new double[Ndof]; dQr_[0] = new double[Ndof];
    dQl_[1] = new double[Ndof]; dQr_[1] = new double[Ndof];

    double Qf_[4] = {0.,0.,0.,0.};

    for(j=0; j<Nbndfaces; j++){

        nx = grid_->facelist[j].nx;
        ny = grid_->facelist[j].ny;
        iL= grid_->facelist[j].Lcell;
        iR= grid_->facelist[j].Rcell;

        Compute_left_right_boundfacesol(j,&Ql[0],&Qr[0]);

        Compute_common_inviscidflux(&Ql[0],&Qr[0],nx,ny, &flux_com[j][0]);

        for(k=0; k<Ndof; k++){
            dQl_[0][k] = dQdx[iL][k]; dQl_[1][k] = dQdy[iL][k];
            dQr_[0][k] = 0.0; dQr_[1][k] = 0.0;
        }

        Compute_common_facesol(&Ql[0],&Qr[0],&Qf_[0]);

        Compute_common_bound_viscousflux(j,&Qc[iL][0],&Qc[iR][0],&Qf_[0]
                , dQl_,dQr_,&visc_flux_com[j][0]);
    }

    // Interior faces
    for(j=Nbndfaces; j<Nfaces; j++){

        nx = grid_->facelist[j].nx;
        ny = grid_->facelist[j].ny;
        iL= grid_->facelist[j].Lcell;
        iR= grid_->facelist[j].Rcell;

        Compute_left_right_facesol(j,&Ql[0],&Qr[0]);

        Compute_common_inviscidflux(&Ql[0],&Qr[0],nx,ny, &flux_com[j][0]);

        for(k=0; k<Ndof; k++){
            dQl_[0][k] = dQdx[iL][k]; dQl_[1][k] = dQdy[iL][k];
            dQr_[0][k] = dQdx[iR][k]; dQr_[1][k] = dQdy[iR][k];
        }

        Compute_common_facesol(&Ql[0],&Qr[0],&Qf_[0]);

        Compute_common_viscousflux(j,&Qc[iL][0],&Qc[iR][0],&Qf_[0]
                , dQl_,dQr_,&visc_flux_com[j][0]);
    }

    emptyarray(2,dQl_);
    emptyarray(2,dQr_);

    /* Face loop to update the Residuals in each respective cell
     * Two cells are updated at a time
     */

    // Set Residuals to zero before computing:
    Residulas_setZero(Resid_);

    double resid_temp[4]={0.,0.,0.,0.};

    // Bound Faces
    for(i=0; i<Nbndfaces; i++){

        iL = grid_->facelist[i].Lcell;

        //printf("\nboundFace: %d, type: %d\n",i, grid_->facelist[i].bnd_type);

        compute_resid_OneFace(grid_->facelist[i].Af, &flux_com[i][0]
                , &visc_flux_com[i][0], &resid_temp[0]);

        for(j=0; j<Ndof; j++)
            Resid_[iL][j] += (resid_temp[j]);
    }

    // Interior Faces:

    for(i=Nbndfaces; i<Nfaces; i++){

        iL = grid_->facelist[i].Lcell;
        iR = grid_->facelist[i].Rcell;

        //printf("\nInteriorFace: %d, type: %d\n",i, grid_->facelist[i].bnd_type);

        compute_resid_OneFace(grid_->facelist[i].Af, &flux_com[i][0]
                , &visc_flux_com[i][0], &resid_temp[0]);

        for(j=0; j<Ndof; j++){
            Resid_[iL][j] +=  ( resid_temp[j]);
            Resid_[iR][j] +=  (-resid_temp[j]);
        }
    }

    return;
}

void NS2DSolver::SetGhostVariables(){

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

void NS2DSolver::compute_GhostSol_WallBC(const double& nx, const double& ny
                                                    ,double *Ql, double *Qr){
    // Viscous No Slip Boundary Condition:
    //-----------------------------------------
    double ul=0.0,vl=0.0,ur=0.0,vr=0.0;

    ul = Ql[1]/Ql[0];
    vl = Ql[2]/Ql[0];

    ur = -ul;
    vr = -vl;

    Qr[0] = Ql[0];
    Qr[1] = Qr[0] * ur;
    Qr[2] = Qr[0] * vr;
    Qr[3] = Ql[3];   // since they are of same velocity magnitude and Pr=Pl

    return;
}

void NS2DSolver::compute_GhostSol_SymmetryBC(const double& nx, const double& ny
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

void NS2DSolver::compute_GhostSol_farfieldBC(const double& nx, const double& ny
                                                ,double *Ql, double *Qr){

    if(simdata_->FarFieldBC=="Extrapolation"){
        int j=0;
        for(j=0; j<Ndof; j++) Qr[j] = Ql[j];

    }else if(simdata_->FarFieldBC=="Fixed"){

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

void NS2DSolver::Compute_Inlet_charBC(const double &nx, const double &ny
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

void NS2DSolver::Compute_Exit_charBC(const double &nx, const double &ny
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

void NS2DSolver::Compute_left_right_facesol(const int& fID
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

void NS2DSolver::Compute_left_right_boundfacesol(const int& fID
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

void NS2DSolver::Compute_common_facesol(const double *Ql_, const double *Qr_, double *Qf_){

    int k=0;

    // Simple Averaging:
    //-----------------------
    for(k=0; k<Ndof; k++)
        Qf_[k] = 0.5 * (Ql_[k]+Qr_[k]);

    return;
}

void NS2DSolver::Compute_common_face_gradient(const int &fID, const double *Ql, const double *Qr
                                              ,double **dQl_,double **dQr_, double **dQf_){

    // Unbiased Averaging:
    //---------------------------

    int k,iL,iR;
    double Xcl,Xcr,Ycl,Ycr,lx,ly,dl;

    iL = grid_->facelist[fID].Lcell;
    iR = grid_->facelist[fID].Rcell;

    Xcl = grid_->elemlist[iL].Xc;
    Ycl = grid_->elemlist[iL].Yc;
    Xcr = grid_->elemlist[iR].Xc;
    Ycr = grid_->elemlist[iR].Yc;

    lx = (Xcr-Xcl);
    ly = (Ycr-Ycl);
    dl = sqrt( pow(lx,2) + pow(ly,2) );
    lx = lx/dl; ly = ly/dl;

    if(simdata_->ViscFlux=="ZJWang"){
        int v0,v1;
        double tx,ty,dtau;
        v0 = grid_->facelist[fID].v0 ;
        v1 = grid_->facelist[fID].v1 ;

        tx = grid_->Xn[v1] - grid_->Xn[v0] ;
        ty = grid_->Yn[v1] - grid_->Yn[v0] ;
        dtau = sqrt( pow(tx,2) + pow(ty,2) );
        tx = tx/dtau; ty = ty/dtau;

        double dQdl,dQl_dtau,dQr_dtau,dQ_dtau;
        for(k=0; k<Ndof; k++){

            dQdl = (Qr[k] - Ql[k] )/ dl;

            if(grid_->facelist[fID].bnd_type==-1){

                dQ_dtau = 0.0;

            }else if(grid_->facelist[fID].bnd_type==0){

                dQl_dtau = dQl_[0][k] * tx + dQl_[1][k] * ty;
                dQr_dtau = dQr_[0][k] * tx + dQr_[1][k] * ty;
                dQ_dtau = 0.5 * (dQl_dtau+dQr_dtau);
            }

            double DD = lx*ty - ly *tx; if(DD==0) FatalError_exit("\nProblem in Face Gradient \n");
            dQf_[0][k] = ( (dQdl *ty) - (dQ_dtau * ly) ) / DD;
            dQf_[1][k] = ( (-dQdl *tx) + (dQ_dtau * lx) ) / DD;

            if(std::isnan(dQf_[0][k]) || std::isnan(dQf_[1][k])){
                printf("\nTruobled face: %d, bndtype:%d gradients are nan: %e %e"
                       ,fID, grid_->facelist[fID].bnd_type, dQf_[0][k], dQf_[1][k]);
                std::cin.get();
            }
        }

    }else if(simdata_->ViscFlux=="CarlGooch"){

        double ex,ey,dQx,dQy,dQdl;

        ex = -ly; ey = lx;

        for(k=0; k<Ndof; k++){
            dQx = 0.5 *(dQl_[0][k] + dQr_[0][k]);
            dQy = 0.5 *(dQl_[1][k] + dQr_[1][k]);

            dQdl = (Qr[k] - Ql[k] )/ dl;

            dQf_[0][k] = dQdl * lx + (dQx*ex+dQy*ey)*ex;
            dQf_[1][k] = dQdl * ly + (dQx*ex+dQy*ey)*ey;
        }
    }else{
        _notImplemented("Viscous Flux treatment is not implemented");
    }

    return;
}

void NS2DSolver::Compute_common_viscousflux(const int &fID, const double *Ql
                                            , const double *Qr, const double *Qf_
                                            , double **dQl_,double **dQr_, double *flux_){
    double nx = grid_->facelist[fID].nx;
    double ny = grid_->facelist[fID].ny;

    double Ux,Uy,Vx,Vy,Tx,Ty,q_x,q_y;
    double Tau_xx,Tau_xy,Tau_yy,Tau;
    //double Qf_[4]={0.,0.,0.,0.};
    double **dQf_=nullptr;
    dQf_ = new double*[2];
    dQf_[0] = new double[Ndof];
    dQf_[1] = new double[Ndof];

//    Compute_common_facesol(Ql,Qr,Qf_);

    Compute_common_face_gradient(fID,Ql,Qr,dQl_,dQr_,dQf_);

    Ux = ( dQf_[0][1] - ( Qf_[1] * dQf_[0][0]/ pow(Qf_[0],2) ) ) / Qf_[0] ;
    Uy = ( dQf_[1][1] - ( Qf_[1] * dQf_[1][0]/ pow(Qf_[0],2) ) ) / Qf_[0] ;
    Vx = ( dQf_[0][2] - ( Qf_[2] * dQf_[0][0]/ pow(Qf_[0],2) ) ) / Qf_[0] ;
    Vy = ( dQf_[1][2] - ( Qf_[2] * dQf_[1][0]/ pow(Qf_[0],2) ) ) / Qf_[0] ;

    double gama = gasdata_->gama, Mu_i = gasdata_->Mu_inf;
    double rho_ = Qf_[0],u=Qf_[1]/Qf_[0],v=Qf_[2]/Qf_[0],t1,t2,t3;

    t1 = rho_ * dQf_[0][3];
    t2 = ( Qf_[3] - (pow(Qf_[1],2)+pow(Qf_[2],2))/rho_ ) *dQf_[0][0];
    t3 = Qf_[1] * dQf_[0][1] + Qf_[2] * dQf_[0][2];
    Tx = gama*(gama-1.0) * ( t1 - t2 - t3 ) / pow(rho_,2) ;

    t1 = rho_ * dQf_[1][3];
    t2 = ( Qf_[3] - (pow(Qf_[1],2)+pow(Qf_[2],2))/rho_ ) *dQf_[1][0];
    t3 = Qf_[1] * dQf_[1][1] + Qf_[2] * dQf_[1][2];
    Ty = gama*(gama-1.0) * ( t1 - t2 - t3 ) / pow(rho_,2) ;

    q_x = Mu_i * Tx / ((gama-1.0) * gasdata_->Prndtl);
    q_y = Mu_i * Ty / ((gama-1.0) * gasdata_->Prndtl);

    Tau = Ux + Vy;
    Tau_xx = Mu_i * ( 2*Ux + Lambda_stokes * Tau);
    Tau_yy = Mu_i * ( 2*Vy + Lambda_stokes * Tau);
    Tau_xy = Mu_i * ( Uy + Vx );

    flux_[0] = 0.0;
    flux_[1] = Tau_xx * nx + Tau_xy * ny;
    flux_[2] = Tau_xy * nx + Tau_yy * ny;
    flux_[3] = ( u*Tau_xx + v*Tau_xy + q_x ) * nx + ( u*Tau_xy + v*Tau_yy + q_y ) * ny ;

    emptyarray(2,dQf_);

    return;
}

void NS2DSolver::Compute_common_bound_viscousflux(const int &fID, const double *Ql
                                                  , const double *Qr, const double *Qf_
                                                  , double **dQl_,double **dQr_, double *flux_){

    if(simdata_->FarFieldBC=="Extrapolation"){

        flux_[0] = 0.0 ;
        flux_[1] = 0.0 ;
        flux_[2] = 0.0 ;
        flux_[3] = 0.0 ;

    }else{
        if(grid_->facelist[fID].bnd_type==-2){  // FarField B.C
            flux_[0] = 0.0 ;
            flux_[1] = 0.0 ;
            flux_[2] = 0.0 ;
            flux_[3] = 0.0 ;

        }else if(grid_->facelist[fID].bnd_type==-1) { // Viscous Wall B.C.
            double nx = grid_->facelist[fID].nx;
            double ny = grid_->facelist[fID].ny;
            double Ux,Uy,Vx,Vy;
            double Tau_xx,Tau_xy,Tau_yy,Tau;

            double **dQf_=nullptr;
            dQf_ = new double*[2];
            dQf_[0] = new double[4];
            dQf_[1] = new double[4];

            Compute_common_face_gradient(fID,Ql,Qr,dQl_,dQr_,dQf_);

            Ux = ( dQf_[0][1] ) / Qf_[0] ; // d(RhoU)_dx
            Uy = ( dQf_[1][1] ) / Qf_[0] ; // d(RhoU)_dy
            Vx = ( dQf_[0][2] ) / Qf_[0] ; // d(RhoV)_dx
            Vy = ( dQf_[1][2] ) / Qf_[0] ; // d(RhoV)_dy

            double Mu_i = gasdata_->Mu_inf;

            Tau = Ux + Vy;
            Tau_xx = Mu_i * ( 2*Ux + Lambda_stokes * Tau);
            Tau_yy = Mu_i * ( 2*Vy + Lambda_stokes * Tau);
            Tau_xy = Mu_i * ( Uy + Vx );

            flux_[0] = 0.0;
            flux_[1] = Tau_xx * nx + Tau_xy * ny;
            flux_[2] = Tau_xy * nx + Tau_yy * ny;
            flux_[3] = 0.0 ;

            emptyarray(2,dQf_);
        }else{
            _notImplemented("This Boundary Condition is not implemented");
        }
    }

    return;
}

void NS2DSolver::Reconstruct_sol(){

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

void NS2DSolver::Residulas_setZero(double** Residuals){
    register int i,j;
    for (i=0; i<Nelem; i++)
        for(j=0; j<Ndof; j++)
            Residuals[i][j]=0.0;
    return;
}

void NS2DSolver::compute_resid_OneFace(const double &face_area, double* inv_flux_
                                       , double* visc_flux_, double* resid_){

    int j=0;

    for(j=0; j<Ndof; j++){
        resid_[j] = ( visc_flux_[j] - inv_flux_[j] ) * face_area ;

        //printf("%d %e %e %e\n",j,inv_flux_[j],visc_flux_[j],resid_[j]);
        //std::cin.get();

    }
        //resid_[j] = (- inv_flux_[j] ) * face_area ;
    // Residual at the RHS of the equation

    return ;
}

void NS2DSolver::Compute_common_inviscidflux(double *Ql, double *Qr
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

void NS2DSolver::Rusanov_flux(double *Ql, double *Qr
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

void NS2DSolver::compute_normal_inViscidFlux(const double& nx, const double& ny
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

void NS2DSolver::evaluate_sol(double &Xp, double &Yp, const int& eID, double *qq_){

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

void NS2DSolver::Compute_vertex_sol(const int& oiter){

    if(comp_vertex_iter<oiter){
        comp_vertex_iter = oiter;

        register int i; int j,eID,iL,fID;

        double rho=0.0,u=0.0,v=0.0,p=0.0,p_inf=0.0,rhoV_inf=0.,V_inf=0.,rho_inf,T,E;
        double qq_[4] ={0.,0.,0.,0.};
        double nx,ny;

        gasdata_->CalculateFlowProperties(rho_inf,u,v,p_inf,T,E);
        V_inf = sqrt(pow(u,2)+pow(v,2));
        rhoV_inf= rho_inf * V_inf * V_inf ;

        for(i=0; i<Nnodes; i++){

            Qv[i][0] =0.;
            Qv[i][1] =0.;
            Qv[i][2] =0.;
            Qv[i][3] =0.;
            Qv[i][4] =0.;

            for(j=0; j<grid_->Nnode_neighElem[i]; j++){

                eID = grid_->node_to_elemlist[i][j];

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

                Qv[i][0] += rho;
                Qv[i][1] += u;
                Qv[i][2] += v;
                Qv[i][3] += p;
            }

            for(j=0; j<Ndof; j++)
                Qv[i][j] = Qv[i][j] / grid_->Nnode_neighElem[i];

            Qv[i][4] = sqrt( pow(Qv[i][1],2) + pow(Qv[i][2],2) )
                    / sqrt(gasdata_->gama * Qv[i][3] /Qv[i][0]);
            Qv[i][0] = Qv[i][0] / rho_inf;
            Qv[i][1] = Qv[i][1] / V_inf;
            Qv[i][2] = Qv[i][2] / V_inf;
            Qv[i][3] = 2.0* (Qv[i][3] - p_inf)/rhoV_inf;
        }

        // Another loop for the extended node list
        // that contains new cells centers:
        //-----------------------------------------

        for(i=Nnodes; i<grid_->Nnodes_postproc; i++){

            eID = grid_->polygon_elem_origID[i-Nnodes];

            rho = Qc[eID][0];
            u = Qc[eID][1]/Qc[eID][0];
            v = Qc[eID][2]/Qc[eID][0];

            p = (gasdata_->gama-1.) * ( Qc[eID][3] - 0.5 * rho * ( pow(u,2) + pow(v,2) ) );

            Qv[i][0] = rho/rho_inf;
            Qv[i][1] = u/V_inf;
            Qv[i][2] = v/V_inf;
            Qv[i][3] = 2.0* (p - p_inf)/rhoV_inf;
            Qv[i][4] = sqrt(pow(u,2)+pow(v,2)) / sqrt(gasdata_->gama * p / rho);
        }
    }

    return;
}

void NS2DSolver::Compute_local_TimeStep(double *dt_cell_){

    double *cell_radii=nullptr; // sum(|Vn|+c * Sf)

    cell_radii = new double[Nelem];

    register int i; int iL,iR;

    double nx,ny,Sf,Vn;
    double rho,u,v,p,E,c;
    double gama_=gasdata_->gama;
    double rho_l,rho_r,Bv1,Bv2;

    double Ql_[4]={0.,0.,0.,0.},Qr_[4]={0.,0.,0.,0.};

    for(i=0; i<Nelem; i++) cell_radii[i]=0.0;

    for(i=0; i<Nbndfaces; i++){

        iL = grid_->facelist[i].Lcell;
        nx = grid_->facelist[i].nx;
        ny = grid_->facelist[i].ny;
        Sf = grid_->facelist[i].Af;

        Compute_left_right_boundfacesol(i,Ql_,Qr_);

        rho = Ql_[0];  rho_l=rho;
        u = Ql_[1]/Ql_[0];
        v = Ql_[2]/Ql_[0];
        E = Ql_[3];

        p = (gama_-1.) * ( E - 0.5 * rho * ( pow(u,2) + pow(v,2) ) );

        c = sqrt(gama_*p/rho);

        Vn = u*nx + v*ny;

        rho = Qr_[0];  rho_r=rho;
        u = Qr_[1]/Qr_[0];
        v = Qr_[2]/Qr_[0];
        E = Qr_[3];

        p = (gama_-1.) * ( E - 0.5 * rho * ( pow(u,2) + pow(v,2) ) );

        c += sqrt(gama_*p/rho);

        Vn += (u*nx + v*ny);

        Vn = fabs(Vn/2.);
        c = c/2.;

        // Viscous contribution terms:
        rho = 0.5 * (rho_l +rho_r);
        Bv1 = (gama_/rho), Bv2 = 4.0/ ( 3.0 * rho) ;
        if(Bv1<Bv2) Bv1 =Bv2;

        Bv1 = Bv1 * gasdata_->Mu_inf * Sf * Sf / gasdata_->Prndtl;

        cell_radii[iL] += ( ((Vn+c)*Sf*0.5) + (Bv1 / grid_->elemlist[iL].Vc ) * simdata_->visc_dt_factor);
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

        rho_l = rho;

        rho = Qr_[0];
        u = Qr_[1]/Qr_[0];
        v = Qr_[2]/Qr_[0];
        E = Qr_[3];

        p = (gama_-1.) * ( E - 0.5 * rho * ( pow(u,2) + pow(v,2) ) );

        c += sqrt(gama_*p/rho);

        Vn += (u*nx + v*ny);

        Vn = fabs(Vn/2.);
        c = c/2.;

        // Viscous contribution terms:
        rho_r = rho; rho = 0.5 * (rho_l +rho_r);
        Bv1 = (gama_/rho), Bv2 = 4.0/ ( 3.0 * rho) ;
        if(Bv1<Bv2) Bv1 =Bv2;

        Bv1 = Bv1 * gasdata_->Mu_inf * Sf * Sf / gasdata_->Prndtl;

        cell_radii[iL] += ( ((Vn+c)*Sf*0.5) +  (Bv1 / grid_->elemlist[iL].Vc ) * simdata_->visc_dt_factor );
        cell_radii[iR] += ( ((Vn+c)*Sf*0.5) +  (Bv1 / grid_->elemlist[iR].Vc ) * simdata_->visc_dt_factor );
    }

    for(i=0; i<Nelem; i++) {

        cell_radii[i] = cell_radii[i] / (grid_->elemlist[i].Vc);

        dt_cell_[i] = CFL_ / cell_radii[i];
    }

    emptyarray(cell_radii);

    return;
}
















