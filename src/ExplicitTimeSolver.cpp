
#include"ExplicitTimeSolver.hpp"

ExplicitTimeSolver::~ExplicitTimeSolver(void){
    Reset_time_solver();
    return;
}

void ExplicitTimeSolver::setupTimeSolver(SpaceSolver *ospace_solver_
                                         , SimData *osimdata_, MeshData* omeshdata){

    space_solver_=ospace_solver_;
    simdata_=osimdata_;
    meshdata_ = omeshdata;

    RK_order_= simdata_->RK_order;
    Ndof = space_solver_->GetNdof();
    Nelem = space_solver_->GetNelem();
    Nelem_extend = space_solver_->GetExtended_Nelem();
    Nfaces = space_solver_->GetNfaces();

    resid = new double* [Nelem];

    register int i;

    for(i=0; i<Nelem; i++){

        resid[i] = new double[Ndof];
    }

    dt_ = space_solver_->GetTimeStep();

    dt_elem_ =  new double [Nelem];

    if(simdata_->use_local_timeStep==1){
        space_solver_->Compute_local_TimeStep(dt_elem_);
    }else{
        for(i=0; i<Nelem; i++)
            dt_elem_[i]=dt_;
    }


    if(RK_order_>1){

        q_temp = new double* [Nelem];

        for(i=0; i<Nelem; i++){

            q_temp[i] = new double[Ndof];
        }
    }

    return;
}

void ExplicitTimeSolver::Reset_time_solver(){

    emptyarray(Nelem,resid);
    emptyarray(Nelem,q_temp);
    emptyarray(dt_elem_);

    return;
}

void ExplicitTimeSolver::SolveOneStep(double **qn_){

    switch (RK_order_) {

    case 1:

        FwdEuler(qn_);

        break;

    case 2:

        SSPRK22(qn_);

        break;

    case 3:

        SSPRK33(qn_);

        break;

    default:

        char *ss=nullptr; ss= new char[100];
        sprintf(ss,"RK order of %d ", RK_order_);
        _notImplemented(ss);
        emptyarray(ss);
        break;
    }

    IterNo_++;

    return;
}

void ExplicitTimeSolver::FwdEuler(double **q_){
    return;
}

void ExplicitTimeSolver::SSPRK22(double **q_){

    register int j;

    int k;

    CopyOldSol(q_temp,q_);  // Copying level n solution and saving it

    // Step1:
    //-----------
    for(j=0; j<Nelem; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] = q_temp[j][k] + (dt_elem_[j] * resid[j][k]/meshdata_->elemlist[j].Vc);

    space_solver_->UpdateResid(resid,q_);

    // Step2:
    //------------
    for(j=0; j<Nelem; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] = 0.5 * ( q_temp[j][k] +  q_[j][k]
                                      + dt_elem_[j] * (resid[j][k]/meshdata_->elemlist[j].Vc) );

    space_solver_->UpdateResid(resid,q_);

    return;
}

void ExplicitTimeSolver::SSPRK33(double **q_){

    register int j;

    int k;

    CopyOldSol(q_temp,q_);  // Copying level n solution and saving it


    // Step1:
    //-----------
    for(j=0; j<Nelem; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] = q_temp[j][k] +  dt_elem_[j] * (resid[j][k]/meshdata_->elemlist[j].Vc);


    space_solver_->UpdateResid(resid,q_);

    // Step2:
    //------------
    for(j=0; j<Nelem; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] =  (0.75 * q_temp[j][k] )
                    + 0.25 * ( q_[j][k] +  dt_elem_[j] * (resid[j][k]/meshdata_->elemlist[j].Vc) );

    space_solver_->UpdateResid(resid,q_);

    // Step3:
    //--------------
    for(j=0; j<Nelem; j++)
        for(k=0; k<Ndof; k++)
            q_[j][k] =  ( q_temp[j][k]/3. )
                    + 2. * ( q_[j][k] +  dt_elem_[j] * (resid[j][k]/meshdata_->elemlist[j].Vc) ) / 3.;

    space_solver_->UpdateResid(resid,q_);

    return;
}

void ExplicitTimeSolver::CopyOldSol(double **q_t_, double **qn_){

    register int j;

    int k;

    for(j=0; j<Nelem; j++)
        for(k=0; k<Ndof; k++)
            q_t_[j][k] = qn_[j][k];

    return;
}

void ExplicitTimeSolver::compute_dt_minmax(){

    register int j;

    dt_min=1e5; dt_max=-1.0;

    for(j=0; j<Nelem; j++){
        if(dt_elem_[j]<dt_min) dt_min = dt_elem_[j];
        if(dt_elem_[j]>dt_max) dt_max = dt_elem_[j];
    }

    return;
}

void ExplicitTimeSolver::update_local_timestep(){

    space_solver_->Compute_local_TimeStep(dt_elem_);

    compute_dt_minmax();

    return;
}

double ExplicitTimeSolver::GetResNorm(){

    double resid_norm_=0.;

    int i=0;
    double *Vol=nullptr; Vol = new double[Nelem];
    for(i=0; i<Nelem; i++) Vol[i] = 1./meshdata_->elemlist[i].Vc;

    resid_norm_ = L2norm(Nelem,Ndof,resid,Vol);

    emptyarray(Vol);

    return resid_norm_;
}

double ExplicitTimeSolver::GetContinuityResNorm(){

    double resid_norm_=0.;

    double *resid_density = new double[Nelem];

    register int i;

    for(i=0; i<Nelem; i++) resid_density[i] = resid[i][0] / meshdata_->elemlist[i].Vc;

    resid_norm_ = L2norm_perdof(Nelem,resid_density);

    emptyarray(resid_density);

    return resid_norm_;
}

double ExplicitTimeSolver::GetMomentumXResNorm(){

    double resid_norm_=0.;

    double *resid_momentum = new double[Nelem];

    register int i;

    for(i=0; i<Nelem; i++) resid_momentum[i] = resid[i][1]/meshdata_->elemlist[i].Vc;

    resid_norm_ = L2norm_perdof(Nelem,resid_momentum);

    emptyarray(resid_momentum);

    return resid_norm_;
}

double ExplicitTimeSolver::GetMomentumYResNorm(){

    double resid_norm_=0.;

    double *resid_momentum = new double[Nelem];

    register int i;

    for(i=0; i<Nelem; i++) resid_momentum [i] = resid[i][2]/ meshdata_->elemlist[i].Vc;

    resid_norm_ = L2norm_perdof(Nelem,resid_momentum);

    emptyarray(resid_momentum);

    return resid_norm_;
}

double ExplicitTimeSolver::GetEnergyResNorm(){

    double resid_norm_=0.;

    double *resid_energy = new double[Nelem];

    register int i;

    for(i=0; i<Nelem; i++) resid_energy [i] = resid[i][3]/ meshdata_->elemlist[i].Vc;

    resid_norm_ = L2norm_perdof(Nelem,resid_energy);

    emptyarray(resid_energy);

    return resid_norm_;
}

void ExplicitTimeSolver::ComputeInitialResid(double** Qn_){

    space_solver_->UpdateResid(resid,Qn_);

    return ;
}













