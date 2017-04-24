
#include"ExplicitTimeSolver.hpp"

ExplicitTimeSolver::~ExplicitTimeSolver(void){
    Reset_time_solver();
    return;
}

void ExplicitTimeSolver::setupTimeSolver(SpaceSolver *ospace_solver_
                                         , SimData *osimdata_){

    space_solver_=ospace_solver_;
    simdata_=osimdata_;

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

    dt_elem_ = new double[Nelem];

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

    space_solver_->UpdateResid(resid,qn_);

    switch (RK_order_) {

    case 1:

        FwdEuler(qn_);

        break;

    case 2:

        SSPRK22(qn_);

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

    unsigned int k;

    CopyOldSol(q_temp,q_);  // Copying level n solution and saving it

    // Step1:
    //-----------
    for(j=0; j<Nelem; j++)
        for(k=0; k<Ndof; k++){
            q_[j][k] = q_temp[j][k] + dt_ * resid[j][k];
        }

    //space_solver->UpdateResid(resid,q_);

    // Step2:
    //------------
    for(j=0; j<Nelem; j++)
        for(k=0; k<Ndof; k++) {
            q_[j][k] = 0.5 * ( q_temp[j][k] +  q_[j][k]
                                      + dt_ * resid[j][k] );
        }

    return;
}

void ExplicitTimeSolver::CopyOldSol(double **q_t_, double **qn_){

    register int j;

    unsigned int k;

    for(j=0; j<Nelem; j++)
        for(k=0; k<Ndof; k++)
            q_t_[j][k] = qn_[j][k];

    return;
}

double ExplicitTimeSolver::GetResNorm(){

    double resid_norm_=0.;

    resid_norm_ = L2norm(Nelem,Ndof,resid);

    return resid_norm_;
}





















