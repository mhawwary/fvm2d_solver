#include"SpaceSolver.hpp"
#include"SimData.hpp"
#include"general_tools.h"
#include"solver_tools.h"

class ExplicitTimeSolver{

public:
    ExplicitTimeSolver(void){}
    ~ExplicitTimeSolver(void);
    void setupTimeSolver(SpaceSolver* ospace_solver_, SimData* osimdata_, MeshData* meshdata);
    void SolveOneStep(double **qn_);

    int GetIter(){
        return IterNo_;
    }

    double GetTimeStep(){
        return dt_;
    }

    double GetPhyTime(){
        return phy_time_;
    }

    void SetIter(const int& iter_){
        IterNo_ = iter_;
        return;
    }

    void SetPhyTime(double& time_){
        phy_time_= time_;
        return;
    }

    void SetTimeStep(double& odt_){
        dt_ = odt_;
        return;
    }

    void   ComputeInitialResid(double** Qn_);
    double GetResNorm();
    double GetContinuityResNorm();
    double GetMomentumXResNorm();
    double GetMomentumYResNorm();
    double GetEnergyResNorm();

    double getdt_min(){ return dt_min; }
    double getdt_max(){ return dt_max; }

    void update_local_timestep();

protected:

    void FwdEuler(double **q_);
    void SSPRK22(double **q_);
    void SSPRK33(double **q_);

    void CopyOldSol(double **q_t_, double **qn_);

    void compute_dt_minmax();

    void Reset_time_solver();

protected:
    int Nelem=1;
    int Nelem_extend=1;
    int Nfaces=1;

    double *dt_elem_=nullptr;
    double **resid=nullptr;
    double **q_temp=nullptr;

    /* The space_solver_ is a pointer passed to the time solver
    * and is not supposed to be freed in this scope
    */

    SpaceSolver *space_solver_=nullptr;
    SimData     *simdata_=nullptr;
    MeshData *meshdata_=nullptr;

    int Ndof=4;

    int IterNo_=0;

    double dt_=0.0;

    double dt_min=1e5,dt_max=-1.0;

    double phy_time_=0.0;

    int RK_order_=1;

};
