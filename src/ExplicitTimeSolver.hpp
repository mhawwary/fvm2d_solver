#include"SpaceSolver.hpp"
#include"SimData.hpp"
#include"general_tools.h"

class ExplicitTimeSolver{

public:
    ExplicitTimeSolver(void){}
    ~ExplicitTimeSolver(void);
    void setupTimeSolver(SpaceSolver* ospace_solver_, SimData* osimdata_);
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

protected:

    void FwdEuler(double **q_);
    void SSPRK22(double **q_);
    //void SSPRK33(double **q_);

    void CopyOldSol(double **q_t_, double **qn_);

    //void UpdateIter();

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

    unsigned int Ndof=1;

    int IterNo_=0;

    double dt_=0.0;

    double phy_time_=0.0;

    int RK_order_=1;

};
