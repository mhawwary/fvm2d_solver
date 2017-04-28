
#ifndef EULER2DSOLVER_H
#define EULER2DSOLVER_H

//#include"MeshData.h"
//#include"SimData.hpp"
//#include"quadrature.h"
//#include"general_tools.h"
//#include"global_var.h"

#include"SpaceSolver.hpp"


class Euler2DSolver:public SpaceSolver{

public:

//  Construction Functions :
    Euler2DSolver(void){}
    virtual ~Euler2DSolver();

    virtual void setup_solver(MeshData*& meshdata_, SimData& simdata_
                             , GasProb& gasdata_);
    virtual void InitSol();

    virtual void UpdateResid(double **Resid_, double **Qn_);
    virtual void ComputeError(){}
    virtual void Compute_vertex_sol();
    virtual void UpdateSolution(double **Qn_){}

    virtual void CalclocalTimeStep(double* dt_cell_);

    virtual void evaluate_sol(double &Xp, double &Yp
                              , const int& eID, double *qq_);

protected:

    void compute_resid_OneFace(const double &Area_face
                               ,double* face_flux_, double* resid_);

    void SetGhostVariables();   // Setting Boundary Conditions

    void Reconstruct_sol();

    void Compute_left_right_facesol(const int& fID, double* Ql_
                                    , double* Qr_);
    void Compute_left_right_boundfacesol(const int& fID, double* Ql_
                                    , double* Qr_);

    void Compute_common_inviscidflux(double* Ql, double* Qr
                                     ,const double& nx, const double& ny
                                     ,double* flux_);

    void Rusanov_flux(double *Ql, double *Qr
                      , const double& nx, const double& ny
                      , double* flux_);
    void Roe_flux(){}

    void compute_normal_inViscidFlux(const double& nx, const double& ny
                                     , double *Q_, double *normInvflux_
                                     , double& Vn, double& c);

    void compute_GhostSol_inviscidWallBC(const double& nx, const double& ny
                                         , double* Ql, double* Qr);

    void compute_GhostSol_farfieldBC(const double& nx, const double& ny,
                                     double* Ql, double* Qr);

    void Residulas_setZero(double** Residuals);

    void Reset_solver();

//protected:


};

#endif
