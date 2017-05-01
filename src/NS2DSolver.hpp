#ifndef NS2DSOLVER_H
#define NS2DSOLVER_H


#include"SpaceSolver.hpp"


class NS2DSolver:public SpaceSolver{

public:

//  Construction Functions :
    NS2DSolver(void){}
    virtual ~NS2DSolver();

    virtual void setup_solver(MeshData*& meshdata_, SimData& simdata_
                             , GasProb& gasdata_);
    virtual void InitSol();

    virtual void UpdateResid(double **Resid_, double **Qn_);
    virtual void ComputeError(){}
    virtual void Compute_vertex_sol(const int& oiter);
    virtual void UpdateSolution(double **Qn_){}

    virtual void Compute_local_TimeStep(double* dt_cell_);

    virtual void evaluate_sol(double &Xp, double &Yp
                              , const int& eID, double *qq_);

    virtual void compute_GhostSol_WallBC(const double& nx, const double& ny
                                         , double* Ql, double* Qr);
    virtual void compute_GhostSol_SymmetryBC(const double& nx, const double& ny
                                             , double* Ql, double* Qr); // Same as Wall B.C. for Inviscid flows

    virtual void compute_GhostSol_farfieldBC(const double& nx, const double& ny,
                                     double* Ql, double* Qr);

    virtual void Compute_Inlet_charBC(const double& nx, const double& ny,
                                      double* Ql, double* Qr);
    virtual void Compute_Exit_charBC(const double& nx, const double& ny,
                                     double* Ql, double* Qr);

    virtual void Residulas_setZero(double** Residuals);

    virtual void SetGhostVariables();   // Setting Boundary Conditions

    virtual void Reconstruct_sol();

    virtual void Compute_left_right_facesol(const int& fID, double* Ql_
                                    , double* Qr_);
    virtual void Compute_left_right_boundfacesol(const int& fID, double* Ql_
                                    , double* Qr_);

    virtual void Compute_common_inviscidflux(double* Ql, double* Qr
                                     ,const double& nx, const double& ny
                                     ,double* flux_);

    virtual void Rusanov_flux(double *Ql, double *Qr
                      , const double& nx, const double& ny
                      , double* flux_);
    virtual void Roe_flux(){}

    virtual void compute_normal_inViscidFlux(const double& nx, const double& ny
                                     , double *Q_, double *normInvflux_
                                     , double& Vn, double& c);

protected:

    void Reset_solver();

    void Compute_common_viscousflux(const int &fID, const double *Ql
                                    , const double *Qr, const double *Qf_
                                    , double **dQl_,double **dQr_, double *flux_);

    void Compute_common_bound_viscousflux(const int &fID, const double *Ql
                                                , const double *Qr, const double *Qf_
                                                , double **dQl_,double **dQr_, double *flux_);

    void Compute_common_facesol(const double *Ql_, const double *Qr_, double *Qf_);

    void Compute_common_face_gradient(const int &fID, const double *Ql, const double *Qr
                                      ,double **dQl_,double **dQr_, double **dQf_);

    void compute_resid_OneFace(const double &Area_face
                                    ,double* inv_flux_, double* visc_flux_, double* resid_);

protected:
    double** visc_flux_com=nullptr;


};

#endif
