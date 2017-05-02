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
    virtual void compute_wall_pressure_dist();

    virtual void Get_wall_stress_tensor(double* pp, double* otau_xx_
                                        , double* otau_xy_, double* otau_yy_ ){
        pp = p_wall_;
        otau_xx_ = tau_xx_wall_;
        otau_xy_ = tau_xy_wall_;
        otau_yy_ = tau_yy_wall_;

        return;
    }

    virtual double* Get_tau_xx(){
        return tau_xx_wall_;
    }

    virtual double* Get_tau_xy(){
        return tau_xy_wall_;
    }

    virtual double* Get_tau_yy(){
        return tau_yy_wall_;
    }

    virtual double* Get_wall_skin_friction(){
        return Cf_;
    }

protected:

    void Reset_solver();

    void evaluate_sol_visc(double &Xp, double &Yp, const int& eID, double *qq_);
    void Compute_common_viscousflux(const int &fID, const double *Ql
                                    , const double *Qr, const double *Qf_
                                    , double **dQl_,double **dQr_, double *flux_);

    void Compute_common_bound_viscousflux(const int &fID, const double *Ql
                                                , const double *Qr, const double *Qf_
                                                , double *flux_);

    void Compute_common_facesol(const double *Ql_, const double *Qr_, double *Qf_);

    void Compute_common_face_gradient(const int &fID, const double *Ql, const double *Qr
                                      ,double **dQl_,double **dQr_, double **dQf_);

    void Compute_common_bound_face_gradient(const int &fID, const double *Ql
                                            , const double *Qr, double **dQf_);

    void compute_resid_OneFace(const double &Area_face
                                    ,double* inv_flux_, double* visc_flux_, double* resid_);

    void update_wall_fields(){
        compute_wallShearStress();
        compute_wall_pressure_dist();
        return;
    }
    void compute_wall_velocity_gradients();
    void compute_wallShearStress();

protected:

    double** visc_flux_com=nullptr;

    double* dudx_wall=nullptr;
    double* dudy_wall=nullptr;
    double* dvdx_wall=nullptr;
    double* dvdy_wall=nullptr;

    double* wall_node_nx=nullptr;
    double* wall_node_ny=nullptr;

    double *tau_t_wall_=nullptr;  // tangent shear stress at the wall
    double *Cf_=nullptr;
    double *tau_xx_wall_ = nullptr;
    double *tau_yy_wall_ = nullptr;
    double *tau_xy_wall_ = nullptr;

};

#endif
