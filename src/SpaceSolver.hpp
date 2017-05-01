
#ifndef SPACESOLVER_H
#define SPACESOLVER_H


#include"MeshData.h"
#include"SimData.hpp"
#include"quadrature.h"
#include"general_tools.h"
#include"global_var.h"

class SpaceSolver{

public:

//  Construction Functions :
   SpaceSolver(){}
   virtual ~SpaceSolver(){}

   virtual void setup_solver(MeshData*& meshdata_, SimData& simdata_
                             , GasProb& gasdata_)=0;
   virtual void InitSol()=0;
   virtual void UpdateResid(double **Resid_, double **Qn_)=0;
   virtual void ComputeError()=0;
   virtual void Compute_vertex_sol(const int& oiter)=0;
   virtual void UpdateSolution(double **Qn_)=0;
   virtual void Compute_local_TimeStep(double* dt_cell_)=0;

   double** GetNumSolution(){
       return Qc;
   }

   double** GetVertexSolution(){
       return Qv;
   }

   int GetNdof(){
       return Ndof;
   }

   int GetNelem(){
       return Nelem;
   }

   int GetExtended_Nelem(){
       return Nelem_extend;
   }

   int GetNfaces(){
       return Nfaces;
   }

   void SetIter(const int& oiter){
       iter_ = oiter;
       return;
   }

   int GetIter(void){
       return iter_;
   }

   void SetPhyTime(const double& time_){
       phy_time_ = time_;
   }

   double GetPhyTime(){
       return phy_time_;
   }

   double GetTimeStep(){

       return dt_;
   }

   double GetCFL(){

       return CFL_;
   }

protected:

   virtual void evaluate_sol(double &Xp, double &Yp
                             , const int& eID, double *qq_)=0;

   virtual void compute_GhostSol_WallBC(const double& nx, const double& ny
                                        , double* Ql, double* Qr)=0;
   virtual void compute_GhostSol_SymmetryBC(const double& nx, const double& ny
                                            , double* Ql, double* Qr)=0; // Same as Wall B.C. for Inviscid flows

   virtual void compute_GhostSol_farfieldBC(const double& nx, const double& ny,
                                    double* Ql, double* Qr)=0;

   virtual void Compute_Inlet_charBC(const double& nx, const double& ny,
                                     double* Ql, double* Qr)=0;
   virtual void Compute_Exit_charBC(const double& nx, const double& ny,
                                    double* Ql, double* Qr)=0;

   virtual void Residulas_setZero(double** Residuals)=0;

   virtual void SetGhostVariables()=0;   // Setting Boundary Conditions

   virtual void Reconstruct_sol()=0;

   virtual void Compute_left_right_facesol(const int& fID, double* Ql_
                                   , double* Qr_)=0;
   virtual void Compute_left_right_boundfacesol(const int& fID, double* Ql_
                                   , double* Qr_)=0;

   virtual void Compute_common_inviscidflux(double* Ql, double* Qr
                                    ,const double& nx, const double& ny
                                    ,double* flux_)=0;

   virtual void Rusanov_flux(double *Ql, double *Qr
                     , const double& nx, const double& ny
                     , double* flux_)=0;
   virtual void Roe_flux()=0;

   virtual void compute_normal_inViscidFlux(const double& nx, const double& ny
                                    , double *Q_, double *normInvflux_
                                    , double& Vn, double& c)=0;

protected:
   /* These pointers are passed to the space solver
   * and is not supposed to be freed in this scope
   */
   MeshData *grid_=nullptr;
   SimData *simdata_=nullptr;
   GasProb *gasdata_ = nullptr;

   double phy_time_=0.0;
   double iter_=1e-5;
   double dt_=1e-5;
   double CFL_=0.5;

   int Nelem=1;
   int Nelem_extend=1;
   int Nfaces=1;
   int Nnodes=1;
   int Nbndfaces=1;
   int NbndElem=1;

   int scheme_order = 1;
   int Ndof = 4;

   /* Locally defined arrays and can be freed in this scope */

   double **Qc=nullptr;      // Nelem * Ndof long solution array ( at cell centers )

   double **dQdx=nullptr;
   double **dQdy=nullptr;

   double **Qv=nullptr;       // Solution at the vertices of cells, total no. of vertices long

   double **flux_com=nullptr;  // common interface flux, Nfaces long

   int comp_vertex_iter=0;

};

#endif
