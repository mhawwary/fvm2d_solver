
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
   virtual void UpdateResid(double **Resid_, double **Qn_){}
//   virtual void UpdateSolution(double **Qn_)=0;
//   virtual void ComputeError()=0;
//   virtual void Compute_vertex_sol()=0;
//   virtual void Compute_exact_sol()=0;
//   virtual void UpdatePhyTime(const double& dt_)=0;

//   virtual void SetPhyTime(const double& time_)=0;
//   virtual double GetPhyTime()=0;
//   virtual double GetTimeStep()=0;
//   virtual double GetCFL()=0;

   virtual double** GetNumSolution(){
        return Qn;
    }

//   virtual double* GetVertexNumSol()=0;
//   virtual double* GetExactSolution()=0;

//   virtual void print_num_vertex_sol()=0;
//   virtual void print_exact_sol()=0;
//   virtual void print_exact_average_sol()=0;
//   virtual void print_num_average_sol()=0;

    virtual void CalcTimeStep();

protected:
//   virtual void UpdateResidOneCell(const int& cellid, double* q_, double* resid_);
//   virtual void Reconstruct_sol(const int& face_id, double* ql, double* qr);
//   virtual double Compute_common_flux(const double* ql, const double* qr
//                                      , const double& nx, const double& ny);
//   virtual void Compute_flux_upw();
//   virtual void get_left_right_sol();
//   virtual void Rusanov_flux();
//   virtual void Roe_flux();

   //virtual void CalcTimeStep();
    //virtual void CalcLocalTimeStep();


   void Reset_solver();


//protected:


};

#endif
