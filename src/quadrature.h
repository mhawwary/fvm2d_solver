#ifndef QUADRATURE_H
#define QUADRATURE_H

#include"general_tools.h"

struct GaussQuad {

public:

    unsigned int Nq=1;  // order of quadrature which is accurate for 2(Nq-1)

    double *Gaus_pts=nullptr;  // Gauss eta=xi Points
    //double *Gaus_X=nullptr;    // Gauss X coord
    double *Gaus_wts=nullptr;  // Gauss weights

public:
    //virtual Face& operator =(const Face& RFace)
    //GausQuad(const int& Nq_);
    //~GausQuad(void);

    void setup_quadrature(const int &Nq_){

        Nq = Nq_;

        Gaus_pts = new double[Nq+1];
        Gaus_wts = new double[Nq+1];

        //Gaus_X  = new double[Nq+1];

        define_gauss_quadrature();

        return;

    }

    void Reset_quad(){

        emptyarray(Gaus_pts);
        emptyarray(Gaus_wts);

        return;
    }

protected:
    void define_gauss_quadrature(){

        if(Nq==1){
            Gaus_pts[0] = 0.0 ;

            Gaus_wts[0] = 2.0 ;

        }else if(Nq==2){
            Gaus_pts[0] = -sqrt(1.0/3) ;
            Gaus_pts[1] =  sqrt(1.0/3) ;

            Gaus_wts[0] =  1.0 ;
            Gaus_wts[1] =  1.0 ;

        }else if(Nq==3){
            Gaus_pts[0] = -sqrt(3.0/5) ;
            Gaus_pts[1] = 0.0  ;
            Gaus_pts[2] = sqrt(3.0/5) ;

            Gaus_wts[0] =  5.0/9 ;
            Gaus_wts[1] =  8.0/9 ;
            Gaus_wts[2] =  5.0/9 ;

        }else if(Nq==4){
            Gaus_pts[0] =  -sqrt(3.0/7+2.0/7*sqrt(6.0/5)) ;
            Gaus_pts[1] =  -sqrt(3.0/7-2.0/7*sqrt(6.0/5)) ;
            Gaus_pts[2] =   sqrt(3.0/7-2.0/7*sqrt(6.0/5)) ;
            Gaus_pts[3] =   sqrt(3.0/7+2.0/7*sqrt(6.0/5)) ;

            Gaus_wts[0] =  (18-sqrt(30))/36.0  ;
            Gaus_wts[1] =  (18+sqrt(30))/36.0  ;
            Gaus_wts[2] =  (18+sqrt(30))/36.0  ;
            Gaus_wts[3] =  (18-sqrt(30))/36.0  ;

        }else if(Nq==5){
            Gaus_pts[0] =  -(1.0/3)*sqrt(5+2*sqrt(10.0/7)) ;
            Gaus_pts[1] =  -(1.0/3)*sqrt(5-2*sqrt(10.0/7)) ;
            Gaus_pts[2] =             0              ;
            Gaus_pts[3] =   (1.0/3)*sqrt(5-2*sqrt(10.0/7)) ;
            Gaus_pts[4] =   (1.0/3)*sqrt(5+2*sqrt(10.0/7)) ;

            Gaus_wts[0] =  (322-13*sqrt(70))/900.0 ;
            Gaus_wts[1] =  (322+13*sqrt(70))/900.0 ;
            Gaus_wts[2] =  128.0/225 ;
            Gaus_wts[3] =  (322+13*sqrt(70))/900.0 ;
            Gaus_wts[4] =  (322-13*sqrt(70))/900.0;

        }else{

            std::cout<< "\n Gauss Quadrature with order of:  "<<Nq<<"  is not implemented \n";
        }

        return;
    }
};

#endif

