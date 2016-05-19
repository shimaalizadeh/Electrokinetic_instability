#ifndef TABLEREADING_H
#define TABLEREADING_H

#include<iostream>
#include<fstream>
#include<cmath>
#include<stdio.h>

#include"Friendly_Vector.h"

using namespace std;

class Tablereading{
    
public:
    
    Tablereading(const char* filename="sigma_lambda_star.in", unsigned gx=12, unsigned gy=12, double dy=1./10., double Pe=0.7, double gp=-1./3., unsigned P_BC_L=1, unsigned P_BC_R=1, unsigned Phi_BC_L=0, unsigned Phi_BC_R=0, double Phi_L=0., double Phi_R=0., double P_L=0., double P_R=0.,  unsigned max_sigma_star_counter=48, unsigned max_lambda_star_counter=38, double sigma_star_min=-0.1, double lambda_star_min=0.01);    
   
    ~Tablereading();
    
    int    gx()                  const  {return gx_;}
    int    gy()                  const  {return gy_;} 
    double f_bar(int i,int j)    const  {return f_bar_(i,j);}
    double ge_bar(int i,int j)   const  {return ge_bar_(i,j);}
    double gc_bar(int i,int j)   const  {return gc_bar_(i,j);}
    double gp_m(int i,int j)     const  {return gp_m_(i,j);}
    double ge_m(int i,int j)     const  {return ge_m_(i,j);}
    double gc_m(int i,int j)     const  {return gc_m_(i,j);}
    double gp_p(int i,int j)     const  {return gp_p_(i,j);}
    double ge_p(int i,int j)     const  {return ge_p_(i,j);}
    double gc_p(int i,int j)     const  {return gc_p_(i,j);}
  
    void find_from_table(Friendly_Vector & sigma, Friendly_Vector & lambda, Friendly_Vector & C_bar);
    
    void calculate_A_B_f_RHS(Friendly_Vector &C_bar, Friendly_Vector &C0, Friendly_Vector &Cs, Friendly_Vector &lambda, Friendly_Vector &nondimensional_h, double *x_center, double *x_face);

    void calculate_dC(Friendly_Vector &C_bar, Friendly_Vector &C0, Friendly_Vector &P, Friendly_Vector &Phi,Friendly_Vector &nondimensional_h, Friendly_Vector &U, Friendly_Vector &V, Friendly_Vector &k_, Friendly_Vector &dC, double *x_center, double *x_face, double dt, double b);

    ////some vectors required for future calculations
    Friendly_Vector A1_x_;
    Friendly_Vector A2_x_;
    Friendly_Vector B1_x_;
    Friendly_Vector B2_x_;
    Friendly_Vector f1_x_;
    Friendly_Vector f2_x_;
    Friendly_Vector A1_y_;
    Friendly_Vector A2_y_;
    Friendly_Vector B1_y_;
    Friendly_Vector B2_y_;
    Friendly_Vector f1_y_;
    Friendly_Vector f2_y_;
    Friendly_Vector RHS_1_;
    Friendly_Vector RHS_2_;



private:
    
    Friendly_Vector f_bar_;
    Friendly_Vector ge_bar_;
    Friendly_Vector gc_bar_;
    Friendly_Vector gp_m_;
    Friendly_Vector ge_m_;
    Friendly_Vector gc_m_;
    Friendly_Vector gp_p_;
    Friendly_Vector ge_p_;
    Friendly_Vector gc_p_;

    Friendly_Vector Gp_m_;
    Friendly_Vector Ge_m_;
    Friendly_Vector Gc_m_;

    Friendly_Vector A1_;
    Friendly_Vector A2_;
    Friendly_Vector B1_;
    Friendly_Vector B2_;
    Friendly_Vector f1_;
    Friendly_Vector f2_;
     
    ///used to calculate U and V
    Friendly_Vector a1_x_;
    Friendly_Vector b1_x_;
    Friendly_Vector c1_x_;
    Friendly_Vector a1_y_;
    Friendly_Vector b1_y_;
    Friendly_Vector c1_y_;

    Friendly_Vector adv_x_;
    Friendly_Vector diff_x_;
    Friendly_Vector elec_mig_x_;
    Friendly_Vector adv_y_;
    Friendly_Vector diff_y_;
    Friendly_Vector elec_mig_y_;

    Friendly_Vector Flux_x_;
    Friendly_Vector Flux_y_;
    
    double *table_sigma_star;
    double *table_lambda_star;
    double *table_f_bar;
    double *table_Psi_center;
    double *table_ge_bar;
    double *table_gc_bar;
    double *table_gp_m;
    double *table_ge_m;
    double *table_gc_m;
    double *table_gp_p;
    double *table_ge_p;
    double *table_gc_p;
    
    int     sigma_power_;
    int     lambda_power_;
    double  sigma_star_min_;
    double  lambda_star_min_;
    int     max_sigma_star_counter_;
    int     max_lambda_star_counter_;
    int     gx_;
    int     gy_; //local_size+2
    double  dy_;
    double  Pe_;
    double  gp_;
    int     P_BC_Type_L_ , P_BC_Type_R_;
    int     Phi_BC_Type_L_ , Phi_BC_Type_R_;
    double  Phi_L_ , Phi_R_;
    double  P_L_ , P_R_;
    double  aa;
    double  bb;
    double  dum1;
    double  dum2;

    //prevent copying and assignment
    Tablereading(const Tablereading &);
    Tablereading & operator=(const Tablereading &);

};

#endif
