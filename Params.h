#pragma once

#include<iostream>
#include<fstream>
#include<math.h>
#include<cmath>
#include<stdio.h>

using namespace std;

#define pi 3.141592653589793

class Params{
    
public:
    
    Params(); //Use some default values
    
    ~Params(); //Destructor
    
    void   set_const_parameters();
    
    void   set_ny_dy(int ny) {ny_=ny;
                              dy_=Ly_/ny_;}
   
    double set_dt() {
                       double dt_x=delta_min*delta_min/4.;
                       double dt_y=dy_*dy_/4.;
                       if(dt_x<dt_y)
                          return(dt_x);
                       else
                          return(dt_y);}

    void   set_C_bar_L(double cbar_l) {C_bar_L_=cbar_l;}

    void   set_C_bar_R(double cbar_r) {C_bar_R_=cbar_r;}
    
    void   set_grid_points();
    
    double calculate_sigma_star(double h, double sigma) {return (h/2)*z*e*sigma/(Kb*T*epsilon0*epsilon);}
    
    double g_p_bar()    const {return g_p_bar_;}
    double Lambda_ref() const {return Lambda_ref_;}
    double Pe()         const {return Pe_;}
    double Lx()         const {return Lx_;}
    double Ly()         const {return Ly_;}
    int    nx()         const {return nx_;}
    int    ny()         const {return ny_;}
    int    nz()         const {return nz_;}
    int    gx()         const {return gx_;}
    int    N()          const {return N_;}
    double dy()         const {return dy_;}
    double dz()         const {return dz_;}
    double beta()       const {return beta_;}
    double dx_max()     const {return delta_max;}
    double dx_min()     const {return delta_min;}
    double a(int i)     const {return a_[i];}
    double b(int i)     const {return b_[i];}
    double C_bar_L()    const {return C_bar_L_;}
    double C_bar_R()    const {return C_bar_R_;}

    double *x_center;
    double *x_face;
    double dt;   //time step
    double Tmax;
    double time;
    int    period;
    int    WriteTimeCounter;
    int    confirmation;
    int    restart;  //to show if there is any restart file or not
 
private:
    
    //Electrochemical Parameters
    int    z;                 //valence
    double e;                 //elementary charge(Coulombs!)
    double Kb;                //Boltzmann Constant
    double epsilon0;          //permitivity constant
    double epsilon;           //dielectric relative permitivity
    double mu;                //dynamic viscosity(Pa.s)
    double T;                 //temperature(K)
    double Dif;               //mass Diffusivity (m^2/s),
    double C_ref;             //reference concentration(Molar)
    double g_p_bar_;          //Average Pressure driven flow in a slit
    
    double Lambda_ref_;       //Debye Length
    double Pe_;               //Peclet number
    
    //Problem Parameters
    double Lx_;                //domain length in x
    double Ly_;                //domain length in y
    double Lz_;                //domain length in z
    int    nx_;                //number of cells in x
    int    ny_;                //number of cells in y
    int    nz_;                //number of cells in z
    int     N_;                //total number of computational cells
    int    gx_;                //number of cells in x including the ghost cells
    int    gy_;                //number of cells in y including the ghost cells
    int    number_of_ghost;    //number of ghost cells from each side
    double dy_;                //size of cells in y
    double dz_;                //size of cells in z
    double C_bar_L_;           //left reservoir C_bar (constant)
    double C_bar_R_;           //right reservoir C_bar (constant)
    
    //Non-uniform mesh parameters
    double beta_;
    double ratio_;              //ratio=delta_max/delta_min (we keep it constant)
    double delta_max;
    double delta_min;
    int    refinement;
 
    //RK4 parameters
    double a_[4];
    double b_[4];

};

Params::Params(){
    
    refinement=16;

    beta_=0.1/double(refinement);
    ratio_=10;
    delta_max=0.4/double(refinement);
    delta_min=delta_max/ratio_;


    Lx_=1.;
    Ly_=pi/2.;
    Lz_=0.02;
    
    nx_=400;
    ny_=400;
    nz_=2;
    N_=nx_*ny_;
    
    number_of_ghost=1; 
    
    gx_=nx_+2*number_of_ghost;
    gy_=ny_+2*number_of_ghost;
    
    dy_=Ly_/ny_;
    dz_=Lz_/nz_;
 
    //determining time step for RK4:
    dt=set_dt();
    //dt=dt/2.;
     
    period=400; 
    Tmax=0.2;
    
    x_center =new double[gx_];
    x_face   =new double[gx_-1];

    a_[0]=0. , a_[1]=0.5 , a_[2]=0.5 , a_[3]=1.;
    b_[0]=1./6. , b_[1]=1./3. , b_[2]=1./3. , b_[3]=1./6. ;    
}

Params::~Params(){
    
    delete [] x_center;
    delete [] x_face;
}

void Params::set_const_parameters(){
    
    z=1;                 //valence
    e=1.602e-19;         //elementary charge(Coulombs!)
    Kb=1.38e-23;         //Boltzmann Constant
    epsilon0=8.85e-12;   //permitivity constant
    epsilon=80;          //dielectric relative permitivity
    mu=8.9e-4;           //dynamic viscosity(Pa.s)
    T=300;               //temperature(K)
    Dif=3e-9;            //mass Diffusivity (m^2/s),
    C_ref=0.001;         //reference concentration(Molar)
    g_p_bar_=-1./3.;      //Average Pressure driven flow in a slit
    
    Lambda_ref_=sqrt(epsilon0*epsilon*Kb*T/(2*z*e*e*C_ref*6.02214129e23*1000));
    Pe_=epsilon0*epsilon/(mu*Dif)*(Kb*T/(z*e))*(Kb*T/(z*e));

    ////Note:Manually change Pe and lambda
    Lambda_ref_=5e-7*0.05;
    Pe_=0.177091 ;
    //Pe_=10.;
}

void Params::set_grid_points(){
   
    double ii;
    //x on faces
    for(int i=0;i<nx_+1;i++)
    {
      x_face[i]=i*(Lx_/nx_);
      //x_face[i]=delta_max/beta_*log(delta_min/delta_max*(exp(beta_*i)-1)+1);
    }

    //x on cell center
    for(int i=-1;i<nx_+1;i++)
    {
      x_center[i+1]=(i+0.5)*(Lx_/nx_);
      //ii=(i+i+1)/2.;

      //x_center[i+1]=delta_max/beta_*log(delta_min/delta_max*(exp(beta_*ii)-1)+1);
    }
return; 
}
