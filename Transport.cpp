/*
  Discription:   This code solves height averaged Nernst-Planck transport equation for
                 negative ionic species in a 2-D domain confined with a selective membrane
                 at the end of x-direction.
 
  The governing eqn:   dC_bar/dt + dFx_bar/dx +dFy_bar/dy=0
     
                       Fx_bar=u_bar C_bar - dC_bar/dx + C_bar dPhi/dx

                       Fy_bar=v_bar C_bar - dC_bar/dy + C_bar dPhi/dy

  P Phi Distribution: Pressure and Potential are obtained by solving Poisson equation with 
  variable height using Hypre Linear Algebra Package.
 */

/*
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
 */

#include "mpi.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <cmath>
#include<stdio.h>

#include "Friendly_Vector.h"
#include "Params.h"
#include "TableReading.h"
#include "Communication.h"
#include "Hypre_Poisson_Solver.h"
#include "mpi_write.h"


using namespace std;

#define pi  3.141592653589793

double   sigma=-0.002; //-0.013;           //surface chage (Coulombs)
double   h_ref=1e-6;             //slit reference height (meter)
double   C0_L=0.01;            //Left reservoir nondim C0
double   C0_R=1.+C0_L;          //Right reservoir nondim C0
double   Phi_L=-100.;              //nondim Potential at left
double   Phi_R=0.;               //nondim Potential at right


//////////////////////////////////////////////////////////////////////

void Initialize(int myid, int num_procs, int y_local_size, Params & params, Friendly_Vector &h, Friendly_Vector &nondimensional_h, Friendly_Vector &sigma_star, Friendly_Vector &lambda,
                Friendly_Vector &Cs, Friendly_Vector &C0, Friendly_Vector &P,
                Friendly_Vector &Phi, Friendly_Vector &k_, Friendly_Vector &dC_bar){
 
    if(params.restart==1)
    {
       read_restart_file(params.gx(), params.ny(), y_local_size, 1 , num_procs, num_procs, myid, params.WriteTimeCounter, params.time, C0, P, Phi);
       
       if(myid==0){
       cout<<"final write-counter="<<params.WriteTimeCounter<<endl;
       cout<<"final write-time="<<params.time<<endl;
       }

        
    }

else{
       for(int j=1; j<y_local_size+1; j++)
       {
         for(int i=1; i<params.nx()+1; i++)
         {
                  C0(i,j)=1.;

                  P(i,j)=0.;

                  Phi(i,j)=0.;
         }
            ///to remove negative C0 in the ghost cell I do interpolation in log base
            C0(1,j)=C0_L;
            C0(0,j)= 2*log(C0_L) -log(C0(1,j));
            C0(0,j)= exp(C0(0,j));

            C0(params.nx()+1,j)= 2*log(C0_R) - log(C0(params.nx(),j)) ;
            C0(params.nx()+1,j)=exp(C0(params.nx()+1,j));

            Phi(0,j)=2*Phi_L - Phi(1,j);
            Phi(params.nx()+1,j)= 2*Phi_R - Phi(params.nx(),j);

            P(0,j)=0.;
            P(params.nx()+1,j)=0.;
            
          }

   }//End of else
        
       for(int j=1; j<y_local_size+1; j++)
        {
         for(int i=1; i<params.nx()+1; i++)
          {
            h(i,j)=1e-6;
            
            nondimensional_h(i,j)=h(i,j)/h_ref; 
            
            sigma_star(i,j)=params.calculate_sigma_star(h(i,j),sigma);
            
            lambda(i,j)=params.Lambda_ref()/(h(i,j)/2);

            //Note:remove later!!!!
            lambda(i,j)=0.001;
         
            Cs(i,j)=2*lambda(i,j)*lambda(i,j)*abs(sigma_star(i,j));
        
            //Note:remove later!!!!
            Cs(i,j)=0.01;
  
            dC_bar(i-1,j-1)=0.;

            k_(i-1,j-1)=0.;
        }
     
        //Boundary Conditions in x        
        //Note: We have assumed the same height and sigma_star for ghost cells as their adjacent cells

        h(0,j)=h(1,j);
        h(params.nx()+1,j)=h(params.nx(),j);

        nondimensional_h(0,j)=nondimensional_h(1,j);
        nondimensional_h(params.nx()+1,j)=nondimensional_h(params.nx(),j);

        lambda(0,j)=lambda(1,j);
        lambda(params.nx()+1,j)=lambda(params.nx(),j);

        sigma_star(0,j)= sigma_star(1,j);
        sigma_star(params.nx()+1,j)= sigma_star(params.nx(),j);

        Cs(0,j)=Cs(1,j);
        Cs(params.nx()+1,j)=Cs(params.nx(),j);
    }   

return;
}

void Add_Perturbation(Params &params, int domain_y_ilower,int local_size, double C0_L, double C0_R, Friendly_Vector &C0){

    double thickness=2*(params.x_center[1]-params.x_center[0]);
    double amplitude=5.*thickness;
    double location=0.1;
    double *y=new double[local_size];

    for(int j=0;j<local_size;j++)
    {
      y[j]=(domain_y_ilower+j+0.5)*params.dy();
    }

    for(int j=1;j<local_size+1;j++)
    {
      for(int i=1;i<params.nx()+1;i++)
      {
          C0(i,j)=1./2.*( 1+tanh( (params.x_center[i]-location-amplitude*sin(2*pi/params.Ly()*y[j-1]))/thickness ) ) + C0_L;
      
      }
      
          C0(0,j)= 2*log(C0_L) -log(C0(1,j));
          C0(0,j)= exp(C0(0,j));

          C0(params.nx()+1,j)= 2*log(C0_R) - log(C0(params.nx(),j)) ;
          C0(params.nx()+1,j)=exp(C0(params.nx()+1,j));
     }

delete [] y;

return;
}

void Calculate_Cbar_Initial_Time(int gx, int local_gy, Friendly_Vector &sigma_star, Friendly_Vector &lambda,Friendly_Vector &C0, Friendly_Vector &C_bar){

  Tablereading sigma_lambda0("sigma_lambda0.in", gx, local_gy);

  sigma_lambda0.find_from_table(sigma_star, lambda, C0);

  for(int j=0;j<local_gy;j++)
    {
      for(int i=0;i<gx;i++)
        {
          C_bar(i,j)=C0(i,j)*sigma_lambda0.f_bar(i,j);
        }
    }

  return;
}

void Calculate_C_star(Params & params, Communication & communicate, int y_local_size, int RK4_counter, Friendly_Vector & C_bar, Friendly_Vector & k_, Friendly_Vector & C_star){

       for(int j=1;j<y_local_size+1;j++)
       {
          for(int i=1;i<params.nx()+1;i++)
          {
             C_star(i,j)= C_bar(i,j) + params.a(RK4_counter)*k_(i-1,j-1);
          }

     C_star(0,j)= 2*log(params.C_bar_L()) -log(C_star(1,j));
     C_star(0,j)= exp(C_star(0,j));

     C_star(params.nx()+1,j)= 2*log(params.C_bar_R()) - log(C_star(params.nx(),j)) ;
     C_star(params.nx()+1,j)=exp(C_star(params.nx()+1,j));

       }

      //MPI_Barrier(MPI_COMM_WORLD);
      communicate.SendRecv(C_star);
return;
}



void Calculate_Coefficients_From_Table(Tablereading & sigma_lambda_star, double *x_center, double *x_face, Friendly_Vector &sigma_star, Friendly_Vector &lambda, Friendly_Vector &C_bar, Friendly_Vector &Cs, Friendly_Vector &C0, Friendly_Vector &nondimensional_h){

  sigma_lambda_star.find_from_table(sigma_star, lambda, C_bar);

  for(int j=0;j<sigma_lambda_star.gy();j++)
    {
      for(int i=0;i<sigma_lambda_star.gx();i++)
	{
	  C0(i,j)=C_bar(i,j)/sigma_lambda_star.f_bar(i,j);
	}
    }

  sigma_lambda_star.calculate_A_B_f_RHS(C_bar, C0, Cs, lambda, nondimensional_h, x_center, x_face);

  return;
}

void Solve(Tablereading &sig_lam_star, Communication &communicate, Poisson_Solver &solver, Friendly_Vector &C_star, Friendly_Vector &C0, Friendly_Vector &Cs, Friendly_Vector &nondimensional_h, Friendly_Vector &P, Friendly_Vector &Phi, Friendly_Vector &U, Friendly_Vector &V, Friendly_Vector &dC_bar, Friendly_Vector &k_, double *x_center, double *x_face, double dt, double b, int myid){
//double t1, t2;
  ////Solve for P and Phi
  solver.Hypre_Set_Value(sig_lam_star.A1_x_,sig_lam_star.B1_x_,sig_lam_star.A2_x_,sig_lam_star.B2_x_,sig_lam_star.A1_y_,sig_lam_star.B1_y_,sig_lam_star.A2_y_,sig_lam_star.B2_y_, sig_lam_star.RHS_1_, sig_lam_star.RHS_2_, P, Phi);
  
//t1=MPI_Wtime();
  solver.Hypre_Solve();
  solver.Hypre_Get_Solution(P, Phi);
//t2=MPI_Wtime();

/*
if(myid==0)
{
cout<<"elapse time is  "<<t2-t1<<endl;
}
*/
  //MPI_Barrier(MPI_COMM_WORLD);
  communicate.SendRecv(P);
  communicate.SendRecv(Phi);
  
  ////Calculate k
  sig_lam_star.calculate_dC(C_star, C0, P, Phi, nondimensional_h, U, V, k_, dC_bar, x_center, x_face, dt, b);
  

return;
}

void Update_C_bar(Params &params, int y_local_size, Friendly_Vector &C_bar, Friendly_Vector & dC_bar){

   for(int j=1;j<y_local_size+1;j++)
   {
     for(int i=1;i<params.nx()+1;i++)
     {
        C_bar(i,j)+=dC_bar(i-1,j-1);
        dC_bar(i-1,j-1)=0.;
     }

     ////Enforcing boundary Condition in x
     C_bar(0,j)= 2*log(params.C_bar_L()) -log(C_bar(1,j));
     C_bar(0,j)= exp(C_bar(0,j));

     C_bar(params.nx()+1,j)= 2*log(params.C_bar_R()) - log(C_bar(params.nx(),j)) ;
     C_bar(params.nx()+1,j)=exp(C_bar(params.nx()+1,j));

   }

return;
}

void export_matrix_rhs(Tablereading &sig_lam_star, Poisson_Solver & solver){

       solver.Hypre_Export(sig_lam_star.A1_x_,sig_lam_star.B1_x_,sig_lam_star.A2_x_,sig_lam_star.B2_x_,sig_lam_star.A1_y_,sig_lam_star.B1_y_,sig_lam_star.A2_y_,sig_lam_star.B2_y_, sig_lam_star.RHS_1_, sig_lam_star.RHS_2_ );

return;
}


int main(int argc, char *argv[])
{
    
    if (argc != 2) {
        cerr << "Please clarify if there is any restart file!" << endl;
        exit(1);
    }
      
    //////////////////////////Initialize MPI///////////////////////
    
    int myid, num_procs;   
  
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    ///////////////////////////////////////////////////////////////
    //cout<<setprecision(16);
    
    Params parameters; //Note:later add a file as an input.
    
    parameters.set_const_parameters();
    
    parameters.set_grid_points();
 
    ////// Preliminaries: want at least one processor per y ///////
    if (parameters.ny() < num_procs)
    {
        parameters.set_ny_dy(num_procs);
        parameters.dt=parameters.set_dt();
        parameters.dt/=4.;
    } 
    
    //parameters.dt=parameters.dt/2.;
      parameters.dt = 1e-6;

    if(myid==0)
    {
      cout<<"your simulation Tmax="<<parameters.Tmax<<"  and dt="<<parameters.dt<<endl;
      cout<<"to continue press a number, otherwise push ctrl+c!";
      cin>>parameters.confirmation;
    } 
    //////////////////////////Row Partitioning/////////////////////
    //
    //Each processor knows only elements in its own computational domain.
    //the y-range is denoted by domain_ilower and domain_iupper and x-range
    //is nx for all processors. Note that we partitioning is in y-direction
    // only.
    //
    
    int local_size = parameters.ny()/num_procs;
    
    int domain_y_ilower=(local_size)*myid;
    
    int domain_y_iupper = local_size*(myid+1);
    
        domain_y_iupper=  domain_y_iupper-1;

    ////////////////////////////////////////////////////////////////
    
     ////mesh generation
     Mesh_Generation(parameters.nx(), parameters.ny(), parameters.nz(), local_size, 1, num_procs, 1, num_procs, myid, domain_y_ilower, parameters.x_center, parameters.dy(), parameters.dz());

    ////Define Table reading object
    Tablereading sig_lam_star("sigma_lambda_star.in", parameters.gx(), local_size+2, parameters.dy(), parameters.Pe(), parameters.g_p_bar(), 1, 1, 0, 0, Phi_L, Phi_R);

    ////Define Communication object
    Communication communicate(num_procs, myid, parameters.nx(), local_size);
    
    ////Define Hypre_Poisson_Solver object
    Poisson_Solver solver(0, 1000, 1e-8, myid, domain_y_ilower, domain_y_iupper, local_size, parameters.nx(), parameters.ny(), parameters.dy(), 1, 1, 0, 0, parameters.x_center, parameters.x_face);

    solver.Hypre_Create();

    ////Required vectors
    Friendly_Vector   C_bar(2,parameters.gx(),local_size+2);
    Friendly_Vector   C0(2,parameters.gx(),local_size+2);
    Friendly_Vector   Cs(2,parameters.gx(),local_size+2);
    Friendly_Vector   sigma_star(2,parameters.gx(),local_size+2);
    Friendly_Vector   lambda(2,parameters.gx(),local_size+2);
    Friendly_Vector   h(2,parameters.gx(),local_size+2);
    Friendly_Vector   nondimensional_h(2,parameters.gx(),local_size+2);
    Friendly_Vector   P(2,parameters.gx(),local_size+2);
    Friendly_Vector   Phi(2,parameters.gx(),local_size+2);
    Friendly_Vector   C_star(2,parameters.gx(),local_size+2);
    Friendly_Vector   dC_bar(2,parameters.nx(),local_size);
    Friendly_Vector   k_(2,parameters.nx(),local_size);

    Friendly_Vector   U(2,parameters.nx()+1,local_size);
    Friendly_Vector   V(2,parameters.nx(),local_size+1);


    ////Initialization
    parameters.restart = atoi(argv[1]);
    Initialize(myid, num_procs, local_size , parameters, h, nondimensional_h, sigma_star, lambda,Cs, C0, P, Phi, k_, dC_bar);       

    ////Adding Perturbation
    if(parameters.restart==0){
    Add_Perturbation(parameters, domain_y_ilower, local_size, C0_L, C0_R, C0);

    parameters.time=0.;

    parameters.WriteTimeCounter=0;
    }

    ////Some Initial Communications
    communicate.SendRecv(h);
    communicate.SendRecv(sigma_star);
    communicate.SendRecv(lambda);
    communicate.SendRecv(C0);
    communicate.SendRecv(Cs);
    communicate.SendRecv(P);
    communicate.SendRecv(Phi);

    ////Finding C_bar for initial time using Table 
    Calculate_Cbar_Initial_Time(parameters.gx(), local_size+2, sigma_star, lambda, C0, C_bar);
    
    ////setting C_bar_L and C_bar_R and they remain constant in time
    parameters.set_C_bar_L( exp((log(C_bar(0,1))+log(C_bar(1,1)))/2. ) );
    parameters.set_C_bar_R( exp( (log(C_bar(parameters.nx(),1))+log(C_bar(parameters.nx()+1,1)))/2. ) );
 
if(parameters.restart==0){
 write( parameters.nx(), parameters.ny(), parameters.nz(), local_size, 1, num_procs, 1, num_procs, myid, parameters.WriteTimeCounter, C_bar, U, V, P, Phi); 
}

while(parameters.time<parameters.Tmax){

parameters.time += parameters.dt;

for(int RK4_counter=0; RK4_counter<4; RK4_counter++)
{
  ////C_star calculation
  Calculate_C_star(parameters, communicate, local_size, RK4_counter, C_bar, k_, C_star);
     
  Calculate_Coefficients_From_Table(sig_lam_star, parameters.x_center, parameters.x_face, sigma_star, lambda, C_star, Cs, C0,nondimensional_h);
 
 ////Solve and calculate k
  Solve(sig_lam_star, communicate, solver, C_star, C0, Cs, nondimensional_h, P, Phi, U, V, dC_bar, k_, parameters.x_center, parameters.x_face, parameters.dt, parameters.b(RK4_counter), myid);

 ////Writing
 if(int(parameters.time/parameters.dt)%parameters.period==0 && RK4_counter==0)
      {
          parameters.WriteTimeCounter++;

          write( parameters.nx(), parameters.ny(), parameters.nz(), local_size, 1, num_procs, 1, num_procs, myid, parameters.WriteTimeCounter, C_bar, U, V, P, Phi);

          write_restart_file(parameters.gx(), parameters.ny(), local_size, 1, num_procs, num_procs, myid, parameters.WriteTimeCounter, parameters.time, C0, P, Phi);
      }

}//end of RK4 loop

    ////Update C_bar
    Update_C_bar(parameters, local_size, C_bar, dC_bar);
      
      if(myid==0)
      {
          cout<<parameters.time<<setw(20)<<"C(15,1)="<<C_bar(15,1)<<endl;
          cout<<endl;
      }

}//End of While loop

return(0);
}







 

  
