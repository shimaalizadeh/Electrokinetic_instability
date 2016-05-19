#ifndef HYPRE_POISSON_SOLVER_H
#define HYPRE_POISSON_SOLVER_H

#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<cmath>
#include<stdio.h>

#include"_hypre_utilities.h"
#include"HYPRE_krylov.h"
#include"HYPRE.h"
#include"HYPRE_parcsr_ls.h"

#include"Friendly_Vector.h"

using namespace std;

class Poisson_Solver{

 public:

   Poisson_Solver(unsigned solver_id_, unsigned max_iteration_, double conv_tolerance_,  unsigned my_id, unsigned y_ilower_, unsigned y_iupper_, unsigned y_local_size, unsigned nx_, unsigned ny_, double dy_, unsigned P_BC_L, unsigned P_BC_R, unsigned Phi_BC_L, unsigned Phi_BC_R, double *x_center, double *x_face);  

   ~Poisson_Solver(){ 

            delete [] x_center_;
            delete [] x_face_;
            delete [] rhs_values;
            delete [] x_values;
            delete [] rows;
            delete [] solution;
            delete [] index;
          
   ////Hypre Clean up
   HYPRE_IJMatrixDestroy(A);
   HYPRE_IJVectorDestroy(b);
   HYPRE_IJVectorDestroy(x);
     }

  void    Hypre_Create();

  void    lhs_multiplication(double diag11, double diag12, double diag21, double diag22, double **lhs_1, double **lhs_2);

  void    rhs_multiplication(double diag11, double diag12, double diag21, double diag22, double *rhs_1, double *rhs_2);

  void    Hypre_Set_Value(Friendly_Vector & A1_x, Friendly_Vector & B1_x, Friendly_Vector & A2_x, Friendly_Vector & B2_x, Friendly_Vector & A1_y, Friendly_Vector & B1_y, Friendly_Vector & A2_y, Friendly_Vector & B2_y,Friendly_Vector & RHS_1, Friendly_Vector & RHS_2, Friendly_Vector &P, Friendly_Vector &Phi);

  void    Hypre_Solve();

  void    Hypre_Get_Solution(Friendly_Vector &sol_1, Friendly_Vector &sol_2 );
 
  void    Hypre_Export(Friendly_Vector & A1_x, Friendly_Vector & B1_x, Friendly_Vector & A2_x, Friendly_Vector & B2_x, Friendly_Vector & A1_y, Friendly_Vector & B1_y, Friendly_Vector & A2_y, Friendly_Vector & B2_y,Friendly_Vector & RHS_1, Friendly_Vector & RHS_2);

 private:

  int                solver_id;        //solver id
  int                myid;
  int                max_iteration;    //max. iteration of solver
  double             conv_tolerance;   //convergance tolerance
  int                nnz;              //number of non-zeros
  int                counter;          //matrix row counter
  int                domain_ilower;    //domain ilower for processor myid
  int                domain_iupper;    //domain iupper for processor myid
  int                ilower;           //matrix ilower for processor myid
  int                iupper;           //matrix iupper for processor myid
  int                local_size;       //y local size for each processor
  int                gx;               //total number of cells in x including the ghost cells
  int                ny;               ///total number of cells in y without ghost cells
  int                nx;               // number of cells in x without ghost cells
  int                number_of_ghost;  //number of ghost cells from each side
  int                P_BC_Type_L;
  int                P_BC_Type_R;      // if 0 --> Drichlet if 1 --> Neumann
  int                Phi_BC_Type_L;
  int                Phi_BC_Type_R;
  double             dy;               //mesh size in y direction
  double             values_1[10];     //nonzero elements for each domain cell from eqn.1(max. # of non-zero=10)
  double             values_2[10];     //nonzero elements for each domain cell from eqn.2(max. # of non-zero=10)
  int                cols[10];         //non-zero element column location in the matrix
  double             *rhs_values;      //rhs vector
  double             *x_values;        //vector of solution set to 0.
  int                *rows;            //row of matrix
  double             *solution;
  int                *index;

  double             P0;              //storing value of P(1,1) for all processors

  //for gathering info in processor 0
  double             *P_Phi;

  //Non-uniform mesh parameters
  double             *x_center_;
  double             *x_face_;

  //Hypre Solver Requirements
  HYPRE_IJMatrix A;
  HYPRE_ParCSRMatrix parcsr_A;
  HYPRE_IJVector b;
  HYPRE_ParVector par_b;
  HYPRE_IJVector x;
  HYPRE_ParVector par_x;
  HYPRE_Solver solver, precond;
};
#endif
