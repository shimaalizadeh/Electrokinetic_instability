#include "Hypre_Poisson_Solver.h"

Poisson_Solver::Poisson_Solver(unsigned solver_id_,  unsigned max_iteration_, double conv_tolerance_, unsigned my_id, unsigned y_ilower_, unsigned y_iupper_, unsigned y_local_size, unsigned nx_, unsigned ny_, double dy_, unsigned P_BC_L, unsigned P_BC_R, unsigned Phi_BC_L, unsigned Phi_BC_R, double *x_center, double *x_face){

 ////set some default parameters
 number_of_ghost   =1;
 
 solver_id         =solver_id_;
 max_iteration     =max_iteration_;
 conv_tolerance    =conv_tolerance_;
 myid              =my_id;
 nx                =nx_;
 gx                =nx_ + 2*number_of_ghost;
 ny                =ny_;
 dy                =dy_;
 domain_ilower     =y_ilower_*gx;
 domain_iupper     =(y_iupper_+1)*gx-1;
 ilower            =y_ilower_*gx*2;
 iupper            =(y_iupper_+1)*gx*2-1;
 local_size        =y_local_size;
 P_BC_Type_L       =P_BC_L;
 P_BC_Type_R       =P_BC_R;
 Phi_BC_Type_L     =Phi_BC_L;
 Phi_BC_Type_R     =Phi_BC_R;
 
 ////Memory Allocation
 x_center_  =   new double[gx];
 x_face_    =   new double[gx-1];
 x_values   =   new double[local_size*gx*2];
 rhs_values =   new double[local_size*gx*2];
 rows       =   new int[local_size*gx*2];
 solution   =   new double[local_size*gx*2];
 index      =   new int[local_size*gx*2];
 
 ////Setting grid points
 for(int i=0;i<gx-1;i++)
    {
        x_face_[i]=x_face[i];

        x_center_[i]=x_center[i];
    }

        x_center_[gx-1]=x_center[gx-1];
}

void Poisson_Solver::Hypre_Create(){

     ////creating matrix
     HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
     HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
     HYPRE_IJMatrixInitialize(A);

    ////creating the rhs and the solution vectors
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&b);
    HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(b);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&x);
    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(x);
return;
}

void Poisson_Solver::lhs_multiplication(double diag11, double diag12, double diag21, double diag22, double **lhs_1, double **lhs_2){

    //int tt;
    double det;
    det=diag11 * diag22- diag12 * diag21;
   
    lhs_2[0][0]=1./det *(diag22 * lhs_1[0][0] - diag12 * lhs_1[1][0]);

    lhs_2[0][1]=1./det *(diag22 * lhs_1[0][1] - diag12 * lhs_1[1][1]);

    lhs_2[1][0]=1./det *(-diag21 * lhs_1[0][0] + diag11 * lhs_1[1][0]);

    lhs_2[1][1]=1./det *(-diag21 * lhs_1[0][1] + diag11 * lhs_1[1][1]);

    //cout<<diag11<<setw(20)<<diag12<<setw(20)<<diag21<<setw(20)<<diag22<<endl;
    //cout<<det<<setw(20)<<lhs_2[0][0]<<setw(20)<<lhs_2[0][1]<<setw(20)<<lhs_2[1][0]<<setw(20)<<lhs_2[1][1]<<endl;
//cin>>tt;

return;
}

 
void Poisson_Solver::rhs_multiplication(double diag11, double diag12, double diag21, double diag22, double *rhs_1, double *rhs_2){
   
    double det;
    det=diag11 * diag22 - diag12 * diag21;

    rhs_2[0]=1./det *(diag22 * rhs_1[0] -diag12 * rhs_1[1]);

    rhs_2[1]=1./det *(-diag21 * rhs_1[0] + diag11 * rhs_1[1]);

return;
}

void  Poisson_Solver::Hypre_Set_Value(Friendly_Vector & A1_x, Friendly_Vector & B1_x, Friendly_Vector & A2_x, Friendly_Vector & B2_x, Friendly_Vector & A1_y, Friendly_Vector & B1_y, Friendly_Vector & A2_y, Friendly_Vector & B2_y, Friendly_Vector & RHS_1, Friendly_Vector & RHS_2, Friendly_Vector & P, Friendly_Vector & Phi ){

    ////Setting all vectors to zeros
    for(int i=0; i<10; i++)
    {
        values_1[i]=0.;
        values_2[i]=0.;
        cols[i]=0.;
    }


    //I use these variables to store the diagonal blocks
    double diag11, diag12, diag21, diag22;
    double  **lhs_1=new double* [2];
    double  **lhs_2=new double* [2];
for(int i=0;i<2;i++)
{
    lhs_1[i]=new double [2];
    lhs_2[i]=new double [2];
}
    double  *rhs_1=new double [2];
    double  *rhs_2= new double [2];


   ////Setting Boundary Condition in x////
   
   //Note: in one case I have assumed Neumann for P and Drichlet for Phi
   if(P_BC_Type_L==1 &&  P_BC_Type_R==1 && Phi_BC_Type_L==0 && Phi_BC_Type_R==0)
   {      
      for(int j=0; j<local_size;j++)
      {
          ////Left Boundary
          
          //storing diagonal elements
          diag11=-A1_x(0,j);
          diag12=-B1_x(0,j);
          diag21=0.;
          diag22=1.;
          
          nnz=0;
          
          //Diagonal Block: P(0,j+1) Phi(0,j+1)
          lhs_1[0][0]=-A1_x(0,j);
          lhs_1[0][1]=-B1_x(0,j);
          lhs_1[1][0]=0.;
          lhs_1[1][1]=1.;
          lhs_multiplication(diag11, diag12, diag21, diag22, lhs_1, lhs_2);
         
          cols[nnz]=ilower+j*gx*2;
          values_1[nnz]=1.;
          values_2[nnz]=0.;
          
          nnz++;
 
          cols[nnz]=ilower+j*gx*2+1;
          values_1[nnz]=0.;
          values_2[nnz]=1.;
          
          nnz++;
          
          //Right Block: P(1,j+1) Phi(1,j+1)
          lhs_1[0][0]=A1_x(0,j);
          lhs_1[0][1]=B1_x(0,j);
          lhs_1[1][0]=0.;
          lhs_1[1][1]=1.;
          lhs_multiplication(diag11, diag12, diag21, diag22, lhs_1, lhs_2);

          cols[nnz]=ilower+(j*gx+1)*2;
          values_1[nnz]=lhs_2[0][0];
          values_2[nnz]=lhs_2[1][0];
          
          nnz++;

          cols[nnz]=ilower+(j*gx+1)*2+1;
          values_1[nnz]=lhs_2[0][1];
          values_2[nnz]=lhs_2[1][1];
         
          nnz++;

          //Mass Eqn.
          counter=ilower+j*gx*2;
          HYPRE_IJMatrixSetValues(A, 1, &nnz, &counter, cols, values_1);
          x_values[j*gx*2]=P(0,j+1);
          rows[j*gx*2]=counter;

          //Potential Drichlet BC. Eqn.
          counter=ilower+j*gx*2+1;
          HYPRE_IJMatrixSetValues(A, 1, &nnz, &counter, cols, values_2);
          x_values[j*gx*2+1]=Phi(0,j+1);
          rows[j*gx*2+1]=counter;

          //rhs:
          rhs_1[0]= RHS_1(0,j)*(x_center_[1]-x_center_[0]);
          rhs_1[1]= RHS_2(0,j);
          rhs_multiplication(diag11, diag12, diag21, diag22, rhs_1, rhs_2);
          rhs_values[j*gx*2]=rhs_2[0];
          rhs_values[j*gx*2+1]=rhs_2[1];
          /////////////////////////////

          ////Right Boundary
          
          //storing diagonal elements
          diag11=A1_x(nx,j);
          diag12=B1_x(nx,j);
          diag21=0.;
          diag22=1.;
         
          nnz=0;

          //Left Block: P(nx,j+1) Phi(nx,j+1)
          lhs_1[0][0]=-A1_x(nx,j);
          lhs_1[0][1]=-B1_x(nx,j);
          lhs_1[1][0]=0.;
          lhs_1[1][1]=1.;
          lhs_multiplication(diag11, diag12, diag21, diag22, lhs_1, lhs_2);

          cols[nnz]=ilower+((j+1)*gx-2)*2;
          values_1[nnz]=lhs_2[0][0];
          values_2[nnz]=lhs_2[1][0];
         
          nnz++;

          cols[nnz]=ilower+((j+1)*gx-2)*2+1;
          values_1[nnz]=lhs_2[0][1];
          values_2[nnz]=lhs_2[1][1];
         
          nnz++;
          
          //Diagonal Block: P(nx+1,j+1) Phi(nx+1,j+1)
          lhs_1[0][0]=A1_x(nx,j);
          lhs_1[0][1]=B1_x(nx,j);
          lhs_1[1][0]=0.;
          lhs_1[1][1]=1.;
          lhs_multiplication(diag11, diag12, diag21, diag22, lhs_1, lhs_2);

          cols[nnz]=ilower+((j+1)*gx-1)*2;
          values_1[nnz]=1.;
          values_2[nnz]=0.;
         
          nnz++;

          cols[nnz]=ilower+((j+1)*gx-1)*2+1;
          values_1[nnz]=0.;
          values_2[nnz]=1.;
         
          nnz++;

          //Mass Eqn.
          counter=ilower+((j+1)*gx-1)*2;
          HYPRE_IJMatrixSetValues(A, 1, &nnz, &counter, cols, values_1);
          x_values[((j+1)*gx-1)*2]=P(nx+1,j+1);
          rows[((j+1)*gx-1)*2]=counter;

          //Potential Drichlet BC. Eqn.
          counter=ilower+((j+1)*gx-1)*2+1;
          HYPRE_IJMatrixSetValues(A, 1, &nnz, &counter, cols, values_2);
          x_values[((j+1)*gx-1)*2+1]=Phi(nx+1,j+1);
          rows[((j+1)*gx-1)*2+1]=counter;

          //rhs:
         rhs_1[0]=RHS_1(nx+1,j)*(x_center_[nx+1]-x_center_[nx]);
         rhs_1[1]=RHS_2(nx+1,j);
         rhs_multiplication(diag11, diag12, diag21, diag22, rhs_1, rhs_2);
         rhs_values[((j+1)*gx-1)*2]=rhs_2[0];
         rhs_values[((j+1)*gx-1)*2+1]=rhs_2[1];
         ///////////////////

      }//End of for(j...)
   }//End of if

   //Note: For now the second option is Drichlet for both P and Phi
   else
   {
       //Note: for this case because all blocs are Identity block, everything is already set
       // and no need for lhs and rhs multiplication
       for(int j=0; j<local_size;j++)
      { 

      ////Left Boundary
      //Pressure Drichlet BC.
      nnz=0;
          
      cols[nnz]=ilower+j*gx*2;
      values_1[nnz]=1.;
      nnz++;

      cols[nnz]=ilower+(j*gx+1)*2;
      values_1[nnz]=1.;
      nnz++;

      counter=ilower+j*gx*2;
      HYPRE_IJMatrixSetValues(A, 1, &nnz, &counter, cols, values_1);

      rhs_values[j*gx*2]=RHS_1(0,j);
      x_values[j*gx*2]=P(0,j+1);
      rows[j*gx*2]=counter;

      //Potential Drichlet BC.
      nnz=0;

      cols[nnz]=ilower+j*gx*2+1;
      values_2[nnz]=1.;
      nnz++;

      cols[nnz]=ilower+(j*gx+1)*2+1;
      values_2[nnz]=1.;
      nnz++;

      counter=ilower+j*gx*2+1;
      HYPRE_IJMatrixSetValues(A, 1, &nnz, &counter, cols, values_2);

      rhs_values[j*gx*2+1]=RHS_2(0,j);
      x_values[j*gx*2+1]=Phi(0,j+1);
      rows[j*gx*2+1]=counter;
      ///////////////////

      ////Right Boundary
      //Pressure Drichlet BC.
      nnz=0;

      cols[nnz]=ilower+((j+1)*gx-2)*2;
      values_1[nnz]=1.;
      nnz++;

      cols[nnz]=ilower+((j+1)*gx-1)*2;
      values_1[nnz]=1.;
      nnz++;

      counter=ilower+((j+1)*gx-1)*2;
      HYPRE_IJMatrixSetValues(A, 1, &nnz, &counter, cols, values_1);

      rhs_values[((j+1)*gx-1)*2]=RHS_1(nx+1,j);
      x_values[((j+1)*gx-1)*2]=P(nx+1,j+1);
      rows[((j+1)*gx-1)*2]=counter;

      //Potential Drichlet BC.
      nnz=0;

      cols[nnz]=ilower+((j+1)*gx-2)*2+1;
      values_2[nnz]=1.;
      nnz++;

      cols[nnz]=ilower+((j+1)*gx-1)*2+1;
      values_2[nnz]=1.;
      nnz++;

      counter=ilower+((j+1)*gx-1)*2+1;
      HYPRE_IJMatrixSetValues(A, 1, &nnz, &counter, cols, values_2);
      
      rhs_values[((j+1)*gx-1)*2+1]=RHS_2(nx+1,j);
      x_values[((j+1)*gx-1)*2+1]=Phi(nx+1,j+1);
      rows[((j+1)*gx-1)*2+1]=counter;

      ////////////////

     }//End of for(j...)
  }//End of else
 
   /////Now for the rest of the cells////
    
         for(int i=0; i<10; i++)
         {
             values_1[i]=0.;
             values_2[i]=0.;
             cols[i]=0.;
         }

   for(int j=0; j<local_size; j++)
   {
      for(int i=1;i<nx+1;i++)
      {
          //storing diagonal elements
          diag11=-1./(x_face_[i]-x_face_[i-1])*(A1_x(i-1,j)/(x_center_[i]-x_center_[i-1]) + A1_x(i,j)/(x_center_[i+1]-x_center_[i])) -(A1_y(i-1,j)+A1_y(i-1,j+1))/(dy*dy);
          diag11=diag11*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
          
          diag21=-1./(x_face_[i]-x_face_[i-1])*(A2_x(i-1,j)/(x_center_[i]-x_center_[i-1]) + A2_x(i,j)/(x_center_[i+1]-x_center_[i])) -(A2_y(i-1,j)+A2_y(i-1,j+1))/(dy*dy);
          diag21=diag21*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
          
          diag12=-1./(x_face_[i]-x_face_[i-1])*(B1_x(i-1,j)/(x_center_[i]-x_center_[i-1]) + B1_x(i,j)/(x_center_[i+1]-x_center_[i])) -(B1_y(i-1,j)+B1_y(i-1,j+1))/(dy*dy);
          diag12=diag12*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
          
          diag22=-1./(x_face_[i]-x_face_[i-1])*(B2_x(i-1,j)/(x_center_[i]-x_center_[i-1]) + B2_x(i,j)/(x_center_[i+1]-x_center_[i])) -(B2_y(i-1,j)+B2_y(i-1,j+1))/(dy*dy);
          diag22=diag22*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
          
          nnz=0;
          
          ///bottom pressure and potential coefficient (P_(i,j-1) & Phi_(i,j-1))///
          
          lhs_1[0][0]=A1_y(i-1,j)/(dy*dy)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
          lhs_1[0][1]=B1_y(i-1,j)/(dy*dy)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
          lhs_1[1][0]=A2_y(i-1,j)/(dy*dy)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
          lhs_1[1][1]=B2_y(i-1,j)/(dy*dy)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
          lhs_multiplication(diag11, diag12, diag21, diag22, lhs_1, lhs_2);
          

          if(domain_ilower +j*gx+i-gx>=0)
          {
            ////Pressure
            cols[nnz]=ilower+(j*gx+i-gx)*2;

              values_1[nnz]=lhs_2[0][0];
              values_2[nnz]=lhs_2[1][0];

              nnz++;

            ///Potential
            cols[nnz]=ilower+(j*gx+i-gx)*2+1;
 
              values_1[nnz]=lhs_2[0][1];
              values_2[nnz]=lhs_2[1][1];

              nnz++;
          }
          else //periodic boundary condition
          {
              ///Pressure
              cols[nnz]=((ny-1)*gx+i)*2;

              values_1[nnz]=lhs_2[0][0];
              values_2[nnz]=lhs_2[1][0];
              
              nnz++;
              
              ///Potential
              cols[nnz]=((ny-1)*gx+i)*2+1;

              values_1[nnz]=lhs_2[0][1];
              values_2[nnz]=lhs_2[1][1];
              
              nnz++;
          }
         
         ///Top pressure and potential coefficient (P_(i,j+1) & Phi_(i,j+1))///
          
          lhs_1[0][0]=A1_y(i-1,j+1)/(dy*dy)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
          lhs_1[0][1]=B1_y(i-1,j+1)/(dy*dy)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
          lhs_1[1][0]=A2_y(i-1,j+1)/(dy*dy)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
          lhs_1[1][1]=B2_y(i-1,j+1)/(dy*dy)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
          lhs_multiplication(diag11, diag12, diag21, diag22, lhs_1, lhs_2);
          
         if(domain_ilower+j*gx+i+gx<ny*gx)
         {
             ///Pressure
             cols[nnz]=ilower+(j*gx+i+gx)*2;

             values_1[nnz]=lhs_2[0][0];
             values_2[nnz]=lhs_2[1][0];
           
             nnz++;

             ///Potential
             cols[nnz]=ilower+(j*gx+i+gx)*2+1;

             values_1[nnz]=lhs_2[0][1];
             values_2[nnz]=lhs_2[1][1];

             nnz++;
         }
         else //periodic boundary condition
         {
             ///Pressure
             cols[nnz]=i*2;
             
             values_1[nnz]=lhs_2[0][0];
             values_2[nnz]=lhs_2[1][0];
           
             nnz++;

             ///Potential
             cols[nnz]=i*2+1;
             
             values_1[nnz]=lhs_2[0][1];
             values_2[nnz]=lhs_2[1][1];

             nnz++;
         }
        
         ///the left pressure and potential coefficients (P_(i-1,j) & Phi_(i-1,j))///
          
          lhs_1[0][0]=A1_x(i-1,j);
          lhs_1[0][1]=B1_x(i-1,j);
          lhs_1[1][0]=A2_x(i-1,j);
          lhs_1[1][1]=B2_x(i-1,j);
          lhs_multiplication(diag11, diag12, diag21, diag22, lhs_1, lhs_2);
          
          ///Pressure
          cols[nnz]=ilower+(j*gx+i-1)*2;

          values_1[nnz]=lhs_2[0][0];
          values_2[nnz]=lhs_2[1][0];
          
          nnz++;

          ///Potential
          cols[nnz]=ilower+(j*gx+i-1)*2+1;

          values_1[nnz]=lhs_2[0][1];
          values_2[nnz]=lhs_2[1][1];
          
          nnz++;
 
          ///the right pressure and potential coefficients (P_(i+1,j) & Phi_(i+1,j))///
          
          lhs_1[0][0]=A1_x(i,j)/(x_center_[i+1]-x_center_[i])*(x_center_[i]-x_center_[i-1]);
          lhs_1[0][1]=B1_x(i,j)/(x_center_[i+1]-x_center_[i])*(x_center_[i]-x_center_[i-1]);
          lhs_1[1][0]=A2_x(i,j)/(x_center_[i+1]-x_center_[i])*(x_center_[i]-x_center_[i-1]);
          lhs_1[1][1]=B2_x(i,j)/(x_center_[i+1]-x_center_[i])*(x_center_[i]-x_center_[i-1]);
          lhs_multiplication(diag11, diag12, diag21, diag22, lhs_1, lhs_2);
        
          ///Pressure
          cols[nnz]=ilower+(j*gx+i+1)*2;
          
          values_1[nnz]=lhs_2[0][0];
          values_2[nnz]=lhs_2[1][0];

          nnz++;

       
          ///Potential
          cols[nnz]=ilower+(j*gx+i+1)*2+1;
          
          values_1[nnz]=lhs_2[0][1];
          values_2[nnz]=lhs_2[1][1];

          nnz++;

          ///Diagonal Elements///
          lhs_1[0][0]=-1./(x_face_[i]-x_face_[i-1])*(A1_x(i-1,j)/(x_center_[i]-x_center_[i-1]) + A1_x(i,j)/(x_center_[i+1]-x_center_[i])) -(A1_y(i-1,j)+A1_y(i-1,j+1))/(dy*dy);
          lhs_1[0][0]=lhs_1[0][0]*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
          
          lhs_1[1][0]=-1./(x_face_[i]-x_face_[i-1])*(A2_x(i-1,j)/(x_center_[i]-x_center_[i-1]) + A2_x(i,j)/(x_center_[i+1]-x_center_[i])) -(A2_y(i-1,j)+A2_y(i-1,j+1))/(dy*dy);
          lhs_1[1][0]=lhs_1[1][0]*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
          
          lhs_1[0][1]=-1./(x_face_[i]-x_face_[i-1])*(B1_x(i-1,j)/(x_center_[i]-x_center_[i-1]) + B1_x(i,j)/(x_center_[i+1]-x_center_[i])) -(B1_y(i-1,j)+B1_y(i-1,j+1))/(dy*dy);
          lhs_1[0][1]=lhs_1[0][1]*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
          
          lhs_1[1][1]=-1./(x_face_[i]-x_face_[i-1])*(B2_x(i-1,j)/(x_center_[i]-x_center_[i-1]) + B2_x(i,j)/(x_center_[i+1]-x_center_[i])) -(B2_y(i-1,j)+B2_y(i-1,j+1))/(dy*dy);
          lhs_1[1][1]=lhs_1[1][1]*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
          
          lhs_multiplication(diag11, diag12, diag21, diag22, lhs_1, lhs_2);

          ///Pressure
          cols[nnz]=ilower+(j*gx+i)*2;

          values_1[nnz]=lhs_2[0][0];;
          values_2[nnz]=lhs_2[1][0];;
     
          nnz++;

          ///Potential
          cols[nnz]=ilower+(j*gx+i)*2+1;

          values_1[nnz]=lhs_2[0][1];;
          values_2[nnz]=lhs_2[1][1];;

          nnz++;

       ////Setting Values////
          
          //Mass Eqn.
          counter=ilower+(j*gx+i)*2;
          HYPRE_IJMatrixSetValues(A, 1, &nnz, &counter, cols, values_1);
    
          x_values[(j*gx+i)*2]=P(i,j+1);
          rows[(j*gx+i)*2]=counter;

      
          //Current Eqn.
          counter=ilower+(j*gx+i)*2+1;
          HYPRE_IJMatrixSetValues(A, 1, &nnz, &counter, cols, values_2);

          x_values[(j*gx+i)*2+1]=Phi(i,j+1);
          rows[(j*gx+i)*2+1]=counter;

          //rhs:
          rhs_1[0]=RHS_1(i,j)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
          rhs_1[1]=RHS_2(i,j)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
          rhs_multiplication(diag11, diag12, diag21, diag22, rhs_1, rhs_2);
          rhs_values[(j*gx+i)*2]=rhs_2[0];
          rhs_values[(j*gx+i)*2+1]=rhs_2[1];
       
      }
   }

    ////Matrix assembly after setting the coefficients
    HYPRE_IJMatrixAssemble(A);

    ////Get the parcsr matrix object to use
    HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);
    //////////////////////////////////////////////////

    ////Setting Values in rhs and solution vectors
    HYPRE_IJVectorSetValues(b, local_size*gx*2, rows, rhs_values);
    HYPRE_IJVectorSetValues(x, local_size*gx*2, rows, x_values);

    ////Vector assembly
    HYPRE_IJVectorAssemble(b);

    HYPRE_IJVectorGetObject(b, (void **) &par_b);
    HYPRE_IJVectorGetObject(x, (void **) &par_x);

    
    ////deleting defined arrays
    for(int i=0;i<2;i++)
    {
        delete [] lhs_1[i];
        delete [] lhs_2[i];
    }
    
    delete [] lhs_1;
    delete [] lhs_2;
    
    delete [] rhs_1;
    delete [] rhs_2;

return;
}

void  Poisson_Solver::Hypre_Export(Friendly_Vector & A1_x, Friendly_Vector & B1_x, Friendly_Vector & A2_x, Friendly_Vector & B2_x, Friendly_Vector & A1_y, Friendly_Vector & B1_y, Friendly_Vector & A2_y, Friendly_Vector & B2_y, Friendly_Vector & RHS_1, Friendly_Vector & RHS_2){
    
    //I use these variables to store the diagonal blocks
    double diag11, diag12, diag21, diag22;
    double  **lhs_1=new double* [2];
    double  **lhs_2=new double* [2];
    for(int i=0;i<2;i++)
    {
        lhs_1[i]=new double [2];
        lhs_2[i]=new double [2];
    }
    double  *rhs_1=new double [2];
    double  *rhs_2= new double [2];


ostringstream filename_writer;
filename_writer<<"M_"<<myid<<".txt";
string filename=filename_writer.str();

ofstream hypre_export;
hypre_export.open(filename.c_str());

//number of rows
hypre_export<<gx*local_size*2<<endl;
//number of nonzeros
hypre_export<<nx*local_size*2*10+local_size*2*4+local_size*2*4<<endl;
//global lower index
hypre_export<<ilower<<endl;
//global upper index
hypre_export<<iupper<<endl;


    for(int i=0; i<10; i++)
    {
        values_1[i]=0.;
        values_2[i]=0.;
        cols[i]=0.;
    }


//Boundary element:
//this function just export for the first type of boundary condition
    
for(int j=0; j<local_size;j++)
{
    ////////Left Boundary:
    
    //storing diagonal elements
    diag11=-A1_x(0,j);
    diag12=-B1_x(0,j);
    diag21=0.;
    diag22=1.;
    
    nnz=0;
    
    //Diagonal Block: P(0,j+1) Phi(0,j+1)
    lhs_1[0][0]=-A1_x(0,j);
    lhs_1[0][1]=-B1_x(0,j);
    lhs_1[1][0]=0.;
    lhs_1[1][1]=1.;
    lhs_multiplication(diag11, diag12, diag21, diag22, lhs_1, lhs_2);
    
    cols[nnz]=ilower+j*gx*2;
    values_1[nnz]=lhs_2[0][0];
    values_2[nnz]=lhs_2[1][0];
    
    nnz++;
    
    cols[nnz]=ilower+j*gx*2+1;
    values_1[nnz]=lhs_2[0][1];
    values_2[nnz]=lhs_2[1][1];
    
    nnz++;
    
    //Right Block: P(1,j+1) Phi(1,j+1)
    lhs_1[0][0]=A1_x(0,j);
    lhs_1[0][1]=B1_x(0,j);
    lhs_1[1][0]=0.;
    lhs_1[1][1]=1.;
    lhs_multiplication(diag11, diag12, diag21, diag22, lhs_1, lhs_2);
    
    cols[nnz]=ilower+(j*gx+1)*2;
    values_1[nnz]=lhs_2[0][0];
    values_2[nnz]=lhs_2[1][0];
    
    nnz++;
    
    cols[nnz]=ilower+(j*gx+1)*2+1;
    values_1[nnz]=lhs_2[0][1];
    values_2[nnz]=lhs_2[1][1];
    
    nnz++;
    
    //rhs:
    rhs_1[0]= RHS_1(0,j)*(x_center_[1]-x_center_[0]);
    rhs_1[1]= RHS_2(0,j);
    rhs_multiplication(diag11, diag12, diag21, diag22, rhs_1, rhs_2);
    rhs_values[j*gx*2]=rhs_2[0];
    rhs_values[j*gx*2+1]=rhs_2[1];
    
    //Mass Eqn.
    counter=ilower+j*gx*2;
    
    hypre_export<<counter<<setw(20)<<nnz<<setw(20);
    for(int ii=0;ii<nnz;ii++)
    {
        hypre_export<<cols[ii]<<setw(20)<<values_1[ii]<<setw(20);
    }
    hypre_export<<setw(30)<<rhs_values[j*gx*2]<<endl;
    
    //Potential Drichlet BC. Eqn.
    counter=ilower+j*gx*2+1;

    hypre_export<<counter<<setw(20)<<nnz<<setw(20);
    for(int ii=0;ii<nnz;ii++)
    {
        hypre_export<<cols[ii]<<setw(20)<<values_2[ii]<<setw(20);
    }
    hypre_export<<setw(30)<<rhs_values[j*gx*2+1]<<endl;
    /////////////////////////////
    
    ////Right Boundary
    
    //storing diagonal elements
    diag11=A1_x(nx,j);
    diag12=B1_x(nx,j);
    diag21=0.;
    diag22=1.;
    
    nnz=0;
    
    //Left Block: P(nx,j+1) Phi(nx,j+1)
    lhs_1[0][0]=-A1_x(nx,j);
    lhs_1[0][1]=-B1_x(nx,j);
    lhs_1[1][0]=0.;
    lhs_1[1][1]=1.;
    lhs_multiplication(diag11, diag12, diag21, diag22, lhs_1, lhs_2);
    
    cols[nnz]=ilower+((j+1)*gx-2)*2;
    values_1[nnz]=lhs_2[0][0];
    values_2[nnz]=lhs_2[1][0];
    
    nnz++;
    
    cols[nnz]=ilower+((j+1)*gx-2)*2+1;
    values_1[nnz]=lhs_2[0][1];
    values_2[nnz]=lhs_2[1][1];
    
    nnz++;
    
    //Diagonal Block: P(nx+1,j+1) Phi(nx+1,j+1)
    lhs_1[0][0]=A1_x(nx,j);
    lhs_1[0][1]=B1_x(nx,j);
    lhs_1[1][0]=0.;
    lhs_1[1][1]=1.;
    lhs_multiplication(diag11, diag12, diag21, diag22, lhs_1, lhs_2);
    
    cols[nnz]=ilower+((j+1)*gx-1)*2;
    values_1[nnz]=lhs_2[0][0];
    values_2[nnz]=lhs_2[1][0];
    
    nnz++;
    
    cols[nnz]=ilower+((j+1)*gx-1)*2+1;
    values_1[nnz]=lhs_2[0][1];
    values_2[nnz]=lhs_2[1][1];
    
    nnz++;
    
    //rhs:
    rhs_1[0]=RHS_1(nx+1,j)*(x_center_[nx+1]-x_center_[nx]);
    rhs_1[1]=RHS_2(nx+1,j);
    rhs_multiplication(diag11, diag12, diag21, diag22, rhs_1, rhs_2);
    rhs_values[((j+1)*gx-1)*2]=rhs_2[0];
    rhs_values[((j+1)*gx-1)*2+1]=rhs_2[1];
    
        //Mass Eqn.
        counter=ilower+((j+1)*gx-1)*2;
    
        hypre_export<<counter<<setw(20)<<nnz<<setw(20);
        for(int ii=0;ii<nnz;ii++)
        {
            hypre_export<<cols[ii]<<setw(20)<<values_1[ii]<<setw(20);
        }
        hypre_export<<setw(30)<<rhs_values[((j+1)*gx-1)*2]<<endl;
        
        
        //Potential Drichlet BC.
        counter=ilower+((j+1)*gx-1)*2+1;
    
        hypre_export<<counter<<setw(20)<<nnz<<setw(20);
        for(int ii=0;ii<nnz;ii++)
        {
            hypre_export<<cols[ii]<<setw(20)<<values_2[ii]<<setw(20);
        }
        hypre_export<<setw(30)<<rhs_values[((j+1)*gx-1)*2+1]<<endl;
        
}//End of for(j...)
   

    //rest of the cells
    
   for(int i=0; i<10; i++)
    {
        values_1[i]=0.;
        values_2[i]=0.;
        cols[i]=0.;
    }
    
    for(int j=0; j<local_size; j++)
    {
        for(int i=1;i<nx+1;i++)
        {
            //storing diagonal elements
            diag11=-1./(x_face_[i]-x_face_[i-1])*(A1_x(i-1,j)/(x_center_[i]-x_center_[i-1]) + A1_x(i,j)/(x_center_[i+1]-x_center_[i])) -(A1_y(i-1,j)+A1_y(i-1,j+1))/(dy*dy);
            diag11=diag11*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
            
            diag21=-1./(x_face_[i]-x_face_[i-1])*(A2_x(i-1,j)/(x_center_[i]-x_center_[i-1]) + A2_x(i,j)/(x_center_[i+1]-x_center_[i])) -(A2_y(i-1,j)+A2_y(i-1,j+1))/(dy*dy);
            diag21=diag12*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
            
            diag12=-1./(x_face_[i]-x_face_[i-1])*(B1_x(i-1,j)/(x_center_[i]-x_center_[i-1]) + B1_x(i,j)/(x_center_[i+1]-x_center_[i])) -(B1_y(i-1,j)+B1_y(i-1,j+1))/(dy*dy);
            diag12=diag21*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
            
            diag22=-1./(x_face_[i]-x_face_[i-1])*(B2_x(i-1,j)/(x_center_[i]-x_center_[i-1]) + B2_x(i,j)/(x_center_[i+1]-x_center_[i])) -(B2_y(i-1,j)+B2_y(i-1,j+1))/(dy*dy);
            diag22=diag22*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
            
            nnz=0;
            
            ///bottom pressure and potential coefficient (P_(i,j-1) & Phi_(i,j-1))///
            
            lhs_1[0][0]=A1_y(i-1,j)/(dy*dy)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
            lhs_1[0][1]=B1_y(i-1,j)/(dy*dy)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
            lhs_1[1][0]=A2_y(i-1,j)/(dy*dy)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
            lhs_1[1][1]=B2_y(i-1,j)/(dy*dy)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
            lhs_multiplication(diag11, diag12, diag21, diag22, lhs_1, lhs_2);
            
            
            if(domain_ilower +j*gx+i-gx>=0)
            {
                ////Pressure
                cols[nnz]=ilower+(j*gx+i-gx)*2;
                
                values_1[nnz]=lhs_2[0][0];
                values_2[nnz]=lhs_2[1][0];
                
                nnz++;
                
                ///Potential
                cols[nnz]=ilower+(j*gx+i-gx)*2+1;
                
                values_1[nnz]=lhs_2[0][1];
                values_2[nnz]=lhs_2[1][1];
                
                nnz++;
            }
            else //periodic boundary condition
            {
                ///Pressure
                cols[nnz]=((ny-1)*gx+i)*2;
                
                values_1[nnz]=lhs_2[0][0];
                values_2[nnz]=lhs_2[1][0];
                
                nnz++;
                
                ///Potential
                cols[nnz]=((ny-1)*gx+i)*2+1;
                
                values_1[nnz]=lhs_2[0][1];
                values_2[nnz]=lhs_2[1][1];
                
                nnz++;
            }
            
            ///Top pressure and potential coefficient (P_(i,j+1) & Phi_(i,j+1))///
            
            lhs_1[0][0]=A1_y(i-1,j+1)/(dy*dy)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
            lhs_1[0][1]=B1_y(i-1,j+1)/(dy*dy)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
            lhs_1[1][0]=A2_y(i-1,j+1)/(dy*dy)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
            lhs_1[1][1]=B2_y(i-1,j+1)/(dy*dy)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
            lhs_multiplication(diag11, diag12, diag21, diag22, lhs_1, lhs_2);
            
            if(domain_ilower+j*gx+i+gx<ny*gx)
            {
                ///Pressure
                cols[nnz]=ilower+(j*gx+i+gx)*2;
                
                values_1[nnz]=lhs_2[0][0];
                values_2[nnz]=lhs_2[1][0];
                
                nnz++;
                
                ///Potential
                cols[nnz]=ilower+(j*gx+i+gx)*2+1;
                
                values_1[nnz]=lhs_2[0][1];
                values_2[nnz]=lhs_2[1][1];
                
                nnz++;
            }
            else //periodic boundary condition
            {
                ///Pressure
                cols[nnz]=i*2;
                
                values_1[nnz]=lhs_2[0][0];
                values_2[nnz]=lhs_2[1][0];
                
                nnz++;
                
                ///Potential
                cols[nnz]=i*2+1;
                
                values_1[nnz]=lhs_2[0][1];
                values_2[nnz]=lhs_2[1][1];
                
                nnz++;
            }
            
            ///the left pressure and potential coefficients (P_(i-1,j) & Phi_(i-1,j))///
            
            lhs_1[0][0]=A1_x(i-1,j);
            lhs_1[0][1]=B1_x(i-1,j);
            lhs_1[1][0]=A2_x(i-1,j);
            lhs_1[1][1]=B2_x(i-1,j);
            lhs_multiplication(diag11, diag12, diag21, diag22, lhs_1, lhs_2);
            
            ///Pressure
            cols[nnz]=ilower+(j*gx+i-1)*2;
            
            values_1[nnz]=lhs_2[0][0];
            values_2[nnz]=lhs_2[1][0];
            
            nnz++;
            
            ///Potential
            cols[nnz]=ilower+(j*gx+i-1)*2+1;
            
            values_1[nnz]=lhs_2[0][1];
            values_2[nnz]=lhs_2[1][1];
            
            nnz++;
            
            ///the right pressure and potential coefficients (P_(i+1,j) & Phi_(i+1,j))///
            
            lhs_1[0][0]=A1_x(i,j)/(x_center_[i+1]-x_center_[i])*(x_center_[i]-x_center_[i-1]);
            lhs_1[0][1]=B1_x(i,j)/(x_center_[i+1]-x_center_[i])*(x_center_[i]-x_center_[i-1]);
            lhs_1[1][0]=A2_x(i,j)/(x_center_[i+1]-x_center_[i])*(x_center_[i]-x_center_[i-1]);
            lhs_1[1][1]=B2_x(i,j)/(x_center_[i+1]-x_center_[i])*(x_center_[i]-x_center_[i-1]);
            lhs_multiplication(diag11, diag12, diag21, diag22, lhs_1, lhs_2);
            
            ///Pressure
            cols[nnz]=ilower+(j*gx+i+1)*2;
            
            values_1[nnz]=lhs_2[0][0];
            values_2[nnz]=lhs_2[1][0];
            
            nnz++;
            
            
            ///Potential
            cols[nnz]=ilower+(j*gx+i+1)*2+1;
            
            values_1[nnz]=lhs_2[0][1];
            values_2[nnz]=lhs_2[1][1];
            
            nnz++;
            
            ///Diagonal Elements///
            lhs_1[0][0]=-1./(x_face_[i]-x_face_[i-1])*(A1_x(i-1,j)/(x_center_[i]-x_center_[i-1]) + A1_x(i,j)/(x_center_[i+1]-x_center_[i])) -(A1_y(i-1,j)+A1_y(i-1,j+1))/(dy*dy);
            lhs_1[0][0]=lhs_1[0][0]*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
            
            lhs_1[1][0]=-1./(x_face_[i]-x_face_[i-1])*(A2_x(i-1,j)/(x_center_[i]-x_center_[i-1]) + A2_x(i,j)/(x_center_[i+1]-x_center_[i])) -(A2_y(i-1,j)+A2_y(i-1,j+1))/(dy*dy);
            lhs_1[1][0]=lhs_1[1][0]*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
            
            lhs_1[0][1]=-1./(x_face_[i]-x_face_[i-1])*(B1_x(i-1,j)/(x_center_[i]-x_center_[i-1]) + B1_x(i,j)/(x_center_[i+1]-x_center_[i])) -(B1_y(i-1,j)+B1_y(i-1,j+1))/(dy*dy);
            lhs_1[0][1]=lhs_1[0][1]*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
            
            lhs_1[1][1]=-1./(x_face_[i]-x_face_[i-1])*(B2_x(i-1,j)/(x_center_[i]-x_center_[i-1]) + B2_x(i,j)/(x_center_[i+1]-x_center_[i])) -(B2_y(i-1,j)+B2_y(i-1,j+1))/(dy*dy);
            lhs_1[1][1]=lhs_1[1][1]*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
            
            lhs_multiplication(diag11, diag12, diag21, diag22, lhs_1, lhs_2);
            
            ///Pressure
            cols[nnz]=ilower+(j*gx+i)*2;
            
            values_1[nnz]=lhs_2[0][0];;
            values_2[nnz]=lhs_2[1][0];;
            
            nnz++;
            
            ///Potential
            cols[nnz]=ilower+(j*gx+i)*2+1;
            
            values_1[nnz]=lhs_2[0][1];;
            values_2[nnz]=lhs_2[1][1];;
            
            nnz++;
            
            //rhs:
            rhs_1[0]=RHS_1(i,j)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
            rhs_1[1]=RHS_2(i,j)*(x_face_[i]-x_face_[i-1])*(x_center_[i]-x_center_[i-1]);
            rhs_multiplication(diag11, diag12, diag21, diag22, rhs_1, rhs_2);
            rhs_values[(j*gx+i)*2]=rhs_2[0];
            rhs_values[(j*gx+i)*2+1]=rhs_2[1];

            
            ////Setting Values////
            
            //Mass Eqn.
            counter=ilower+(j*gx+i)*2;
           
            hypre_export<<counter<<setw(20)<<nnz<<setw(20);
            for(int ii=0;ii<nnz;ii++)
            {
                hypre_export<<cols[ii]<<setw(20)<<values_1[ii]<<setw(20);
            }
            hypre_export<<setw(30)<<rhs_values[(j*gx+i)*2]<<endl;
            
            //Current Eqn.
            counter=ilower+(j*gx+i)*2+1;
            
            hypre_export<<counter<<setw(20)<<nnz<<setw(20);
            for(int ii=0;ii<nnz;ii++)
            {
                hypre_export<<cols[ii]<<setw(20)<<values_2[ii]<<setw(20);
            }
            hypre_export<<setw(30)<<rhs_values[(j*gx+i)*2+1]<<endl;
        }
    }
   hypre_export.close();
    
    ////deleting defined arrays
    for(int i=0;i<2;i++)
    {
        delete [] lhs_1[i];
        delete [] lhs_2[i];
    }
    
    delete [] lhs_1;
    delete [] lhs_2;
    
    delete [] rhs_1;
    delete [] rhs_2;
    
return;
}

void  Poisson_Solver::Hypre_Solve(){

      ////BiCGSTAB with AMG Preconditioner
    if(solver_id==0)
    {  
      int num_iterations;
      double final_res_norm;

      //Create solver
      HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, &solver);
     
      //HYPRE_ParCSRBiCGSTABSetPrintLevel(solver,3);
      HYPRE_ParCSRBiCGSTABSetMaxIter(solver, max_iteration); 
      HYPRE_ParCSRBiCGSTABSetTol(solver, conv_tolerance);
      HYPRE_ParCSRBiCGSTABSetLogging(solver,1);

      //AMG preconditioner
      HYPRE_BoomerAMGCreate(&precond);
      HYPRE_BoomerAMGSetCoarsenType(precond, 0);
      HYPRE_BoomerAMGSetRelaxType(precond, 3);
      HYPRE_BoomerAMGSetNumSweeps(precond, 1);
      HYPRE_BoomerAMGSetTol(precond, 0.0);
      HYPRE_BoomerAMGSetMaxIter(precond, 3); 
      HYPRE_ParCSRBiCGSTABSetPrecond(solver, (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSolve,
                          (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSetup, precond);


      //Now setup and solve!
      HYPRE_ParCSRBiCGSTABSetup(solver, parcsr_A, par_b, par_x);
      HYPRE_ParCSRBiCGSTABSolve(solver, parcsr_A, par_b, par_x);

      // Run info - needed logging turned on
      HYPRE_ParCSRBiCGSTABGetNumIterations(solver, &num_iterations);
      HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);

      
      if(myid==0)
      {
      //printf("BiCGSTABwAMG\n");
      printf("Iterations = %d" , num_iterations);
      printf("   , Final Relative Residual Norm = %e\n", final_res_norm);
      printf("\n");
      }
      

      // Destroy solver and preconditioner
      HYPRE_ParCSRBiCGSTABDestroy(solver);
      HYPRE_BoomerAMGDestroy(precond);
    }
  
    ////PCG with AMG preconditioner
    else if(solver_id==1)
    {
      int num_iterations;
      double final_res_norm;

      HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

      HYPRE_PCGSetMaxIter(solver, max_iteration); 
      HYPRE_PCGSetTol(solver, conv_tolerance); 
      HYPRE_PCGSetTwoNorm(solver, 1);
      HYPRE_PCGSetPrintLevel(solver, 2); 
      HYPRE_PCGSetLogging(solver, 1); 

      //AMG preconditioner
      HYPRE_BoomerAMGCreate(&precond);
      HYPRE_BoomerAMGSetPrintLevel(precond, 1); 
      HYPRE_BoomerAMGSetCoarsenType(precond, 6);
      HYPRE_BoomerAMGSetRelaxType(precond, 6);
      HYPRE_BoomerAMGSetNumSweeps(precond, 1);
      HYPRE_BoomerAMGSetTol(precond, 0.0);
      HYPRE_BoomerAMGSetMaxIter(precond, 1);
      HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                          (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);

      // Now setup and solve!
      HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
      HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);

      //Run info - needed logging turned on
      HYPRE_PCGGetNumIterations(solver, &num_iterations);
      HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (myid == 0)
      {
         //printf("\n");
         printf("Iterations = %d\n", num_iterations);
         printf("Final Relative Residual Norm = %e\n", final_res_norm);
         printf("\n");
      }

      // Destroy solver and preconditioner
      HYPRE_ParCSRPCGDestroy(solver);
      HYPRE_BoomerAMGDestroy(precond);
   }
 
   ////PCG with Parasails Preconditioner
   else if(solver_id==2)
   {
      int    num_iterations;
      double final_res_norm;
      int      sai_max_levels = 1;
      double   sai_threshold = 0.1;
      double   sai_filter = 0.05;
      int      sai_sym = 1;

      //Create solver
      HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

      HYPRE_PCGSetMaxIter(solver, max_iteration); 
      HYPRE_PCGSetTol(solver, conv_tolerance); 
      HYPRE_PCGSetTwoNorm(solver, 1);
      HYPRE_PCGSetPrintLevel(solver, 2);
      HYPRE_PCGSetLogging(solver, 1);

      //ParaSails preconditioner
      HYPRE_ParaSailsCreate(MPI_COMM_WORLD, &precond);
      HYPRE_ParaSailsSetParams(precond, sai_threshold, sai_max_levels);
      HYPRE_ParaSailsSetFilter(precond, sai_filter);
      HYPRE_ParaSailsSetSym(precond, sai_sym);
      HYPRE_ParaSailsSetLogging(precond, 3);
      HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
                          (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);

      //Now setup and solve!
      HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
      HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);

      // Run info - needed logging turned on
      HYPRE_PCGGetNumIterations(solver, &num_iterations);
      HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

      if (myid == 0)
      {
         printf("\n");
         printf("Iterations = %d\n", num_iterations);
         printf("Final Relative Residual Norm = %e\n", final_res_norm);
         printf("\n");
      }

      //Destory solver and preconditioner
      HYPRE_ParCSRPCGDestroy(solver);
      HYPRE_ParaSailsDestroy(precond);
    }
    
    ////PCG
   else if(solver_id==3)
   {
      int num_iterations;
      double final_res_norm;

      //Create solver
      HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
     
      HYPRE_PCGSetMaxIter(solver, max_iteration);
      HYPRE_PCGSetTol(solver, conv_tolerance);
      HYPRE_PCGSetTwoNorm(solver, 1);
      HYPRE_PCGSetPrintLevel(solver, 2);
      HYPRE_PCGSetLogging(solver, 1); 

      // Now setup and solve!
      HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
      HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);

      // Run info - needed logging turned on
      HYPRE_PCGGetNumIterations(solver, &num_iterations);
      HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
      
      if (myid == 0)
      {
         printf("\n");
         printf("Iterations = %d\n", num_iterations);
         printf("Final Relative Residual Norm = %e\n", final_res_norm);
         printf("\n");
      }
      // Destroy solver
      HYPRE_ParCSRPCGDestroy(solver);
   }

   ////AMG
   else if(solver_id==4)
   {
      int num_iterations;
      double final_res_norm;

      /* Create solver */
      HYPRE_BoomerAMGCreate(&solver);
     
      HYPRE_BoomerAMGSetPrintLevel(solver, 3);  // print solve info + parameters
      HYPRE_BoomerAMGSetCoarsenType(solver, 6); // Falgout coarsening
      HYPRE_BoomerAMGSetRelaxType(solver, 3);   // G-S/Jacobi hybrid relaxation
      HYPRE_BoomerAMGSetNumSweeps(solver, 1);   // Sweeeps on each level
      HYPRE_BoomerAMGSetMaxLevels(solver, 20);  // maximum number of levels
      HYPRE_BoomerAMGSetTol(solver, conv_tolerance); //convergance tolerance     

      // Now setup and solve!
      HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
      HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);

      //Run info - needed logging turned on
      HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
      HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
      if (myid == 0)
      {
         printf("\n");
         printf("Iterations = %d\n", num_iterations);
         printf("Final Relative Residual Norm = %e\n", final_res_norm);
         printf("\n");
      }
      //Destroy solver
      HYPRE_BoomerAMGDestroy(solver);
   }

   ////GMRES with AMG preconditioner 
   else if(solver_id==5)
    {
      int num_iterations;
      double final_res_norm;

      //Create solver
      HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);

      HYPRE_ParCSRGMRESSetMaxIter(solver, max_iteration);
      HYPRE_ParCSRGMRESSetTol(solver, conv_tolerance);
      HYPRE_ParCSRGMRESSetLogging(solver,1);

      //AMG preconditioner
      HYPRE_BoomerAMGCreate(&precond);
      HYPRE_BoomerAMGSetCoarsenType(precond, 0);
      HYPRE_BoomerAMGSetRelaxType(precond, 3);
      HYPRE_BoomerAMGSetNumSweeps(precond, 1);
      HYPRE_BoomerAMGSetTol(precond, 0.0);
      HYPRE_BoomerAMGSetMaxIter(precond, 3);
      HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                          (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);

      //Now setup and solve!
      HYPRE_ParCSRGMRESSetup(solver, parcsr_A, par_b, par_x);
      HYPRE_ParCSRGMRESSolve(solver, parcsr_A, par_b, par_x);

      // Run info - needed logging turned on
      HYPRE_ParCSRGMRESGetNumIterations(solver, &num_iterations);
      HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

      
      if(myid==0)
      {
      //printf("GMRESwAMG\n");
      printf("Iterations = %d" , num_iterations);
      printf("   , Final Relative Residual Norm = %e\n", final_res_norm);
      printf("\n");
      }
      

      // Destroy solver and preconditioner
      HYPRE_ParCSRGMRESDestroy(solver);
      //HYPRE_BoomerAMGDestroy(precond);
    }
   else
   {
      if (myid ==0) printf("Invalid solver id specified.\n");
   }

return;
}

void Poisson_Solver::Hypre_Get_Solution(Friendly_Vector &sol_1, Friendly_Vector &sol_2){

    ////Note: vectors sol_1 and sol_2 have local_size+2 row!

   for (int j=0;j<local_size*gx;j++)
   {
       index[2*j] = ilower + 2*j;

       index[2*j+1] = ilower + 2*j+1;
   }

   HYPRE_IJVectorGetValues(x,local_size*gx*2,index,solution);


   for(int j=0;j<local_size;j++)
   {
      for(int i=0;i<gx;i++)
      {
         sol_1(i,j+1)=solution[(j*gx+i)*2];
    
         sol_2(i,j+1)=solution[(j*gx+i)*2+1];
      }
   }

   if(myid==0) P0=sol_1(0,1);

   //broadcasting sol_1(0,1) to all processors
   MPI_Bcast(&P0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   //subtarcting P0 from all sol_1
   for(int j=0; j<local_size;j++)
   {
      for(int i=0;i<gx;i++)
      {
         sol_1(i,j+1)=sol_1(i,j+1)-P0;
      }
   }

return;
}
 
