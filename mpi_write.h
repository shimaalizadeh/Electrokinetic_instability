#include<iostream>
#include<fstream>
#include<iomanip>
#include<sstream>
#include<stdio.h>
using namespace std;

#include"mpi.h"
#include"Friendly_Vector.h"

void ToBuf(int nx, int ny, int nz, Friendly_Vector & vector, double *buffer);

void VelocityToBuf(int nx, int ny, int nz, Friendly_Vector & vector, double *buffer);

void Mesh_Generation(int nx_tot, int ny_tot, int nz_tot, int local_ny, int num_proc_x, int num_proc_y, int num_proc_z,int num_procs, int myid, int domain_y_ilower, double *x, double dy, double dz);

void write(int nx_tot, int ny_tot, int nz_tot, int local_ny, int num_proc_x, int num_proc_y, int num_proc_z,int num_procs, int myid, int TimeStepCounter, Friendly_Vector &C_bar, Friendly_Vector &U, Friendly_Vector &V, Friendly_Vector &P, Friendly_Vector &Phi);
////////Declerations/////////

void ToBuf(int nx, int ny, int nz, Friendly_Vector & vector, double *buffer){

////Note:I have already assumed one ghost cell from each side

  int counter=0;

  for(int k=0;k<nz;k++)
  {
     for(int j=1;j<ny+1;j++)
     {
        for(int i=1;i<nx+1;i++)
        {
           buffer[counter]=vector(i,j);
           counter++;
        }
     }
  }
return;
}

void ToBuf_2D(int gx, int ny, Friendly_Vector & vector, double *buffer){

////Note:I have already assumed one ghost cell from each side
////Note: ghost cells in x-directions are also written-->gx=nx+2
////Note: Current simulation is 2D.
  int counter=0;

     for(int j=1;j<ny+1;j++)
     {
        for(int i=0;i<gx;i++)
        {
           buffer[counter]=vector(i,j);
           counter++;
        }
     }

return;
}

void FromBuf_2D(int gx, int ny, Friendly_Vector & vector, double *buffer){

////Note:I have already assumed one ghost cell from each side
////Note: Current simulation is 2D.

  int counter=0;

     for(int j=1;j<ny+1;j++)
     {
        for(int i=0;i<gx;i++)
        {
           vector(i,j)=buffer[counter];
           counter++;
        }
     }

return;
}

void VelocityToBuf(int nx, int ny, int nz, Friendly_Vector & vector, double *buffer){

////Note:I have already assumed one ghost cell from each side

  int counter=0;

  for(int k=0;k<nz;k++)
  {
     for(int j=0;j<ny;j++)
     {
        for(int i=0;i<nx;i++)
        {
           buffer[counter]=vector(i,j);
           counter++;
        }
     }
  }
return;
}

void Mesh_Generation(int nx_tot, int ny_tot, int nz_tot, int local_ny, int num_proc_x, int num_proc_y, int num_proc_z, int num_procs, int myid, int domain_y_ilower, double *x, double dy, double dz){

//Note: I have domain partitioning in y direction ONLY.

MPI_File Myfile_Mesh;
MPI_Datatype filetype;
MPI_Status status;

double *buffer=new double[nx_tot*local_ny*nz_tot];
int domain_num=1;
int counter=0;

int Num_Phyiscal_Cells=nx_tot * local_ny * nz_tot; //for each processor

int Total_Num_Domain_Bytes=nx_tot*ny_tot*nz_tot*sizeof(double); //for the whole domain

int gsizes[3]={nz_tot, ny_tot, nx_tot};

int distribs[3]={MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_BLOCK};

int dargs[3]={MPI_DISTRIBUTE_DFLT_DARG,MPI_DISTRIBUTE_DFLT_DARG,MPI_DISTRIBUTE_DFLT_DARG};

int psizes[3]={num_proc_z, num_proc_y, num_proc_x};

MPI_Type_create_darray(num_procs, myid, 3, gsizes, distribs, dargs, psizes, MPI_ORDER_C, MPI_DOUBLE, &filetype);

MPI_Type_commit(&filetype);

//writing to file
int header[4]={domain_num, nx_tot, ny_tot, nz_tot};

MPI_File_open(MPI_COMM_WORLD, "Mesh.bin", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &Myfile_Mesh);

//writing header by processor 0
if(myid==0) 
{
  MPI_File_write(Myfile_Mesh, header, 4, MPI_INT, &status);
}

int disp=4*sizeof(int);

////X
MPI_File_set_view(Myfile_Mesh, disp+0*Total_Num_Domain_Bytes, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

for(int k=0;k<nz_tot;k++)
{
   for(int j=0;j<local_ny;j++)
   {
     for(int i=1;i<nx_tot+1;i++)
     {
        buffer[counter]=x[i];
        counter++;
     }
   }
}
counter=0;

MPI_File_write_all(Myfile_Mesh, buffer, Num_Phyiscal_Cells, MPI_DOUBLE,&status); 

////Y
MPI_File_set_view(Myfile_Mesh, disp+1*Total_Num_Domain_Bytes, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

for(int k=0;k<nz_tot;k++)
{
   for(int j=0;j<local_ny;j++)
   {
     for(int i=1;i<nx_tot+1;i++)
     {
        buffer[counter]=(domain_y_ilower+j+0.5)*dy;
        counter++;
     }
   }
}
counter=0;

MPI_File_write_all(Myfile_Mesh, buffer, Num_Phyiscal_Cells, MPI_DOUBLE,&status);

////Z
MPI_File_set_view(Myfile_Mesh, disp+2*Total_Num_Domain_Bytes, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

for(int k=0;k<nz_tot;k++)
{
   for(int j=0;j<local_ny;j++)
   {
     for(int i=1;i<nx_tot+1;i++)
     {
        buffer[counter]=(k+0.5)*dz;
        counter++;
     }
   }
}
counter=0;
MPI_File_write_all(Myfile_Mesh, buffer, Num_Phyiscal_Cells, MPI_DOUBLE,&status);

delete [] buffer;
MPI_File_close(&Myfile_Mesh);
MPI_Type_free(&filetype);

return;
}


void write(int nx_tot, int ny_tot, int nz_tot, int local_ny, int num_proc_x, int num_proc_y, int num_proc_z,int num_procs, int myid, int TimeStepCounter, Friendly_Vector &C_bar, Friendly_Vector &U, Friendly_Vector &V, Friendly_Vector &P, Friendly_Vector &Phi){

MPI_File Myfile_Data_Tec;
MPI_Datatype filetype;
MPI_Status status;

int domain_num=1;
int variable_num=6;
double *buffer=new double[nx_tot*local_ny*nz_tot];
int counter=0;

int Num_Phyiscal_Cells=nx_tot * local_ny * nz_tot; //for each processor

int Total_Num_Domain_Bytes=nx_tot*ny_tot*nz_tot*sizeof(double); //for the whole domain

int gsizes[3]={nz_tot, ny_tot, nx_tot};

int distribs[3]={MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_BLOCK};

int dargs[3]={MPI_DISTRIBUTE_DFLT_DARG,MPI_DISTRIBUTE_DFLT_DARG,MPI_DISTRIBUTE_DFLT_DARG};

int psizes[3]={num_proc_z, num_proc_y, num_proc_x};

MPI_Type_create_darray(num_procs, myid, 3, gsizes, distribs, dargs, psizes, MPI_ORDER_C, MPI_DOUBLE, &filetype);

MPI_Type_commit(&filetype);

////creating file name
string filename_base_Tec="Tecplot";
ostringstream filename_out_Data_Tec;

if(TimeStepCounter<10){
filename_out_Data_Tec<<filename_base_Tec<<"_000"<<TimeStepCounter<<".bin";
}
else if(TimeStepCounter>9 && TimeStepCounter<100){
filename_out_Data_Tec<<filename_base_Tec<<"_00"<<TimeStepCounter<<".bin";
}
else{
filename_out_Data_Tec<<filename_base_Tec<<"_0"<<TimeStepCounter<<".bin";
}

string filename_Tec=filename_out_Data_Tec.str();

int header[5]={domain_num, nx_tot, ny_tot, nz_tot, variable_num};

MPI_File_open(MPI_COMM_WORLD, (char*)(filename_Tec.c_str()), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &Myfile_Data_Tec);

//writing header by processor 0
if(myid==0)
{
  MPI_File_write(Myfile_Data_Tec, header, 5, MPI_INT, &status);
}

////writing C_bar into file
int disp=5*sizeof(int);

MPI_File_set_view(Myfile_Data_Tec, disp+0*Total_Num_Domain_Bytes, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

ToBuf(nx_tot, local_ny, nz_tot, C_bar, buffer);

MPI_File_write_all(Myfile_Data_Tec, buffer, Num_Phyiscal_Cells, MPI_DOUBLE,&status);

////writing U into file
MPI_File_set_view(Myfile_Data_Tec, disp+1*Total_Num_Domain_Bytes, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

VelocityToBuf(nx_tot, local_ny, nz_tot, U, buffer);

MPI_File_write_all(Myfile_Data_Tec, buffer, Num_Phyiscal_Cells, MPI_DOUBLE,&status);

////writing V into file
MPI_File_set_view(Myfile_Data_Tec, disp+2*Total_Num_Domain_Bytes, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

VelocityToBuf(nx_tot, local_ny, nz_tot, V, buffer);

MPI_File_write_all(Myfile_Data_Tec, buffer, Num_Phyiscal_Cells, MPI_DOUBLE,&status);

////writing W=0 into file
Friendly_Vector W(2,nx_tot, local_ny);
for(int i=0; i<nx_tot; i++)
{
  for(int j=0;j<local_ny; j++)
  {
    W(i,j)=0.0;
  }
}

MPI_File_set_view(Myfile_Data_Tec, disp+3*Total_Num_Domain_Bytes, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

VelocityToBuf(nx_tot, local_ny, nz_tot, W, buffer);

MPI_File_write_all(Myfile_Data_Tec, buffer, Num_Phyiscal_Cells, MPI_DOUBLE,&status);

////writing Pressure
MPI_File_set_view(Myfile_Data_Tec, disp+4*Total_Num_Domain_Bytes, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

ToBuf(nx_tot, local_ny, nz_tot, P, buffer);

MPI_File_write_all(Myfile_Data_Tec, buffer, Num_Phyiscal_Cells, MPI_DOUBLE,&status);

////write Potential
MPI_File_set_view(Myfile_Data_Tec, disp+5*Total_Num_Domain_Bytes, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

ToBuf(nx_tot, local_ny, nz_tot, Phi, buffer);

MPI_File_write_all(Myfile_Data_Tec, buffer, Num_Phyiscal_Cells, MPI_DOUBLE,&status);

delete [] buffer;
MPI_File_close(&Myfile_Data_Tec);
MPI_Type_free(&filetype);

return;
}

void write_restart_file(int gx_tot, int ny_tot, int local_ny, int num_proc_x, int num_proc_y, int num_procs, int myid, int TimeStepCounter, double time, Friendly_Vector &C0, Friendly_Vector &P, Friendly_Vector &Phi){

//Note: ghost cells in x are also written -->gx=nx+2
MPI_File Myfile_Data_restart;
MPI_Datatype filetype;
MPI_Status status;

double *buffer=new double[gx_tot*local_ny];
int counter=0;

int Num_Phyiscal_Cells=gx_tot * local_ny; //for each processor including ghost

int Total_Num_Domain_Bytes=gx_tot*ny_tot*sizeof(double); //for the whole domain including ghost

int gsizes[2]={ny_tot, gx_tot};

int distribs[2]={MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_BLOCK};

int dargs[2]={MPI_DISTRIBUTE_DFLT_DARG,MPI_DISTRIBUTE_DFLT_DARG};

int psizes[2]={num_proc_y, num_proc_x};

MPI_Type_create_darray(num_procs, myid, 2, gsizes, distribs, dargs, psizes, MPI_ORDER_C, MPI_DOUBLE, &filetype);

MPI_Type_commit(&filetype);

int header_0[2] = {gx_tot,ny_tot};
int header_1[1] = {TimeStepCounter};
double header_2[1] = {time};

MPI_File_open(MPI_COMM_WORLD, "restart.bin", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &Myfile_Data_restart);

//writing header by processor 0
if(myid==0)
{
  MPI_File_write(Myfile_Data_restart, header_0, 2, MPI_INT, &status);
  MPI_File_write(Myfile_Data_restart, header_1, 1, MPI_INT, &status);
  MPI_File_write(Myfile_Data_restart, header_2, 1, MPI_DOUBLE, &status);
}

int disp=3*sizeof(int) + 1*sizeof(double);

////writing C0 into file
MPI_File_set_view(Myfile_Data_restart, disp+0*Total_Num_Domain_Bytes, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

ToBuf_2D(gx_tot, local_ny, C0, buffer);

MPI_File_write_all(Myfile_Data_restart, buffer, Num_Phyiscal_Cells, MPI_DOUBLE,&status);

////writing P into file
MPI_File_set_view(Myfile_Data_restart, disp+1*Total_Num_Domain_Bytes, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

ToBuf_2D(gx_tot, local_ny, P, buffer);

MPI_File_write_all(Myfile_Data_restart, buffer, Num_Phyiscal_Cells, MPI_DOUBLE,&status);

////writing Phi into file
MPI_File_set_view(Myfile_Data_restart, disp+2*Total_Num_Domain_Bytes, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

ToBuf_2D(gx_tot, local_ny, Phi, buffer);

MPI_File_write_all(Myfile_Data_restart, buffer, Num_Phyiscal_Cells, MPI_DOUBLE,&status);

delete [] buffer;
MPI_File_close(&Myfile_Data_restart);
MPI_Type_free(&filetype);

return;
}

void read_restart_file(int gx_tot, int ny_tot, int local_ny, int num_proc_x, int num_proc_y, int num_procs, int myid, int &TimeStepCounter, double &time, Friendly_Vector &C0, Friendly_Vector &P, Friendly_Vector &Phi){

//Note: ghost cells in x are also written -->gx=nx+2

MPI_File Myfile_Data_restart;
MPI_Datatype filetype;
MPI_Status status;

double *buffer=new double[gx_tot*local_ny];
int counter=0;

int Num_Phyiscal_Cells=gx_tot * local_ny; //for each processor

int Total_Num_Domain_Bytes=gx_tot*ny_tot*sizeof(double); //for the whole domain

int gsizes[2]={ny_tot, gx_tot};

int distribs[2]={MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_BLOCK};

int dargs[2]={MPI_DISTRIBUTE_DFLT_DARG,MPI_DISTRIBUTE_DFLT_DARG};

int psizes[2]={num_proc_y, num_proc_x};

MPI_Type_create_darray(num_procs, myid, 2, gsizes, distribs, dargs, psizes, MPI_ORDER_C, MPI_DOUBLE, &filetype);

MPI_Type_commit(&filetype);

MPI_File_open(MPI_COMM_WORLD,"restart.bin",MPI_MODE_RDONLY,MPI_INFO_NULL,&Myfile_Data_restart);

int dum_int[2]={};

if (myid==0)
  {
    MPI_File_read(Myfile_Data_restart, dum_int, 2, MPI_INT,&status);
    MPI_File_read(Myfile_Data_restart, &TimeStepCounter, 1, MPI_INT,&status);
    MPI_File_read(Myfile_Data_restart, &time, 1, MPI_DOUBLE,&status);
  }

MPI_Bcast(&TimeStepCounter,1,MPI_INT,0,MPI_COMM_WORLD);
MPI_Bcast(&time,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

int disp=3*sizeof(int) + 1*sizeof(double);
	    
////reading C0
MPI_File_set_view(Myfile_Data_restart,disp+0*Total_Num_Domain_Bytes,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);
MPI_File_read_all(Myfile_Data_restart, buffer, Num_Phyiscal_Cells, MPI_DOUBLE, &status);
FromBuf_2D(gx_tot, local_ny, C0, buffer);

////reading P
MPI_File_set_view(Myfile_Data_restart,disp+1*Total_Num_Domain_Bytes,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);
MPI_File_read_all(Myfile_Data_restart, buffer, Num_Phyiscal_Cells, MPI_DOUBLE, &status);
FromBuf_2D(gx_tot, local_ny, P, buffer);

////reading Phi
MPI_File_set_view(Myfile_Data_restart,disp+2*Total_Num_Domain_Bytes,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL);
MPI_File_read_all(Myfile_Data_restart, buffer, Num_Phyiscal_Cells, MPI_DOUBLE, &status);
FromBuf_2D(gx_tot, local_ny, Phi, buffer);

delete [] buffer;
MPI_File_close(&Myfile_Data_restart);
MPI_Type_free(&filetype);

return;
}

