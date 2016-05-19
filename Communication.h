#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include<iostream>
#include<math.h>
#include<cmath>
#include<stdio.h>
#include"mpi.h"
#include"Friendly_Vector.h"

using namespace std;

class Communication{

public:

Communication(unsigned num_procs=1, unsigned myid=0, unsigned nx=10, unsigned local_ny=10){
 
              num_procs_=num_procs;
              myid_=myid;
              nx_=nx;
              ny_=local_ny;

              //memory allocation
              sender_bottom_=new double[nx_];

              sender_top_=new double[nx_];

              receiver_bottom_=new double[nx_];

              receiver_top_=new double[nx_];
}

~Communication(){
              
               delete [] sender_bottom_;
               delete [] sender_top_;
               delete [] receiver_bottom_;
               delete [] receiver_top_;
 }

void SendRecv(Friendly_Vector & vector){

           if(num_procs_==1)
           {     
                 for(int i=1; i<nx_+1; i++)
                 { 
                   vector(i,0)=vector(i,ny_);
        
                   vector(i,ny_+1)=vector(i,1);
                 }
           }
           else
           {
                top_=(myid_+1)%num_procs_;
                bottom_=myid_-1;
                if(bottom_<0) bottom_=num_procs_-1;

                for(int i=0;i<nx_;i++)
                {
                   sender_bottom_[i]=vector(i+1,1);
                   sender_top_[i]=vector(i+1,ny_);
                }

                ////Communication
                MPI_Irecv(receiver_bottom_, nx_, MPI_DOUBLE, bottom_, 100, MPI_COMM_WORLD, &request_recv_[0]);
                MPI_Isend(sender_top_, nx_, MPI_DOUBLE, top_, 100, MPI_COMM_WORLD, &request_send_[0]);

                MPI_Irecv(receiver_top_, nx_, MPI_DOUBLE, top_, 200, MPI_COMM_WORLD, &request_recv_[1]);
                MPI_Isend(sender_bottom_, nx_, MPI_DOUBLE, bottom_, 200, MPI_COMM_WORLD, &request_send_[1]);

               ////wait for the communications to be done
               MPI_Waitall(2, &request_send_[0],&status[0]);
               MPI_Waitall(2, &request_recv_[0],&status[0]);

               
               for(int i=0;i<nx_;i++)
               {
                    vector(i+1,0)=receiver_bottom_[i];

                    vector(i+1,ny_+1)=receiver_top_[i];
               }

            }//End of else
          
return;
}


private:

int     num_procs_;
int     myid_;
int     nx_;
int     ny_;
double  *sender_bottom_;
double  *sender_top_;
double  *receiver_bottom_;
double  *receiver_top_; 

int      bottom_;
int      top_;
MPI_Request  request_recv_[2], request_send_[2];
MPI_Status   status[2];

};
#endif
