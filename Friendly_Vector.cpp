#include<iostream>
#include<cmath>
#include<stdio.h>

#include"Friendly_Vector.h"

using namespace std;

Friendly_Vector::Friendly_Vector(unsigned Dim, unsigned Nx, unsigned Ny, unsigned Nz){
    
    Nx_=Nx;
    Ny_=Ny;
    Nz_=Nz;
    
    vector_dimension_=Dim;
    
    
    switch(vector_dimension_){
       
    //making 1D vector:
        case(1):
             vector_=new double[Nx_];
             break;
   
    
    //making 2D vector:
        case(2):
        
        vector_=new double[Nx_*Ny_];
            break;

    
    //making 3D vector:
        case(3):
            
        vector_=new double[Nx_*Ny_*Nz_];
            break;
    }
}

Friendly_Vector::~Friendly_Vector(){

    delete [] vector_;
}

void Friendly_Vector::allocate_memory(unsigned Dim, unsigned Nx, unsigned Ny, unsigned Nz){
    
    Nx_=Nx;
    Ny_=Ny;
    Nz_=Nz;
    
    vector_dimension_=Dim;
    
    switch(vector_dimension_){
            
            //making 1D vector:
        case(1):
            
            delete [] vector_;
            vector_=new double[Nx_];
            break;
            
            
            //making 2D vector:
        case(2):
            
            delete [] vector_;
            vector_=new double[Nx_*Ny_];
            break;
            
            
            //making 3D vector:
        case(3):
            
            delete [] vector_;
            vector_=new double[Nx_*Ny_*Nz_];
            break;
    }
    return;
    
}


double& Friendly_Vector::operator()(unsigned nx, unsigned ny, unsigned nz){
    
    switch(vector_dimension_){
            
        case(1):
            return(vector_[nx]);
            
        case(2):
            return(vector_[ny*Nx_+nx]);
            
        case(3):
            return(vector_[nz*Nx_*Ny_+ny*Nx_+nx]);
    }
    
}

double Friendly_Vector::operator()(unsigned nx, unsigned ny, unsigned nz) const{
    
    switch(vector_dimension_){
            
        case(1):
            return(vector_[nx]);
            
        case(2):
            return(vector_[ny*Nx_+nx]);
            
        case(3):
            return(vector_[nz*Nx_*Ny_+ny*Nx_+nx]);
    }
    
}
