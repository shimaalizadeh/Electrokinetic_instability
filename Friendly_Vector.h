#ifndef FRIENDLY_VECTOR_H
#define FRIENDLY_VECTOR_H

#include<iostream>
#include<cmath>
#include<stdio.h>

using namespace std;

class Friendly_Vector{
    
public:
 
    Friendly_Vector(unsigned Dim=2, unsigned Nx=0, unsigned Ny=0, unsigned Nz=0); //Constructor
    
    ~Friendly_Vector(); //Destructor
    
    void allocate_memory(unsigned Dim=2, unsigned Nx=10, unsigned Ny=10, unsigned Nz=0); //Use in case you dont allocate number at the beginning

    double& operator() (unsigned nx=0, unsigned ny=0, unsigned nz=0); //operator ()
    
    double operator() (unsigned nx=0, unsigned ny=0, unsigned nz=0) const; //operator () pair
    
    double Nx() const {return Nx_;}
    
    double Ny()  const{return Ny_;}
    
private:
    unsigned Nx_, Ny_, Nz_;
    unsigned vector_dimension_;
    double* vector_;
    
    //prevent copying and assignment 
    Friendly_Vector(const Friendly_Vector &);
    Friendly_Vector& operator=(const Friendly_Vector &);
};

#endif



