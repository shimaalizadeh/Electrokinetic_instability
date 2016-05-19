#include "TableReading.h"

using namespace std;

Tablereading::Tablereading(const char* filename,  unsigned gx, unsigned gy, double dy,  double Pe, double gp, unsigned P_BC_L, unsigned P_BC_R, unsigned Phi_BC_L, unsigned Phi_BC_R, double Phi_L, double Phi_R, double P_L, double P_R,  unsigned max_sigma_star_counter, unsigned max_lambda_star_counter, double sigma_star_min, double lambda_star_min){ 
   
    //set some values:
    max_sigma_star_counter_=max_sigma_star_counter;
    
    max_lambda_star_counter_=max_lambda_star_counter;
   
    Pe_=Pe;
 
    gp_=gp;

    gx_=gx;

    gy_=gy;

    dy_=dy;
 
    P_BC_Type_L_=P_BC_L;
   
    P_BC_Type_R_=P_BC_R;
    
    Phi_BC_Type_L_=Phi_BC_L;
 
    Phi_BC_Type_R_=Phi_BC_R;
    
    Phi_L_= Phi_L;

    Phi_R_=Phi_R;
    
    P_L_=P_L;

    P_R_=P_R;
   
    sigma_star_min_=sigma_star_min;
    
    lambda_star_min_=lambda_star_min;
    
    
    //memory allocation:
  
    table_sigma_star=new double[max_sigma_star_counter*max_lambda_star_counter];
    table_lambda_star=new double[max_sigma_star_counter*max_lambda_star_counter];
    table_f_bar=new double[max_sigma_star_counter*max_lambda_star_counter];
    table_Psi_center=new double[max_sigma_star_counter*max_lambda_star_counter];
    table_ge_bar=new double[max_sigma_star_counter*max_lambda_star_counter];
    table_gc_bar=new double[max_sigma_star_counter*max_lambda_star_counter];
    table_gp_m=new double[max_sigma_star_counter*max_lambda_star_counter];
    table_ge_m=new double[max_sigma_star_counter*max_lambda_star_counter];
    table_gc_m=new double[max_sigma_star_counter*max_lambda_star_counter];
    table_gp_p=new double[max_sigma_star_counter*max_lambda_star_counter];
    table_ge_p=new double[max_sigma_star_counter*max_lambda_star_counter];
    table_gc_p=new double[max_sigma_star_counter*max_lambda_star_counter];
 
    //Note:values in cell centers
    f_bar_.allocate_memory(2,gx_, gy_);
    ge_bar_.allocate_memory(2,gx_, gy_);
    gc_bar_.allocate_memory(2,gx_, gy_);
    gp_m_.allocate_memory(2,gx_, gy_);
    ge_m_.allocate_memory(2,gx_, gy_);
    gc_m_.allocate_memory(2,gx_, gy_);
    gp_p_.allocate_memory(2,gx_, gy_);
    ge_p_.allocate_memory(2,gx_, gy_);
    gc_p_.allocate_memory(2,gx_, gy_);
   
    //Note:values in cell centers
    Gp_m_.allocate_memory(2,gx_, gy_);
    Ge_m_.allocate_memory(2,gx_, gy_);
    Gc_m_.allocate_memory(2,gx_, gy_); 
   
    //Note:values in cell centers
    A1_.allocate_memory(2,gx_, gy_);
    A2_.allocate_memory(2,gx_, gy_);
    B1_.allocate_memory(2,gx_, gy_);
    B2_.allocate_memory(2,gx_, gy_);
    f1_.allocate_memory(2,gx_, gy_);
    f2_.allocate_memory(2,gx_, gy_);
     
    RHS_1_.allocate_memory(2,gx_,gy_-2); //Note: does not include 2 rows in y
    RHS_2_.allocate_memory(2,gx_,gy_-2); //Note: does not include 2 rows in y

     //Note:values in faces
     A1_x_.allocate_memory(2,gx_-1, gy_-2);
     A2_x_.allocate_memory(2,gx_-1, gy_-2);
     B1_x_.allocate_memory(2,gx_-1, gy_-2);
     B2_x_.allocate_memory(2,gx_-1, gy_-2);
     f1_x_.allocate_memory(2,gx_-1, gy_-2);
     f2_x_.allocate_memory(2,gx_-1, gy_-2);

     a1_x_.allocate_memory(2,gx_-1, gy_-2);
     b1_x_.allocate_memory(2,gx_-1, gy_-2);
     c1_x_.allocate_memory(2,gx_-1, gy_-2);

     adv_x_.allocate_memory(2,gx_-1, gy_-2);
     diff_x_.allocate_memory(2,gx_-1, gy_-2);
     elec_mig_x_.allocate_memory(2,gx_-1, gy_-2);
     Flux_x_.allocate_memory(2,gx_-1, gy_-2);
     
     A1_y_.allocate_memory(2,gx_-2, gy_-1);
     A2_y_.allocate_memory(2,gx_-2, gy_-1);
     B1_y_.allocate_memory(2,gx_-2, gy_-1);
     B2_y_.allocate_memory(2,gx_-2, gy_-1);
     f1_y_.allocate_memory(2,gx_-2, gy_-1);
     f2_y_.allocate_memory(2,gx_-2, gy_-1);

     a1_y_.allocate_memory(2,gx_-2, gy_-1);
     b1_y_.allocate_memory(2,gx_-2, gy_-1);
     c1_y_.allocate_memory(2,gx_-2, gy_-1);

     adv_y_.allocate_memory(2,gx_-2, gy_-1);
     diff_y_.allocate_memory(2,gx_-2, gy_-1);
     elec_mig_y_.allocate_memory(2,gx_-2, gy_-1);
     Flux_y_.allocate_memory(2,gx_-2, gy_-1);
    
     
    //Importing from input file: 
    ifstream importing(filename);
    
    for(int i=0; i<max_sigma_star_counter*max_lambda_star_counter;i++)
    {
        importing>>table_sigma_star[i]>> table_lambda_star[i]>>table_f_bar[i]>>
        
        table_Psi_center[i]>>table_ge_bar[i]>> table_gc_bar[i]>>table_gp_m[i]>>
        
        table_ge_m[i]>>table_gc_m[i]>>table_gp_p[i]>>table_ge_p[i]>>table_gc_p[i];
        
    }
    
    importing.close();
}

Tablereading::~Tablereading(){
    
    //Destroying the object
    delete [] table_sigma_star;
    delete [] table_lambda_star;
    delete [] table_f_bar;
    delete [] table_Psi_center;
    delete [] table_ge_bar;
    delete [] table_gc_bar;
    delete [] table_gp_m;
    delete [] table_ge_m;
    delete [] table_gc_m;
    delete [] table_gp_p;
    delete [] table_ge_p;
    delete [] table_gc_p;
    
}

void Tablereading::find_from_table(Friendly_Vector & sigma, Friendly_Vector & lambda, Friendly_Vector & C_bar){
   
    double lambda_star;
    
    if(floor(abs(sigma(1,1))/abs(sigma_star_min_))==0)
    {   
       
        for(int j=0;j<gy_;j++)
        {
           for(int i=0;i<gx_;i++)
           {
            f_bar_(i,j)=1.;
            ge_bar_(i,j)=0.;
            gc_bar_(i,j)=0.;
            gp_m_(i,j)=0.;
            ge_m_(i,j)=0.;
            gc_m_(i,j)=0.;
            gp_p_(i,j)=0.;
            ge_p_(i,j)=0.;
            gc_p_(i,j)=0.;
            }
        }
    }
    
    else{

         for(int j=0;j<gy_;j++)
         {
           for(int i=0;i<gx_;i++)
           {         
             if(i%(gx_-1)==0 && j%(gy_-1)==0)
              {
                continue;
              }
             else
              {
        /*      
        lambda_star=lambda(i,j)/sqrt(C_bar(i,j));

        sigma_power_=int(10*log10(sigma(i,j)/sigma_star_min_));
        
        lambda_power_=int(10*log10(lambda_star/lambda_star_min_));
        
        aa=(log10(abs(sigma(i,j)))-log10(abs(table_sigma_star[max_lambda_star_counter_*sigma_power_+lambda_power_])))
        
        /(log10(abs(table_sigma_star[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_]))-
          
          log10(abs(table_sigma_star[max_lambda_star_counter_*sigma_power_+lambda_power_])));
        
        bb=(log10(lambda_star)-log10(table_lambda_star[max_lambda_star_counter_*sigma_power_+lambda_power_]))/
        
        (log10(table_lambda_star[max_lambda_star_counter_*sigma_power_+lambda_power_+1])-
         
         log10(table_lambda_star[max_lambda_star_counter_*sigma_power_+lambda_power_]));
        
        //f_bar
        dum1=aa*(log10(abs(table_f_bar[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_]))-
                 
                 log10(abs(table_f_bar[max_lambda_star_counter_*sigma_power_+lambda_power_])))+
        
        log10(abs(table_f_bar[max_lambda_star_counter_*sigma_power_+lambda_power_]));
        
        dum2=aa*(log10(abs(table_f_bar[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_+1]))
                 
                 -log10(abs(table_f_bar[max_lambda_star_counter_*sigma_power_+lambda_power_+1])))+
        
        log10(abs(table_f_bar[max_lambda_star_counter_*sigma_power_+lambda_power_+1]));
        
        f_bar_(i,j)=bb*(dum2-dum1)+dum1;
        
        f_bar_(i,j)=pow(10,f_bar_(i,j));
        

        //ge_bar
        dum1=aa*(log10(abs(table_ge_bar[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_]))-
                 
                 log10(abs(table_ge_bar[max_lambda_star_counter_*sigma_power_+lambda_power_])))+
        
        log10(abs(table_ge_bar[max_lambda_star_counter_*sigma_power_+lambda_power_]));
        
        dum2=aa*(log10(abs(table_ge_bar[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_+1]))-
                 
                 log10(abs(table_ge_bar[max_lambda_star_counter_*sigma_power_+lambda_power_+1])))+
        
        log10(abs(table_ge_bar[max_lambda_star_counter_*sigma_power_+lambda_power_+1]));
        
        ge_bar_(i,j)=bb*(dum2-dum1)+dum1;
        
        ge_bar_(i,j)=-pow(10,ge_bar_(i,j));
   
       
        //gc_bar
        dum1=aa*(log10(abs(table_gc_bar[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_]))-
                 
                 log10(abs(table_gc_bar[max_lambda_star_counter_*sigma_power_+lambda_power_])))+
        
        log10(abs(table_gc_bar[max_lambda_star_counter_*sigma_power_+lambda_power_]));
        
        
        dum2=aa*(log10(abs(table_gc_bar[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_+1]))-
                 
                 log10(abs(table_gc_bar[max_lambda_star_counter_*sigma_power_+lambda_power_+1])))+
        
        log10(abs(table_gc_bar[max_lambda_star_counter_*sigma_power_+lambda_power_+1]));
        
        gc_bar_(i,j)=bb*(dum2-dum1)+dum1;
        
        gc_bar_(i,j)=-pow(10,gc_bar_(i,j));
        
        
        //gp_m
        dum1=aa*(log10(abs(table_gp_m[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_]))-
                 
                 log10(abs(table_gp_m[max_lambda_star_counter_*sigma_power_+lambda_power_])))+
        
        log10(abs(table_gp_m[max_lambda_star_counter_*sigma_power_+lambda_power_]));
        
        dum2=aa*(log10(abs(table_gp_m[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_+1]))
                 
                 -log10(abs(table_gp_m[max_lambda_star_counter_*sigma_power_+lambda_power_+1])))+
        
        log10(abs(table_gp_m[max_lambda_star_counter_*sigma_power_+lambda_power_+1]));
        
        gp_m_(i,j)=bb*(dum2-dum1)+dum1;
        
        gp_m_(i,j)=pow(10,gp_m_(i,j));
        
        
        //ge_m
        dum1=aa*(log10(abs(table_ge_m[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_]))-
                 
                 log10(abs(table_ge_m[max_lambda_star_counter_*sigma_power_+lambda_power_])))+
        
        log10(abs(table_ge_m[max_lambda_star_counter_*sigma_power_+lambda_power_]));
        
        dum2=aa*(log10(abs(table_ge_m[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_+1]))-
                 
                 log10(abs(table_ge_m[max_lambda_star_counter_*sigma_power_+lambda_power_+1])))+
        
        log10(abs(table_ge_m[max_lambda_star_counter_*sigma_power_+lambda_power_+1]));
        
        ge_m_(i,j)=bb*(dum2-dum1)+dum1;
        
        ge_m_(i,j)=pow(10,ge_m_(i,j));
        
        
        //gc_m
        dum1=aa*(log10(abs(table_gc_m[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_]))-
                 
                 log10(abs(table_gc_m[max_lambda_star_counter_*sigma_power_+lambda_power_])))+
        
        log10(abs(table_gc_m[max_lambda_star_counter_*sigma_power_+lambda_power_]));
        
        dum2=aa*(log10(abs(table_gc_m[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_+1]))-
                 
                 log10(abs(table_gc_m[max_lambda_star_counter_*sigma_power_+lambda_power_+1])))+
        
        log10(abs(table_gc_m[max_lambda_star_counter_*sigma_power_+lambda_power_+1]));
        
        gc_m_(i,j)=bb*(dum2-dum1)+dum1;
        
        gc_m_(i,j)=pow(10,gc_m_(i,j));
        
        //gp_p
        dum1=aa*(log10(abs(table_gp_p[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_]))
                 
                 -log10(abs(table_gp_p[max_lambda_star_counter_*sigma_power_+lambda_power_])))+
        
        log10(abs(table_gp_p[max_lambda_star_counter_*sigma_power_+lambda_power_]));
        
        dum2=aa*(log10(abs(table_gp_p[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_+1]))-
                 
                 log10(abs(table_gp_p[max_lambda_star_counter_*sigma_power_+lambda_power_+1])))+
        
        log10(abs(table_gp_p[max_lambda_star_counter_*sigma_power_+lambda_power_+1]));
        
        gp_p_(i,j)=bb*(dum2-dum1)+dum1;
        
        gp_p_(i,j)=-pow(10,gp_p_(i,j));
        
        //ge_p
        dum1=aa*(log10(abs(table_ge_p[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_]))-
                 
                 log10(abs(table_ge_p[max_lambda_star_counter_*sigma_power_+lambda_power_])))+
        
        log10(abs(table_ge_p[max_lambda_star_counter_*sigma_power_+lambda_power_]));
        
        dum2=aa*(log10(abs(table_ge_p[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_+1]))-
                 
                 log10(abs(table_ge_p[max_lambda_star_counter_*sigma_power_+lambda_power_+1])))+
        
        log10(abs(table_ge_p[max_lambda_star_counter_*sigma_power_+lambda_power_+1]));
        
        ge_p_(i,j)=bb*(dum2-dum1)+dum1;
        
        ge_p_(i,j)=-pow(10,ge_p_(i,j));
        
        //gc_p
        dum1=aa*(log10(abs(table_gc_p[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_]))-
                 
                 log10(abs(table_gc_p[max_lambda_star_counter_*sigma_power_+lambda_power_])))+
        
        log10(abs(table_gc_p[max_lambda_star_counter_*sigma_power_+lambda_power_]));
        
        dum2=aa*(log10(abs(table_gc_p[max_lambda_star_counter_*(sigma_power_+1)+lambda_power_+1]))-
                 
                 log10(abs(table_gc_p[max_lambda_star_counter_*sigma_power_+lambda_power_+1])))+
        
        log10(abs(table_gc_p[max_lambda_star_counter_*sigma_power_+lambda_power_+1]));
        
        gc_p_(i,j)=bb*(dum2-dum1)+dum1;
        
        gc_p_(i,j)=-pow(10,gc_p_(i,j));
        */
        
        ///for simplified system
        f_bar_(i,j)=1.;
        ge_bar_(i,j)=-50.*(1.+ C_bar(i,j));
        gc_bar_(i,j)=0.;
        gp_m_(i,j)=0.;
        ge_m_(i,j)=0.;
        gc_m_(i,j)=0.;
        gp_p_(i,j)=0.;
        ge_p_(i,j)=0.;
        gc_p_(i,j)=0.;
        //////////////////////////

            }//End of inner else

       }//End of i for loop
     }//End of j for  loop
 }//End of else

    return;
}

void Tablereading::calculate_A_B_f_RHS(Friendly_Vector &C_bar, Friendly_Vector &C0, Friendly_Vector &Cs, Friendly_Vector &lambda, Friendly_Vector &nondimensional_h, double *x_center, double *x_face){

     for(int j=0;j<gy_;j++)
     {
        for(int i=0;i<gx_;i++)
        {
           A1_(i,j)=Pe_/(2*lambda(i,j)*lambda(i,j))*gp_;

           B1_(i,j)=Pe_*ge_bar_(i,j);

           f1_(i,j)=-Pe_/(2*lambda(i,j)*lambda(i,j))*gc_bar_(i,j);

           A2_(i,j)=Pe_/(2*lambda(i,j)*lambda(i,j))*Cs(i,j)*(gp_+gp_p_(i,j)-gp_m_(i,j));
       
           B2_(i,j)=-(2*C_bar(i,j)+Cs(i,j))+ Pe_*Cs(i,j)*(ge_bar_(i,j)+ge_p_(i,j)-ge_m_(i,j));

           f2_(i,j)=-Pe_/(2*lambda(i,j)*lambda(i,j))*Cs(i,j)*(gc_bar_(i,j)+gc_p_(i,j)-gc_m_(i,j));
           
           Gp_m_(i,j)=Pe_/(2*lambda(i,j)*lambda(i,j))*gp_m_(i,j)*Cs(i,j);

           Ge_m_(i,j)=Pe_*ge_m_(i,j)*Cs(i,j);

           Gc_m_(i,j)=Pe_/(2*lambda(i,j)*lambda(i,j))*gc_m_(i,j)*Cs(i,j);

        }
     }


    ////x faces
    for(int j=0;j<gy_-2;j++)
    {
       for(int i=0;i<gx_-1;i++)
       {

           a1_x_(i,j)=(A1_(i,j+1)+A1_(i+1,j+1))/2.;

           b1_x_(i,j)=(B1_(i,j+1)+B1_(i+1,j+1))/2.;

           c1_x_(i,j)=(f1_(i,j+1)+f1_(i+1,j+1))/2.*(C0(i+1,j+1)-C0(i,j+1))/(x_center[i+1]-x_center[i]) ;

           A1_x_(i,j)=(A1_(i,j+1)+A1_(i+1,j+1))/2.*(nondimensional_h(i,j+1)+nondimensional_h(i+1,j+1) )/2.;

           B1_x_(i,j)=(B1_(i,j+1)+B1_(i+1,j+1))/2.*(nondimensional_h(i,j+1)+nondimensional_h(i+1,j+1) )/2.;

           f1_x_(i,j)=(f1_(i,j+1)+f1_(i+1,j+1))/2.*(C0(i+1,j+1)-C0(i,j+1))/(x_center[i+1]-x_center[i]) *(nondimensional_h(i,j+1)+nondimensional_h(i+1,j+1) )/2.;

           A2_x_(i,j)=(A2_(i,j+1)+A2_(i+1,j+1))/2.*(nondimensional_h(i,j+1)+nondimensional_h(i+1,j+1) )/2.;

           B2_x_(i,j)=(B2_(i,j+1)+B2_(i+1,j+1))/2.*(nondimensional_h(i,j+1)+nondimensional_h(i+1,j+1) )/2.;

           f2_x_(i,j)=(f2_(i,j+1)+f2_(i+1,j+1))/2.*(C0(i+1,j+1)-C0(i,j+1))/(x_center[i+1]-x_center[i])+(Cs(i,j+1)+Cs(i+1,j+1))/2.*(log(C0(i+1,j+1))-log(C0(i,j+1)))/(x_center[i+1]-x_center[i]);
           f2_x_(i,j)=f2_x_(i,j)*(nondimensional_h(i,j+1)+nondimensional_h(i+1,j+1) )/2.*0.0;

       }
    } 

    ////y faces
    for(int j=0;j<gy_-1;j++)
    {
       for(int i=0;i<gx_-2;i++)
       {
          a1_y_(i,j)=(A1_(i+1,j)+A1_(i+1,j+1))/2.;

          b1_y_(i,j)=(B1_(i+1,j)+B1_(i+1,j+1))/2.;

          c1_y_(i,j)=(f1_(i+1,j)+f1_(i+1,j+1))/2.*(C0(i+1,j+1)-C0(i+1,j))/dy_;


          A1_y_(i,j)=(A1_(i+1,j)+A1_(i+1,j+1))/2.*(nondimensional_h(i+1,j)+nondimensional_h(i+1,j+1))/2.;

          B1_y_(i,j)=(B1_(i+1,j)+B1_(i+1,j+1))/2.*(nondimensional_h(i+1,j)+nondimensional_h(i+1,j+1))/2.;

          f1_y_(i,j)=(f1_(i+1,j)+f1_(i+1,j+1))/2.*(C0(i+1,j+1)-C0(i+1,j))/dy_*(nondimensional_h(i+1,j)+nondimensional_h(i+1,j+1))/2.;

          A2_y_(i,j)=(A2_(i+1,j)+A2_(i+1,j+1))/2.*(nondimensional_h(i+1,j)+nondimensional_h(i+1,j+1))/2.;

          B2_y_(i,j)=(B2_(i+1,j)+B2_(i+1,j+1))/2.*(nondimensional_h(i+1,j)+nondimensional_h(i+1,j+1))/2.;

          f2_y_(i,j)=(f2_(i+1,j)+f2_(i+1,j+1))/2.*(C0(i+1,j+1)-C0(i+1,j))/dy_+(Cs(i+1,j)+Cs(i+1,j+1))/2.*(log(C0(i+1,j+1))-log(C0(i+1,j)))/dy_;
          f2_y_(i,j)=f2_y_(i,j)*(nondimensional_h(i+1,j)+nondimensional_h(i+1,j+1))/2.*0.0;
     }
   }
 
  ////RHS calculation
 
  //Boundary elements
  if(P_BC_Type_L_==1 &&  P_BC_Type_R_==1 && Phi_BC_Type_L_==0 && Phi_BC_Type_R_==0 )
  {
    for(int j=0;j<gy_-2;j++)
    {
       RHS_1_(0,j)=f1_x_(0,j);
       RHS_2_(0,j)=2*Phi_L_;

       RHS_1_(gx_-1,j)=f1_x_(gx_-2,j);
       RHS_2_(gx_-1,j)=2*Phi_R_;
    }
  }
  else //drichlet for both P and Phi
  {
    for(int j=0;j<gy_-2;j++)
    {
      RHS_1_(0,j)=2*P_L_;
      RHS_2_(0,j)=2*Phi_L_;

      RHS_1_(gx_-1,j)=2*P_R_;
      RHS_2_(gx_-1,j)=2*Phi_R_;
    }
  }

  //Rest of the elements
  for(int j=0;j<gy_-2;j++)
  {
    for(int i=1;i<gx_-1;i++)
    {
      RHS_1_(i,j)=1./(x_face[i]-x_face[i-1])*(f1_x_(i,j)-f1_x_(i-1,j)) + (f1_y_(i-1,j+1)-f1_y_(i-1,j))/(dy_);
     
      RHS_2_(i,j)=1./(x_face[i]-x_face[i-1])*(f2_x_(i,j)-f2_x_(i-1,j)) + (f2_y_(i-1,j+1)-f2_y_(i-1,j))/(dy_);
    }
  }  
   
return;
}

void Tablereading::calculate_dC(Friendly_Vector &C_bar, Friendly_Vector &C0, Friendly_Vector &P, Friendly_Vector &Phi,Friendly_Vector &nondimensional_h, Friendly_Vector &U, Friendly_Vector &V, Friendly_Vector &k_, Friendly_Vector &dC, double *x_center, double *x_face, double dt, double b){

   ////Velocity and flux in x
   for(int j=0;j<gy_-2;j++)
   {
     for(int i=0;i<gx_-1;i++)
     {
        U(i,j)=a1_x_(i,j)*(P(i+1,j+1)-P(i,j+1))/(x_center[i+1]-x_center[i])+ b1_x_(i,j)*(Phi(i+1,j+1)-Phi(i,j+1))/(x_center[i+1]-x_center[i])-c1_x_(i,j);

        adv_x_(i,j)=U(i,j)*(C_bar(i+1,j+1)+C_bar(i,j+1))/2. + (Gp_m_(i+1,j+1)+Gp_m_(i,j+1))/2.*(P(i+1,j+1)-P(i,j+1))/(x_center[i+1]-x_center[i]) + (Ge_m_(i+1,j+1)+Ge_m_(i,j+1))/2.*(Phi(i+1,j+1)-Phi(i,j+1))/(x_center[i+1]-x_center[i]) + (Gc_m_(i+1,j+1)+Gc_m_(i,j+1))/2.*(C0(i+1,j+1)-C0(i,j+1))/(x_center[i+1]-x_center[i]);
    
        diff_x_(i,j)=-(C_bar(i+1,j+1)+C_bar(i,j+1))/2.*(log(C0(i+1,j+1))-log(C0(i,j+1)))/(x_center[i+1]-x_center[i]);

        elec_mig_x_(i,j)=(C_bar(i+1,j+1)+C_bar(i,j+1))/2.*(Phi(i+1,j+1)-Phi(i,j+1))/(x_center[i+1]-x_center[i]);

        Flux_x_(i,j)=adv_x_(i,j)+diff_x_(i,j)+elec_mig_x_(i,j);

        Flux_x_(i,j)=Flux_x_(i,j)*(nondimensional_h(i+1,j+1)+nondimensional_h(i,j+1))/2.;
     }
  }  

  ////Velocity and flux in y
  for(int j=0;j<gy_-1;j++)
  {
     for(int i=0;i<gx_-2;i++)
     {
       V(i,j)=a1_y_(i,j)*(P(i+1,j+1)-P(i+1,j))/dy_ + b1_y_(i,j)*(Phi(i+1,j+1)-Phi(i+1,j))/dy_ -c1_y_(i,j);

       adv_y_(i,j)=V(i,j)*(C_bar(i+1,j+1)+C_bar(i+1,j))/2. + (Gp_m_(i+1,j+1)+Gp_m_(i+1,j))/2.*(P(i+1,j+1)-P(i+1,j))/dy_ + (Ge_m_(i+1,j+1)+Ge_m_(i+1,j))/2.*(Phi(i+1,j+1)-Phi(i+1,j))/dy_ + (Gc_m_(i+1,j+1)+Gc_m_(i+1,j))/2.*(C0(i+1,j+1)-C0(i+1,j))/dy_;

       diff_y_(i,j)=-(C_bar(i+1,j+1)+C_bar(i+1,j))/2.*(log(C0(i+1,j+1))-log(C0(i+1,j)))/dy_;

      elec_mig_y_(i,j)=(C_bar(i+1,j+1)+C_bar(i+1,j))/2.*(Phi(i+1,j+1)-Phi(i+1,j))/dy_;

      Flux_y_(i,j)=adv_y_(i,j)+diff_y_(i,j)+elec_mig_y_(i,j);

      Flux_y_(i,j)=Flux_y_(i,j)*(nondimensional_h(i+1,j+1)+nondimensional_h(i+1,j))/2.;
    }
 }

  /////k and dC calculation
  for(int j=0;j<gy_-2;j++)
  {
     for(int i=0;i<gx_-2;i++)
     {

      k_(i,j)=-dt*(Flux_x_(i+1,j)-Flux_x_(i,j))/(x_face[i+1]-x_face[i]) -dt*(Flux_y_(i,j+1)-Flux_y_(i,j))/dy_;

      k_(i,j)=k_(i,j)/nondimensional_h(i+1,j+1);

      dC(i,j)+=b*k_(i,j);

     }
  }
  
return;
}

