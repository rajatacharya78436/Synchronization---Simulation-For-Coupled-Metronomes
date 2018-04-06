/* SL oscillator with coupling ***************/



#include<stdio.h>

#include<math.h>

#include<stdlib.h>

#include<limits.h>

#include<time.h>

#define pi 4.0*atan(1.0)

#define RANDMAX 2147483647

void rk41(double t,double x[],double h, double coupx[],double coupy[],double coupz[]);

double F1(double t,double x[],int i,double coupx[],double coupy[],double coupz[]);

int NSYS=100;

int M=3;

int n;

double dt = 0.01;

double C,r,x11,x12,x21,x22;

double  e,g,et,om,a,sig,rm;

int j,ip1,ip2,P;



void main()

{ 
    
  FILE *a_out0,*a_out1,*a_out2;
 
   a_out0=fopen("sl_t","w");
 
  a_out1=fopen("sl_s","w");
 
   a_out2=fopen("sl_f","w");
  
  
 long int i,k,na,l,ip,is,iq,ic,NP,im,ij,l1,j1,k1,k2,k3;
 
   na=300000;
    
   n=(NSYS*M);
  
   double M,M1; 

   double t = 0.0,x[n],z[n],coupx[n],coupy[n],coupz[n],temp1[NSYS],temp2[NSYS],temp3[NSYS],th[NSYS],count[NSYS],freq[NSYS],phase[NSYS],insp[NSYS],stphi[NSYS]; 

   double X,X1; 
 
  double min1,min2,min3,step,min,max;
  
   
/* initial condition */
   
   min=-1; max=1;
   
   
 for(i=1;i<=n;i++)
 
    {
     
  
      x[i]=min+(max-min)*(double)(rand())/(double)(RANDMAX);
 
  
  }


/** parameters **********/
  

   sig=4;
 
   rm=1.35;

   e=((sig*(rm*rm))/2);

   g=(-sig/4);

   et=0.1;
 
   om=3;

   a=0.8;



     for(j=1;j<=NSYS;j++)

     {  
     
      temp1[j]=0;
 
      temp2[j]=0;

      temp3[j]=0;
 
      phase[j]=0;
 
      insp[j]=0;
 
      stphi[j]=0;
 
      count[j]=0;

     }    

  
/* vary the value of r and C and find the dynamics ****/
  


    r=0.3;
  
    C=5.0;  

    P=(r*NSYS);
 
 

     for(i=1;i<=na;i++)
 
      { 
 
    
             /*   coupling denotes no of oscillator connected in each sie of the ring ********/
  
          for(l=1;l<=NSYS;l++)      
          {

           ip1=(l-P);
   
           ip2=(l+P); 
   
           coupx[(l-1)*3+1]=0.00;
 
           coupy[(l-1)*3+2]=0.00;

           coupz[(l-1)*3+3]=0.00;
 
                
                for(j=ip1;j<=ip2;j++) 
   
                {
           
                  if(j<=0)
            
                      { 
              
                       j1=(NSYS+j);
  
                      }
   
                  else if(j>NSYS)

                      {
   
                       j1=(j-NSYS);
  
                      }
        
                  else 
       
                      {
         
                       j1=j;
      
                      } 
        
   
        coupx[(l-1)*3+1]=coupx[(l-1)*3+1]+(x[(j1-1)*3+1]-x[(l-1)*3+1]);
  
           coupy[(l-1)*3+2]=coupy[(l-1)*3+2]+(x[(j1-1)*3+2]-x[(l-1)*3+2]);
  
           coupz[(l-1)*3+3]=coupz[(l-1)*3+3]+(x[(j1-1)*3+3]-x[(l-1)*3+3]);
 
                
                 }
 
 
           
           coupx[(l-1)*3+1]=(coupx[(l-1)*3+3]/(2*P));
   
           coupy[(l-1)*3+2]=(coupy[(l-1)*3+3]/(2*P));
  
           coupz[(l-1)*3+3]=(coupz[(l-1)*3+3]/(2*P)); 

         
 }  
       
        

       
	 
     rk41(t,x,dt,coupx,coupy,coupz);
       
          
   for(l1=1;l1<=NSYS;l1++)

               {
    
                 temp1[l1]=temp2[l1];
  
                 temp2[l1]=temp3[l1]; 
    
                 temp3[l1]=x[(l1-1)*3+1];
   
               }
           

   
                 
/***************snapshot ***************/
    

            if(i==299500)
         
             {
           
               for(j=1;j<=NSYS;j++)
   
                  { 
  
 
                   fprintf(a_out0,"%d\t %f\t  %f\t %f\n",j,x[(j-1)*3+1],x[(j-1)*3+2],x[(j-1)*3+3]);
  
                  }
   
         
    }
 
     
  	
             
                     /* space plot **************/
    
               
	
             if(i>=290000)
   
                {
       
                  for(j=1;j<=NSYS;j++)
    
   	            {
                
       	
                       fprintf(a_out1,"%f\t %d\t %f\t %f\t %f\n",t,j,x[(j-1)*3+1],x[(j-1)*3+2],x[(j-1)*3+3]);
     
      
       	    } 
     
 
                }     
 
                     
 /***Calculate Frequency*************************/
   
       
             if(i>=200000)
    
                {
    
                  for(j=1;j<=NSYS;j++)
  
                     {
    
                       if((temp2[j]>temp1[j])&&(temp2[j]>temp3[j]))  
    
                        { 
     
                               count[j] += 1;
   
       
                     }
   
                     }
      
                } 

                               /***********************************/
  
     
       
    
    
         
      t=t+dt;   
             
      
     }
 
        
     
   
      for(ic=1;ic<=NSYS;ic++)

             {
  
               freq[ic] =(2*pi*count[ic])/(double)(100000*dt);
 
               fprintf(a_out2,"%ld\t  %lf\t %lf\n",ic,phase[ic], freq[ic]);

             } 
 
  

}
   
       
          

 
/*********************-Equation------------------------*/


double F1(double t,double x[],int i,double coupx[],double coupy[],double coupz[])
 
  {
   
   
   
 
   for(j=1;j<=n;j=j+3)
  
      {
           

            if(i==j)
                 return((x[j+2]*x[j]-om*x[j+1]+2*pow(x[j],3)+2*x[j]*pow(x[j+1],2)-(e*x[j+1]*pow(x[j],2))-(e*pow(x[j+1],3))-pow(x[j],5)-(x[j]*pow(x[j+1],4))-(2*pow(x[j],3))*pow(x[j+1],2)-(g*x[j+1]*pow(x[j],4))-(g*pow(x[j+1],5))-(2*g*pow(x[j],2)*pow(x[j+1],3)))+C*coupx[j]);
 
            if(i==j+1) 
                 return((x[j+2]*x[j+1]+om*x[j]+(2*x[j+1]*pow(x[j],2))+2*pow(x[j+1],3)+(e*pow(x[j],3))+(e*x[j]*pow(x[j+1],2))-(x[j+1]*pow(x[j],4))-pow(x[j+1],5)-(2*pow(x[j],2)*pow(x[j+1],3))+g*pow(x[j],5)+g*x[j]*pow(x[j+1],4)+2*g*pow(x[j],3)*pow(x[j+1],2))+C*coupy[j+1]);   
   
            if(i==j+2)
                 return(et*(a-(pow(x[j],2)+pow(x[j+1],2))));
   
      }  








  }




/******************** --RK4 FUNCTION -------------------------------------------*/



void rk41(double t,double x[],double step,double coupx[],double coupy[],double coupz[])


  {
  
    double h = step/2.0;

    double t1[n],t2[n],t3[n],k1[n],k2[n],k3[n],k4[n];
   
 
   int i;
 
     for(i=1;i<=n;i++) 
          t1[i] = x[i] + 0.5*(k1[i] = step*F1(t,x,i,coupx,coupy,coupz));
  
     for(i=1;i<=n;i++) 
          t2[i] = x[i] + 0.5*(k2[i] = step*F1(t+h,t1,i,coupx,coupy,coupz));
 
    for(i=1;i<=n;i++)
          t3[i] = x[i] + (k3[i] = step*F1(t+h,t2,i,coupx,coupy,coupz));

     for(i=1;i<=n;i++)
          k4[i] =             dt*F1(t+step,t3,i,coupx,coupy,coupz);
 
 
    for(i=1;i<=n;i++) 
          x[i] += (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.0;
 
}

 

/***********************************************************************/

   



















