#pragma once
struct  time_frequency
{     	int N_f;
      	int N_f_thz;
	double f_max;
	double df;
	double omega_0; 
	double dt; 
	Eigen::ArrayXd omega_thz;
        Eigen::ArrayXd omega_b;
        Eigen::ArrayXd t;  
	time_frequency(int& N_f_,int& N_f_thz_,double& f_max_,double& df_,double& omega_0_): N_f(N_f_),N_f_thz(N_f_thz_),f_max(f_max_),df(df_),omega_0(omega_0_)
        {               double pi = M_PI;
			dt=1/(df*N_f);
			omega_b=Eigen::ArrayXd::Zero(N_f);
			t=Eigen::ArrayXd::Zero(N_f);
			omega_thz=Eigen::ArrayXd::Zero(N_f_thz);
			
           		for(int i=0;i<N_f_thz;i++)
    			omega_thz[i]= 2*pi*df*(i+1); 
 			
		     	for(int i=0;i<N_f;i++)
			{
   		 	omega_b[i]= omega_0-2*pi*N_f*df/2+2*pi*df*i;   
     		 	t[i]= -(N_f-1)/(2*N_f*df)+i*dt; 
			}

        }
	~time_frequency(){} 
 	
};



