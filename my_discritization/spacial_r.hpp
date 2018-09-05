#pragma once
struct  spacial_r
{       double dr;
	Eigen::VectorXd r_0;
	int N_r; 
	spacial_r(const double& dr_,const int& N_r_ ): dr(dr_),N_r(N_r_)
        {   
            r_0=Eigen::VectorXd::Zero(N_r+2);
	
	
            for (int j=0;j<N_r+2;j++)r_0[j]=j*dr+3e-4;     
	
            //avoid 0 point!!! Coordinate array radius direction  r
            //create 2 ghost poins
        }
	~spacial_r(){} 
 	
};



