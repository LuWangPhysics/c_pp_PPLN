#pragma once

struct  spacial_z
{      
        int N_z; 
        double dz;
	Eigen::VectorXd z_0;
        Eigen::VectorXd chi_2_z;

	spacial_z(const double& dppln,const  int&Nppln, const int& nz_ppln,const double& chi_2)
        {   
         
            dz=(dppln/nz_ppln);
            N_z=Nppln*nz_ppln;
            chi_2_z=Eigen::VectorXd::Zero(N_z);
            z_0=Eigen::VectorXd::Zero(N_z);
            
            for(int i=0;i<Nppln;i++)
            {
                for(int j=0;j<nz_ppln;j++)
                {
                    chi_2_z[i*nz_ppln+j]=pow(-1,i)*chi_2;
                    z_0[i*nz_ppln+j]=(i*nz_ppln+j)*dz;
                }

            }
        };
	~spacial_z(){} ;
 	
};
