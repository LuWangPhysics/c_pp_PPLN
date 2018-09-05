#pragma once

struct Material
{    	
	Eigen::ArrayXd n_ir_b;
	Eigen::ArrayXd n_thz;
	Eigen::ArrayXd modify_c;
	Eigen::ArrayXd alpha_thz;
        Eigen::ArrayXd k_ir;
        Eigen::ArrayXd k_thz;
        double f_c_thz;
        double n_ir_0;
        double n_g;
        double n_thz_c;
	double df;
	Material(const MESH::time_frequency& tf, const double& f_c_thz_, const double& df_)
	{	f_c_thz=f_c_thz_;
                double c=3e8;
		df=df_;
		n_ir_b=Eigen::ArrayXd::Zero(tf.N_f);
                n_thz=Eigen::ArrayXd::Zero(tf.N_f_thz);
                modify_c=Eigen::ArrayXd::Zero(tf.N_f_thz);
                alpha_thz=Eigen::ArrayXd::Zero(tf.N_f_thz);
                k_thz=Eigen::ArrayXd::Zero(tf.N_f_thz);
                k_ir=Eigen::ArrayXd::Zero(tf.N_f);
                
		LiNb_80k_ir(tf.omega_b,n_ir_b);
                LiNb_80k_thz(tf.omega_thz,n_thz,alpha_thz);
                for(int i=0;i<tf.N_f_thz;i++)
		{
                    modify_c[i]= ((c/(n_thz[i]*tf.omega_thz[i]))>4e-5?4e-5:(c/(n_thz[i]*tf.omega_thz[i])));   //2e-5 for 0.3Thz  
		}
                n_ir_0=n_ir_b[(int) tf.N_f/2];
                n_g=Gro_v(tf.omega_b,n_ir_b);
                n_thz_c=n_thz[(int)(f_c_thz/tf.df+1)];
                k_ir=tf.omega_b*n_ir_b/c; 
                k_thz=tf.omega_thz*n_thz/c; 
	}
	~Material(){};
       
};











