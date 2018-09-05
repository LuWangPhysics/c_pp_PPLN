//all my structures are defined here.
#ifndef MY_STRUCTURES_H
#define MY_STRUCTURES_H
#include<complex.h>
#include "/home/luwang/cpp_ppln/c++_ppln_3d/my_functions/my_print.hpp"


struct P_const
 {
        Eigen::VectorXcd spm; 
        Eigen::VectorXcd phase_ir;
        Eigen::VectorXcd phase_thz; 
        Eigen::VectorXcd dfg;
        Eigen::VectorXcd thz;
      


        
        double df;
        int N_f;
        double chi_2;
        int n_nested;
        std::complex<double> II={{0,1}};
        double c=3e8;
        double eps=8.85e-12;

        P_const(const MESH::time_frequency& tf_,const Material& LiNb_,const double& n_2_) 
	{	 

                 N_f=tf_.N_f;
       		 spm=-0.5*II*eps*LiNb_.n_ir_b*n_2_*tf_.omega_b;
      		 dfg=-0.5*II*tf_.omega_b.square()/(pow(c,2)*LiNb_.k_ir);
      		 thz=-0.5*II*tf_.omega_thz.square()/pow(c,2)*LiNb_.modify_c;  
	}
        ~P_const(){}
        void phase(const double& z_m,const double stage_phase,const double& dz,const Material& LiNb_);
      
}; 

void P_const::phase(const double& z_m,const double stage_phase,const double& dz,const Material& LiNb_)
{       phase_ir=(II*LiNb_.k_ir*(z_m+0.5*dz+stage_phase)).exp();
        phase_thz=(II*LiNb_.k_thz*(z_m+0.5*dz)).exp();
}





struct Output{
    double c=3e8;
    double pi=M_PI;
    double eps =8.85e-12;
    //double df;
   

    Eigen::ArrayXd ir_spectrum;
    Eigen::ArrayXd thz_spectrum;
    Eigen::ArrayXd Energy_ir;
    Eigen::ArrayXd Energy_thz;
    Eigen::ArrayXd ir_spacial;
    Eigen::ArrayXd thz_spacial;
    Eigen::MatrixXcd E_thz_final;
    std::string save_path_;
        Output(std::string save_path) {save_path_=save_path;}
     //Output()
       
        ~Output(){}
        template<typename A,typename B>
        void spectrum_energy(const double df,const Eigen::MatrixXcd& U_ir, const Eigen::MatrixXcd& U_thz,int j, const A& my_r, const B& LiNb);
         template<typename A,typename B,typename C,typename D>
        void save(const int N_sta,const A& my_z,const B& Fmy,const Eigen::MatrixXcd& U_ir, const Eigen::MatrixXcd& U_thz,const MESH::time_frequency& tf,const Eigen::VectorXcd& phase_thz, const C& my_r, const D& LiNb);
};



template<typename A,typename B>
void Output::spectrum_energy(const double df,const Eigen::MatrixXcd& U_ir, const Eigen::MatrixXcd& U_thz,int j, const A& my_r, const B& LiNb)
{       ir_spectrum=pi*c*LiNb.n_ir_0*eps*((U_ir.cwiseAbs2()).matrix()*my_r.r_0).array()*my_r.dr;
        thz_spectrum=pi*c*LiNb.n_thz*eps*((U_thz.cwiseAbs2()).matrix()*my_r.r_0).array()*my_r.dr;
        Energy_ir(j)=ir_spectrum.sum()*df;
        Energy_thz(j)=thz_spectrum.sum()*df; 
};
 template<typename A,typename B,typename C,typename D>
void Output::save(const int N_sta,const A& my_z,const B& Fmy,const Eigen::MatrixXcd& U_ir, const Eigen::MatrixXcd& U_thz,const MESH::time_frequency& tf,const Eigen::VectorXcd& phase_thz, const C& my_r, const D& LiNb)
{
         

        
        ir_spacial=0.5*LiNb.n_ir_0*eps*c*U_ir.cwiseAbs2().colwise().sum()*tf.df;
        thz_spacial=0.5*LiNb.n_thz_c*eps*c*U_thz.cwiseAbs2().colwise().sum()*tf.df;
        
        //get the pahse of the E_thz
        #pragma omp parallel for
        for (int i=0;i<my_r.N_r;i++)
        { 
           E_thz_final.col(i).array()=U_thz.col(i+1).array()*(phase_thz.conjugate().array());
        }

          
        save_array(my_r.N_r,my_r.r_0.data(),save_path_+"l_x.txt" );
        save_array(tf.N_f,tf.omega_b.data(),save_path_+"frequency_ir.txt" );
        save_array(tf.N_f_thz,tf.omega_thz.data(),save_path_+"frequency_thz.txt" );
        save_array(my_z.N_z-1,my_z.z_0.data(),save_path_+"l_z.txt" );
        save_array(tf.N_f,Fmy.E_t_thz,save_path_+"E_t_thz_"+std::to_string(N_sta)+".txt");
        save_array(tf.N_f,ir_spectrum.data(),save_path_+"ir_spectrum_"+std::to_string(N_sta)+".txt" );
        save_array(tf.N_f_thz,thz_spectrum.data(),save_path_+"thz_spectrum_"+std::to_string(N_sta)+".txt" );
        save_array(my_z.N_z-1,Energy_ir.data(),save_path_+"ir_energy_"+std::to_string(N_sta)+".txt" );
        save_array(my_z.N_z-1,Energy_thz.data(),save_path_+"thz_energy_"+std::to_string(N_sta)+".txt" );
        save_array(my_r.N_r+2,ir_spacial.data(),save_path_+"ir_spacial_"+std::to_string(N_sta)+".txt");  
        save_array(my_r.N_r+2,thz_spacial.data(),save_path_+"thz_spacial_"+std::to_string(N_sta)+".txt");  
        save_array(tf.N_f_thz,my_r.N_r,E_thz_final,save_path_+"E_thz_f_"+std::to_string(N_sta)+".txt" );
        save_array(tf.N_f,my_r.N_r,U_ir,save_path_+"E_ir_f_"+std::to_string(N_sta)+".txt" );

};

#endif
