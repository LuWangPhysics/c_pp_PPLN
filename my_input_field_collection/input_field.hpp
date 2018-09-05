#ifndef INPUT_FIELD_H
#define INPUT_FIELD_H
#include "/home/luwang/cpp_ppln/c++_ppln_3d/my_discritization/name_space.hpp"
#include "/home/luwang/cpp_ppln/c++_ppln_3d/my_material_function/LiNb_final.hpp"
#include "/home/luwang/cpp_ppln/c++_ppln_3d/my_input_field_collection/two_lines.hpp"
#include "/home/luwang/cpp_ppln/c++_ppln_3d/my_input_field_collection/two_lines_sinc.hpp"
#include "/home/luwang/cpp_ppln/c++_ppln_3d/my_input_field_collection/many_lines_copa.hpp"
#include "/home/luwang/cpp_ppln/c++_ppln_3d/my_input_field_collection/spatial_phase.hpp"
void input_field(const MESH::spacial_r& my_r,const double& F_peak_comb,const double& F_peak,const double& sigma_in,const int& n_gaussian,const MESH::time_frequency& tf,Eigen::MatrixXcd  &U_ir,const double tau,const Material& LiNb){
    
   const double pi = M_PI;
  const std::complex<double> II={{0,1}};
  double df=(tf.omega_b[2]-tf.omega_b[1])/(2*pi);
    Eigen::ArrayXcd e_w(tf.N_f);
   //many_lines_copa(e_w,F_peak_comb,F_peak, sigma_in, n_gaussian, tf.omega_b, LiNb.n_ir_b, LiNb.n_ir_0, tf.omega_0, tau,LiNb.f_c_thz);
    two_lines(e_w,F_peak, sigma_in, n_gaussian, tf.omega_b, LiNb.n_ir_b, LiNb.n_ir_0, tf.omega_0, tau,LiNb.f_c_thz);
    //two_lines_sinc(e_w,F_peak, sigma_in, n_gaussian, tf.omega_b, LiNb.n_ir_b, LiNb.n_ir_0, tf.omega_0, tau,LiNb.f_c_thz);
    
//------------------------------------------------------------------------------------------
//initializing the input filed and the rk storage matrix
//Eigen save array in colum major
//E=exp(-(r^2/2sigma^2)^n_gaussian)
//------------------------------------------------------------------------------------------
    int n_half=tf.N_f/2+2*log(2)/tau/df;			// 0:n_half covers one strong line
 
    for(int i=0;i<my_r.N_r+2;i++)
    {   //(1+0.05*cos(2.0*pi*4.0*i*my_r.dr/my_r.r_0(my_r.N_r)))   input intensity sin
	double phase_n=spatial_phase(my_r,i);        
	U_ir.col(i).array()+=e_w*exp(-pow(pow(my_r.r_0[i]/sigma_in,2)/2,n_gaussian));     //out product of frequency and position  
	//half intensity sign	
	//U_ir.col(i).head(n_half).array()+=U_ir.col(i).head(n_half).array()*0.05*cos(2.0*pi*4.0*i*my_r.dr/my_r.r_0(my_r.N_r));    
                        
//-----------------------------------------------------------------------------------
//add curved phase front, the lower frequency line has a spacial dependent phase.
//---------------------------------------------------------------------------------
  	 
	//one has sin phase
  	U_ir.col(i).head(n_half).array()=U_ir.col(i).head(n_half).array()*exp(II*phase_n*0.05);     
	//U_ir.col(i).array()=U_ir.col(i).array()*exp(II*phase_n);     
    }


  
  
 
    
}

 #endif
