void many_lines_copa(Eigen::ArrayXcd& e_w,const double& F_peak_comb,const double& F_peak,const double& sigma_in,const int& n_gaussian,const Eigen::ArrayXd& omega_b,const Eigen::ArrayXd& n_ir_b,const double& n_ir_0,const double omega_0,const double tau,const double f_c_thz)
{
     //for 2 lines pump

    double pi=M_PI;
    double eps =8.85e-12;
    double c=3e8;
    double df=(omega_b[2]-omega_b[1])/(2*pi);
    double F_pump_total;
     Eigen::ArrayXd  delta_omega=omega_b-omega_0;


  
//-----------------------------------------------
// spectrum shape directly, calibrated by the total energy.
//one strong line and many weak lines
//------------------------------------------------
    	 e_w=exp(-square(-delta_omega*tau)/4);   							//strong line at 1030
	double f_separation=10e12;//f_c_thz;
	int    n_lines=2;
	double Energy_diff=sqrt(F_peak_comb/F_peak/n_lines);
	for(int m=0;m<n_lines;m++)
		{
			e_w+=Energy_diff*exp(- square((delta_omega+2*pi*(m*f_c_thz+f_separation))*tau)/4);  //spectrum shape 

		}


    

    F_pump_total=0.5*c*eps*(abs2(e_w)*n_ir_b).sum()*df; 

    double ratio=sqrt(F_pump_total/(F_peak+F_peak_comb));
    e_w.array()/=ratio;                                                                               //adjust the amplitude to conserve energy

}
