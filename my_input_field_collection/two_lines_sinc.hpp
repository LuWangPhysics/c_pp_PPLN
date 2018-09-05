Eigen::ArrayXcd my_sinc(const Eigen::ArrayXcd& ff )
{	
	 Eigen::ArrayXcd e_w(ff.size());
	for(int i=0;i<ff.size();i++)
	{	//e_w(i)=std::sin(ff(i))/ff(i);
		e_w(i)=boost::math::sinc_pi(ff(i));
	}	

	return e_w;
}
void two_lines_sinc(Eigen::ArrayXcd& e_w,const double& F_peak,const double& sigma_in,const int& n_gaussian,const Eigen::ArrayXd& omega_b,const Eigen::ArrayXd& n_ir_b,const double& n_ir_0,const double omega_0,const double tau,const double f_c_thz)
{
    double pi=M_PI;
    double eps =8.85e-12;
    double c=3e8;
    double df=(omega_b[2]-omega_b[1])/(2*pi);
    double F_pump_total;
    Eigen::ArrayXd  delta_omega=omega_b-omega_0;
    double tau_fwhm=tau*(sqrt(2*log(2)));
   //sin(pi*f*tau_fwhm)/f
    
     e_w=my_sinc((tau_fwhm/2)*delta_omega.array())+my_sinc((tau_fwhm/2)*(delta_omega-2*pi*f_c_thz).array());  //spectrum shape 

    F_pump_total=0.5*c*eps*(abs2(e_w)*n_ir_b).sum()*df; 

    double ratio=sqrt(F_pump_total/F_peak);
    e_w.array()/=ratio;                                                                               //adjust the amplitude to conserve energy
   
  
}
