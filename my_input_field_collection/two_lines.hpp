void two_lines(Eigen::ArrayXcd& e_w,const double& F_peak,const double& sigma_in,const int& n_gaussian,const Eigen::ArrayXd& omega_b,const Eigen::ArrayXd& n_ir_b,const double& n_ir_0,const double omega_0,const double tau,const double f_c_thz)
{
     //for 2 lines pump

    double pi=M_PI;
    double eps =8.85e-12;
    double c=3e8;
    double df=(omega_b[2]-omega_b[1])/(2*pi);
    double F_pump_total;
     Eigen::ArrayXd  delta_omega=omega_b-omega_0;

//calculated from the energy point of view. everything consistent checked.
   // double e0=sqrt(pow(2,1.0/n_gaussian)*sqrt(2)*energy/(pi*pow(sigma_in,2)*sqrt(pi)*c*eps*n_ir_0*tau*std::tgamma(1.0*(n_gaussian+1)/n_gaussian)));
     //e_w= tau*e0*sqrt(pi)*exp(-square(-delta_omega*tau)/4)+tau*e0*sqrt(pi)*exp(- square((delta_omega-2*pi*f_c_thz)*tau)/4); 
  
//-----------------------------------------------
// spectrum shape directly, calibrated by the total energy.
//------------------------------------------------
     e_w=exp(-square(-delta_omega*tau)/4)+exp(- square((delta_omega-2*pi*f_c_thz)*tau)/4);  //spectrum shape 

    

    F_pump_total=0.5*c*eps*(abs2(e_w)*n_ir_b).sum()*df; 

    double ratio=sqrt(F_pump_total/F_peak);
    e_w.array()/=ratio;                                                                               //adjust the amplitude to conserve energy

}
