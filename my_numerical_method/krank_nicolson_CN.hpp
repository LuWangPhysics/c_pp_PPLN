struct krank_nicolson_CN 
{ 
    int N_f_thz;
    int N_f;
    Eigen::VectorXd r_0;
    Eigen::ArrayXcd kx;
    int N_r;
    Eigen::MatrixXcd rk_ir_mpi;
    Eigen::MatrixXcd rk_thz_mpi;
    Eigen::MatrixXd M_l;
    Eigen::MatrixXd M_r;
    double dr;
     std::complex<double> II={{0,1}};
     krank_nicolson_CN (){}
    ~krank_nicolson_CN (){}
template<typename A,typename B>
void assign(const A& my_r,const B& my_z);
// void method(const int& world_size,const int& n_thread,myfft& Fmy,const P_const& myconst, Eigen::MatrixXcd& U_ir,Eigen::MatrixXcd& P_ir, Eigen::MatrixXcd& U_thz,Eigen::MatrixXcd& P_thz, Material& LiNb ,const double& dz,std::vector<int>& iterations );

    
};
template<typename A,typename B>
void krank_nicolson_CN ::assign(const A& my_r,const B& my_z)
{ 
   

//*****************************************************************
//be careful E_t_irwith the 1, and end element
//finite difference calculation of the 1/r(drA(r)/dr) derivative, assuming
//two ghost points A_0_i A_N_i at boundary condition
//*****************************************************************
   Eigen::ArrayXd r_min=my_r.r_0.segment(0,my_r.N_r)+my_r.r_0.segment(1,my_r.N_r);
   Eigen::ArrayXd r_plus=my_r.r_0.segment(1,my_r.N_r)+my_r.r_0.segment(2,my_r.N_r);
  double c_cons=II*my_z.dz/4/pow(my_r.dr,2);
   M_r=Eigen::MatrixXd::Zero(my_r.N_r,my_r.N_r);
   M_l=Eigen::MatrixXd::Zero(my_r.N_r,my_r.N_r);
//     for(int i=0;i<N_r;i++)
//     {
//          int k=i+1;
//           M_r(0,i)=-(r_0[k]+r_0[k-1])/(r_0[k]*2*pow(dr,2));
//           M_r(1,i)=(r_0[k-1]+2*r_0[k]+r_0[k+1])/(r_0[k]*2*pow(dr,2));
//           M_r(2,i)=-(r_0[k]+r_0[k+1])/(r_0[k]*2*pow(dr,2));
//     }




}
