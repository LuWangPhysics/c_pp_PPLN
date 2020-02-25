#include "/home/luwang/cpp_ppln/c++_ppln_3d/my_discritization/Include.hpp"
#include "/home/luwang/cpp_ppln/c++_ppln_3d/my_functions/my_print.hpp"

struct lowest_storage_RK
{   int N_f_thz;
    int N_f;
    Eigen::VectorXd r_0;
    Eigen::ArrayXcd ekx;
    Eigen::ArrayXcd ekx_t;
    int N_r;
    Eigen::MatrixXcd rk_ir_mpi;
    Eigen::MatrixXcd rk_thz_mpi;
    Eigen::Array3d ARK;
    Eigen::Array3d BRK;
    Eigen::MatrixXd M_r;
    double dr;
     std::complex<double> II={{0,1}};
    lowest_storage_RK(const int& N_f_thz_,const int& N_f_,const MESH::spacial_r& my_r_): N_f_thz(N_f_thz_),N_f(N_f_),r_0(my_r_.r_0),N_r(my_r_.N_r){dr=r_0(1)-r_0(0);}
    ~lowest_storage_RK(){}
    void assign();
    void method(const int& world_size,const int& n_thread,myfft& Fmy,const P_const& myconst, Eigen::MatrixXcd& U_ir,Eigen::MatrixXcd& P_ir, Eigen::MatrixXcd& U_thz,Eigen::MatrixXcd& P_thz, Material& LiNb ,const double& dz,std::vector<int>& iterations );


    void bc_metal(Eigen::MatrixXcd& U_ir, Eigen::MatrixXcd& U_thz);
    template<typename A>
    void bc_dieletric(Eigen::MatrixXcd& U_ir, Eigen::MatrixXcd& U_thz,const A& LiNb);    
    
    template<typename A>
    void bc_transparent(Eigen::MatrixXcd& U_ir, Eigen::MatrixXcd& U_thz,const A& ekx,const A& ekx_t);


};




void lowest_storage_RK::assign()
{ 
   rk_ir_mpi=Eigen::MatrixXcd::Zero(N_f,N_r);
   rk_thz_mpi=Eigen::MatrixXcd::Zero(N_f_thz,N_r);
   ARK<<0,(double)-5/9,(double)-153/128;
   BRK<<(double)1/3,(double)15/16,(double)8/15;

  
   

//*****************************************************************
//be careful E_t_irwith the 1, and end element
//finite difference calculation of the 1/r(drA(r)/dr) derivative, assuming
//A_0_i=A_1_i boundary condition
//*****************************************************************
   

   M_r=Eigen::MatrixXd::Zero(3,N_r);

    for(int i=0;i<N_r;i++)
    {
         int k=i+1;
          M_r(0,i)=-(r_0[k]+r_0[k-1])/(r_0[k]*2*pow(dr,2));
          M_r(1,i)=(r_0[k-1]+2*r_0[k]+r_0[k+1])/(r_0[k]*2*pow(dr,2));
          M_r(2,i)=-(r_0[k]+r_0[k+1])/(r_0[k]*2*pow(dr,2));
    }




}

void lowest_storage_RK::bc_metal(Eigen::MatrixXcd& U_ir, Eigen::MatrixXcd& U_thz)
{   
                    U_ir.col(0)=U_ir.col(1); //cylindrical symmetry at the boundary
                    U_thz.col(0)=U_thz.col(1);
                    //electrif field in metal is zero.
                    U_thz.col(N_r+1).setZero();
                    U_ir.col(N_r+1).setZero(); //metal
}
template<typename A>
void lowest_storage_RK::bc_dieletric(Eigen::MatrixXcd& U_ir, Eigen::MatrixXcd& U_thz,const A& LiNb)
{   
                    U_ir.col(0)=U_ir.col(1); //cylindrical symmetry at the boundary
                    U_thz.col(0)=U_thz.col(1);
                    //electrif field in dieletric E_in*eps_in=E_out*eps_out.
                    U_ir.col(N_r+1)=U_ir.col(N_r).array()*pow(LiNb.n_ir_b.array(),2)*r_0(N_r)/r_0(N_r+1);//dieletric
                    //remove the awful high refractive index at large terahertz frequency	
                    int thz_range =LiNb.n_thz.size()-10*LiNb.f_c_thz/LiNb.df;
		    
                    U_thz.col(N_r+1).array()=U_thz.col(N_r).array()*pow(LiNb.n_thz.array(),2)*r_0(N_r)/r_0(N_r+1);
                    U_thz.col(N_r+1).tail(thz_range).array()=0;
	
}
template<typename A>
void lowest_storage_RK::bc_transparent(Eigen::MatrixXcd& U_ir, Eigen::MatrixXcd& U_thz,const A& ekx,const A& ekx_t)
{   
                    U_ir.col(0)=U_ir.col(1); //cylindrical symmetry at the boundary
                    U_thz.col(0)=U_thz.col(1);
            
                    U_ir.col(N_r+1).array()=U_ir.col(N_r).array()*exp(-ekx)*r_0(N_r)/r_0(N_r+1);
                    U_thz.col(N_r+1).array()=U_thz.col(N_r).array()*exp(-ekx_t)*r_0(N_r)/r_0(N_r+1);
}

void lowest_storage_RK::method(const int& world_size,const int& n_thread,myfft& Fmy,const P_const& myconst, Eigen::MatrixXcd& U_ir,Eigen::MatrixXcd& P_ir, Eigen::MatrixXcd& U_thz,Eigen::MatrixXcd& P_thz, Material& LiNb ,const double& dz,std::vector<int>& iterations )
{ 
                for(int  m=0;m<3;m++)
                {
                
                      
                        
 			
                            for (int k = 0; k<N_r/(n_thread*world_size); k++)
                            { 
				#pragma omp parallel
  				{
                                int thread_rank=omp_get_thread_num();                                      //write in corresponding memory 
                                int N_t_i=thread_rank+world_size*n_thread*k;                               //number of colum in p_ir and p_thz;
             
                                //-----------------------------------
                                //calculate E field in time domain;
                                //--------------------------------------
                                E_time_domain(Fmy,thread_rank,(U_ir.data()+(N_t_i+1)*N_f),(U_thz.data()+(N_t_i+1)*N_f_thz),myconst);
                                //---------------------------------
                                //calculate polarisation
                                //-----------------------------------
                                P_ir_thz(thread_rank, Fmy,(P_ir.data()+N_f*N_t_i),(P_thz.data()+N_f_thz*N_t_i), myconst,LiNb.alpha_thz,(U_thz.data()+N_f_thz*(N_t_i+1)));   
                                //-------------------------------------
                                //add spacial derivative to the polarisation matrix
                                //-------------------------------------------------
                          
                                P_ir.col(N_t_i).array()+=(U_ir.block(0,N_t_i,N_f,3)*M_r.col(N_t_i)).array()*II/(2*LiNb.k_ir.array());
                                rk_ir_mpi.col(N_t_i)*=ARK(m);
                                rk_ir_mpi.col(N_t_i)+=dz*P_ir.col(N_t_i);
                       
                         
                                P_thz.col(N_t_i).array()+=(U_thz.block(0,N_t_i,N_f_thz,3)*M_r.col(N_t_i)).array()*LiNb.modify_c.array()*II/2;
                                rk_thz_mpi.col(N_t_i)*=ARK(m);
                                rk_thz_mpi.col(N_t_i)+=dz*P_thz.col(N_t_i);
                         
                                //count iteration to each thread;
                                iterations[thread_rank]++;                                                              
                         
                                //gather information from all threads to master thread 0;MPI_SEND AND RECEIVE AUTOMATICALLY block
/*
                if(thread_rank!=0){
                    MPI_Send((P_ir_mpi.data()+N_f*thread_rank),N_f,MPI_C_DOUBLE_COMPLEX,0,11, MPI_COMM_WORLD);
                    MPI_Send((P_thz_mpi.data()+N_f_thz*thread_rank),N_f_thz,MPI_C_DOUBLE_COMPLEX,0,22, MPI_COMM_WORLD);
                }
                else{
                    
                    P_ir.col(k*world_size)=P_ir_mpi.col(0);
                    P_thz.col(k*world_size)=P_thz_mpi.col(0);
                    for(int i=1;i<world_size;i++){
                        MPI_Recv((P_ir.data()+i*N_f+k*world_size*N_f),N_f,MPI_C_DOUBLE_COMPLEX,i,11, MPI_COMM_WORLD,&status);
                        MPI_Recv((P_thz.data()+i*N_f_thz+k*world_size*N_f_thz),N_f_thz,MPI_C_DOUBLE_COMPLEX,i,22, MPI_COMM_WORLD,&status);
                    }
                }  

                MPI_Allreduce(P_ir_mpi.data(),(P_ir.data()+k*world_size*N_f),N_f*world_size,MPI_C_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD);
                MPI_Allreduce(P_thz_mpi.data(),(P_thz.data()+k*world_size*N_f_thz),N_f_thz*world_size,MPI_C_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD);

*/ 		
				}
                            }
                          
                    
                    
                    //ekx=(log(U_ir.col(N_r).array())-log(U_ir.col(N_r-1).array()));
                    //ekx_t=(log(U_thz.col(N_r).array())-log(U_thz.col(N_r-1).array()));
                    #pragma omp parallel for
                    for (int i=0;i<N_r;i++)
                    {
                        U_ir.col(i+1)+=BRK(m)*rk_ir_mpi.col(i);
                        U_thz.col(i+1)+=BRK(m)*rk_thz_mpi.col(i);  
                    }
                    //---------------------------------------------------------------------------------------------
                    //boundary conditions ///number 1-N_r is the real matrix we want, 0 and N_r+1 are ghost points!
                    //---------------------------------------------------------------------------------------------
                  // bc_metal(U_ir, U_thz);
                    bc_dieletric(U_ir,U_thz,LiNb);
                    //bc_transparent(U_ir,U_thz,ekx,ekx_t);
                
                    
//                     save_array(N_f,1,ekx,"test_ir"+std::to_string(m)+".txt");               
//                     save_array( N_f_thz,1,ekx_t,"test_thz"+std::to_string(m)+".txt");
    
           }    
         
}
