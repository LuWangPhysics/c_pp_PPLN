
#ifndef POL_H
#define POL_H
#include<complex.h>
#include "/home/luwang/cpp_ppln/c++_ppln_3d/my_discritization/Include.hpp"
#include "/home/luwang/cpp_ppln/c++_ppln_3d/my_structures/Include.hpp"








//*************************************************
//only for even input array
//**************************************************
template<class T1>
  void myfftshift(T1 first1, T1 last1,T1 first2)
{
  while (first1<=last1) {
    std::swap(*first1, *first2);
    ++first1; ++first2;
  }
  
};


  void add_conj_flip(fftw_complex* E_thz ,fftw_complex* in,int N_f_thz)
{
                for(int i=0;i<N_f_thz;i++){
                in[i][0]=E_thz[i][0];
                in[i][1]=E_thz[i][1];
                in[i+N_f_thz][0]=E_thz[N_f_thz-i-1][0];
                in[i+N_f_thz][1]=-E_thz[N_f_thz-i-1][1];
                } 
  
};





//******************************************************
//calculate E field in time domain
// //********************************************************
void E_time_domain(myfft &Fmy,const int& world_rank,const std::complex<double>* U_ir,const std::complex<double>* U_thz,const P_const &myconst){
int N_f=myconst.N_f;
int N_f_thz=N_f/2;


            #pragma omp parallel for num_threads(myconst.n_nested)
               for(int i=0;i<N_f;i++){
                 
                   Fmy.E_ir[i+N_f*world_rank][0]=(U_ir[i]*conj(myconst.phase_ir(i))).real();
                   Fmy.E_ir[i+N_f*world_rank][1]=(U_ir[i]*conj(myconst.phase_ir(i))).imag();
                }
                myfftshift(&Fmy.E_ir[N_f*world_rank],&Fmy.E_ir[N_f*world_rank+N_f_thz-1],&Fmy.E_ir[N_f*world_rank+N_f_thz]);
                
                fftw_execute_dft(Fmy.ifft,&Fmy.E_ir[N_f*world_rank], &Fmy.E_t_ir[N_f*world_rank]); 


                
                //calculate E_t_thz
           #pragma omp parallel for num_threads(myconst.n_nested)

                
                for(int i=0;i<N_f_thz;i++){
                    Fmy.E_thz[1+i+(N_f_thz+1)*world_rank][0]=(U_thz[i]*conj(myconst.phase_thz(i))).real();
                    Fmy.E_thz[1+i+(N_f_thz+1)*world_rank][1]=(U_thz[i]*conj(myconst.phase_thz(i))).imag();
                }
                Fmy.E_thz[(N_f_thz+1)*world_rank][0]=0;
                Fmy.E_thz[(N_f_thz+1)*world_rank][1]=0;
                 fftw_execute_dft_c2r(Fmy.c2r,&Fmy.E_thz[(N_f_thz+1)*world_rank],&Fmy.E_t_thz[N_f*world_rank]);

};



//******************************************************
//calculate polarisation
//********************************************************

void  P_ir_thz(const int& rank,myfft& Fmy,std::complex<double>* P_ir,std::complex<double>* P_thz,const P_const &myconst,const Eigen::ArrayXd& alpha_thz,const std::complex<double>* U_thz){
   
    
double df=myconst.df;
int N_f=myconst.N_f;
int N_f_thz=N_f/2; 
int NN=N_f*rank;

std::complex<double> tmp(0,0);

        //***************************
        //IR spm calculation
        //****************************
           #pragma omp parallel for num_threads(myconst.n_nested)
                for(int i=0;i<N_f;i++){
                Fmy.I_t_ir[i+NN]=(pow(Fmy.E_t_ir[i+NN][0],2)+pow(Fmy.E_t_ir[i+NN][1],2))*pow(df,2);
                Fmy.in[i+NN][0]=Fmy.I_t_ir[i+NN]*df*Fmy.E_t_ir[i+NN][0];
                Fmy.in[i+NN][1]=Fmy.I_t_ir[i+NN]*df*Fmy.E_t_ir[i+NN][1];
                }          
                 
                fftw_execute_dft(Fmy.fft,&Fmy.in[NN],&Fmy.out[NN]);
                myfftshift(&Fmy.out[NN],&Fmy.out[NN+N_f_thz-1],&Fmy.out[NN+N_f_thz]);
             #pragma omp parallel for private(tmp) num_threads(myconst.n_nested) 
                for(int i=0;i<N_f;i++){
                    tmp.real(Fmy.out[NN+i][0]);
                    tmp.imag(Fmy.out[NN+i][1]);
                    P_ir[i]=(1/(N_f*df))*myconst.phase_ir(i)*myconst.spm(i)*tmp;
                }

        
       // ***************************
        //IR dfg calculation
        //****************************
         #pragma omp parallel for num_threads(myconst.n_nested)
               for(int i=0;i<N_f;i++){
                Fmy.in[i+NN][0]=(Fmy.E_t_ir[i+NN][0]*Fmy.E_t_thz[i+NN+1])*pow(df,2);
                Fmy.in[i+NN][1]=(Fmy.E_t_ir[i+NN][1]*Fmy.E_t_thz[i+NN+1])*pow(df,2);
                }
                
              
                fftw_execute_dft(Fmy.fft,&Fmy.in[NN],&Fmy.out[NN]);
                myfftshift(&Fmy.out[NN],&Fmy.out[NN+N_f_thz-1],&Fmy.out[NN+N_f_thz]);
          #pragma omp parallel for private(tmp) num_threads(myconst.n_nested)
                for(int i=0;i<N_f;i++){
                    tmp.real(Fmy.out[NN+i][0]);
                    tmp.imag(Fmy.out[NN+i][1]);
                    P_ir[i]+=(1/(N_f*df))*myconst.phase_ir(i)*myconst.dfg(i)*tmp*myconst.chi_2;
                }

            
       //***************************
       //THz dfg calculation
      // ****************************
             
                fftw_execute_dft_r2c(Fmy.r2c,&Fmy.I_t_ir[NN],&Fmy.out[NN]);
           #pragma omp parallel for private(tmp) num_threads(myconst.n_nested)
               for(int i=0;i<N_f_thz;i++){
                    tmp.real(Fmy.out[NN+i+1][0]);
                    tmp.imag(Fmy.out[NN+i+1][1]);
                    //double n_adjust=0.95;//for thz=0.5THz. 0.3THz is the same
                    P_thz[i]=(1/(N_f*df))*myconst.phase_thz(i)*myconst.thz(i)*tmp*myconst.chi_2+0.5*alpha_thz(i)*U_thz[i];
			              
		}


};



#endif
