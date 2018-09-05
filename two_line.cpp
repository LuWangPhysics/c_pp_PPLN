#include <mpi.h>
#include <iostream>
#include <stdlib.h>                                                                 // include random number
#include <sstream>
#include <fstream>                                                                  //read and write file
#include <math.h>                                                                   // included to use the "fabs" function
#include <Eigen/Dense>
#include <typeinfo>
#include <complex>
#include <fftw3.h>
#include <algorithm>                                                                //include swap function
#include <functional>                                                               // templete of cwise operation of vector
#include <omp.h>
#include <algorithm>
#include <Eigen/LU>
#include <string>	
#include <cmath>
#include<time.h>
#include <array>
#include <vector>
#include <thread>
#include <mutex>
#include </home/luwang/cpp_ppln/c++_ppln_3d/boost_1_65_0/boost/math/special_functions/sinc.hpp>
#include "/home/luwang/cpp_ppln/c++_ppln_3d/my_discritization/Include.hpp"
#include "/home/luwang/cpp_ppln/c++_ppln_3d/my_material_function/Include.hpp"
#include "/home/luwang/cpp_ppln/c++_ppln_3d/my_functions/Include.hpp"
#include "/home/luwang/cpp_ppln/c++_ppln_3d/my_structures/Include.hpp"                                      
#include "/home/luwang/cpp_ppln/c++_ppln_3d/my_numerical_method/Include.hpp"
#include "/home/luwang/cpp_ppln/c++_ppln_3d/my_input_field_collection/Include.hpp"



//---------------------------------------------------------
//define global variables. all constants
//---------------------------------------------------------

int main(int argc, char **argv){ 

srand(time(NULL));                                                                    //seed the random number
const double pi = M_PI;   
const double c=3e8;                                                                  //Speed of light  
 //const double eps =8.85e-12;                                                       //Free space permittivity  


    
mytime time_count;
time_count.start();
std::cout.precision(10);                                                              //change the default prcision std   
    
//***************************
//MPI operation
//******************
int world_size;
int world_rank;
my_mpi_initial(argc,argv, world_size, world_rank);

//***********************
//openmp
//***************************
omp_set_nested(1);
omp_set_dynamic(0);
int n_thread=50;//49 5MM
int n_nested=1;

omp_set_num_threads(n_thread);


//--------------------------------------------------------------------
//Pump pulse Parameters
//--------------------------------------------------------------------

double sigma_ch=(5.0/3.0);
int n_gaussian=5;                                                                    //super gaussian order
double sigma_in=sigma_ch*3e-3;                                                       // Input beam waist
//10/0.02a 3b 5c 7d 11e 13f, where e+f is either 0 or 1
double lambda_0    = 1030e-9;
double tau_fwhm    = 150e-12;//50*N_tau*1e-12;                                        //Input pulse fwhm
double tau = tau_fwhm/(sqrt(2*log(2)));                                              //Converting from FWHM to 1/e^2 value. e^(-t^2/(2*tau^2))
double f_max       = 75e12;    //4stages 120                                                       //2*f_max =frequency range ,120
double df          = 0.0002e12;//0.0002e12 500ps 0.0004e12 150ps;                                          //frequency step size
double omega_0 =  2*pi*c/lambda_0;
double f_c_thz=0.3e12;                                                               //define thz center frequency 
//peak fluence is fixed to damage and the energy is scaled with respect to peak fluence.
double energy=2*pi*pow(sigma_in,2)*std::tgamma(1.0*(n_gaussian+1.0)/n_gaussian)*1e5*sqrt(tau_fwhm/2/1e-8)/pow(2,(1.0/n_gaussian));   //pulse energy J
double F_peak=pow(2,(1.0/n_gaussian))*energy/(2*pi*pow(sigma_in,2)*std::tgamma(1.0*(n_gaussian+1.0)/n_gaussian)); //peak fluence
double energy_comb=100e-3;
double sigma_comb=sigma_in;
double F_peak_comb=pow(2,(1.0/n_gaussian))*energy_comb/(2*pi*pow(sigma_comb,2)*std::tgamma(1.0*(n_gaussian+1.0)/n_gaussian)); //peak fluence

//-----------------------------------------------------------------------
//spacial discritization
//---------------------------------------------------------------------
double dr =sigma_ch*2e-5;  //1.0e-4
double crystal_size=sigma_ch*5e-3;
int N_r;
//rescale the mesh point to fit the thread number
adjust_mesh_r(N_r, world_size, n_thread, crystal_size, dr);

MESH::spacial_r my_r( dr, N_r);

/*%-------------------------------------------------------------------
%Material Parameters
%------------------------------------------------------------------*/                                                                 
double chi_2   =  2*168e-12;                                                         // 2nd order nonlinearity 
double n2      =  1.25e-19;                                                          // Intensity dependent refractive index


//----------------------------------------------------------------------
//define THz and IR frequency +refractive index matrix
//----------------------------------------------1.0----------------------
//IR****************************************************************
int N_f,N_f_thz;
int   NN_f=2*(f_max/df);
N_f=NN_f-(NN_f%4);                                                                   // N_f , N_f_thz even
N_f_thz=(N_f)/2;

MESH::time_frequency tf(N_f, N_f_thz, f_max,df,omega_0);
Material LiNb(tf,f_c_thz,df);

//----------------------------------------------------------------------
//define phase matched period
//----------------------------------------------------------------------

int Nppln=145; //145 for 150ps 190for 500ps stages                                                    // Total number of poling periods                                                     
int nz_ppln=75; //35                                                                 //number of steps in one period
double dppln=c/(2*LiNb.f_c_thz*(LiNb.n_thz_c-LiNb.n_g));                              // Thickness of each poling period
int   Num_stage=1;
MESH::spacial_z my_z(dppln,Nppln,nz_ppln,chi_2);

//std::vector<double> l_sta={0.8,0.8,0.45,0.45};//500ps
std::vector<double> l_sta={0.92,0.75,0.45,0.45};                                             //smaller than 1, crystal length 
//std::vector<double> l_sta={0.92,0.96,0.96,0.76}; //copensate dispersion
double stage_phase=0;
int   N_zz=0;

//--------------------------------------------------          
//define input e_field and the derivative matrix
//e field spectrum vertical direction in each x position horizontal
//--------------------------------------------------------------------
Eigen::MatrixXcd U_thz=Eigen::MatrixXcd::Zero(N_f_thz,N_r+2);
Eigen::MatrixXcd U_ir=Eigen::MatrixXcd::Zero(N_f,N_r+2);
const double delta_safe=1e-16;

U_ir.array()+=delta_safe;
U_thz.array()+=delta_safe;
input_field(my_r,F_peak_comb,F_peak,sigma_in,n_gaussian,tf,U_ir,tau,LiNb);


//---------------------------------------------------
//define calculation constants
//-----------------------------------------------------
P_const myconst(tf,LiNb, n2);

myconst.df=df;
myconst.n_nested=n_nested;

Eigen::MatrixXcd P_ir=Eigen::MatrixXcd::Zero(N_f,N_r);
Eigen::MatrixXcd P_thz=Eigen::MatrixXcd::Zero(N_f_thz,N_r);
std::string save_path="/data/netapp/luwang/2D_150ps_4stages_effSaturate/2D_with_diffraction/2_lines/wave_front_5mm_one_has_sinPhase/0.05_2pi_moudulaton/";

save_array(tf.N_f,my_r.N_r,U_ir,save_path+"E_ir_input.txt" );
Output mydata(save_path);

mydata.Energy_ir=Eigen::ArrayXd::Zero(my_z.N_z-1);
mydata.Energy_thz=Eigen::ArrayXd::Zero(my_z.N_z-1);
mydata.E_thz_final=Eigen::MatrixXcd::Zero(N_f_thz,N_r);


//--------------------------------------------------------------
//numerical method
//---------------------------------------------------------
//**************************************************
//RK method initialization//lowest storage rk method
//**************************************************
lowest_storage_RK RK( N_f_thz,N_f, my_r);
RK.assign();
//krank_nicolson_CN CN;
//CN.assign(my_r,my_z);

//--------------------------------------------------------------
//FFTW initialization in my myfft structures
//--------------------------------------------------------------

fftw_init_threads();
fftw_plan_with_nthreads(n_nested);
//fftw_mpi_init();
myfft Fmy( N_f,N_f_thz, n_thread,n_nested);
Fmy.assign();


std::vector<int> iterations(n_thread,0);



for(int N_sta=0;N_sta<Num_stage;N_sta++)
{

        mydata.Energy_ir.setZero();
        mydata.Energy_thz.setZero();
        U_thz.setZero();
        U_thz.array()+=delta_safe;
        N_zz=my_z.N_z*l_sta[N_sta];

        for(int j=0;j<N_zz-1;j++)//N_zz-1
        {  
            myconst.chi_2=my_z.chi_2_z(j);
            myconst.phase(my_z.z_0(j),stage_phase,my_z.dz,LiNb);
            //-----------------------------------------------------
            //numerical method
            //------------------------------------------------------
            RK.method(world_size,n_thread, Fmy, myconst,  U_ir, P_ir,  U_thz, P_thz,  LiNb ,my_z.dz,iterations );

            //CN.method)();
            //------------------------------------------------------
            if(world_rank==0)
            {   mydata.spectrum_energy(df,U_ir,U_thz,j,my_r,LiNb);
                print(j);  
// 		if((j%10)==0){
// 			save_array(tf.N_f,my_r.N_r,U_ir,save_path+"E_ir_movie_"+std::to_string(j)+".txt" );
// 			save_array(tf.N_f_thz,my_r.N_r,U_thz,save_path+"E_thz_movie_"+std::to_string(j)+".txt" );
// 		}
            }

        }

    
        if(world_rank==0)
        { 
            mydata.save(N_sta,my_z,Fmy,U_ir,U_thz,tf,myconst.phase_thz,my_r,LiNb);
            print(iterations);
            stage_phase+=my_z.z_0(N_zz-1);
        }
        
        
        
}



//****************************************************
//REMOVE all the input pointers from fftw, VERY IMPORTANT
//****************************************************
Fmy.my_fftw_clean();





time_count.end(world_rank);
MPI_Finalize();

return 0;
    
}

