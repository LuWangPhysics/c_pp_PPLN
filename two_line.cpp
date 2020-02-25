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

    
int energy_scale= 100*strtol(argv[1], NULL, 10);                                       //fluence scan    


srand(time(NULL));                                                                     //seed the random number
const double pi = M_PI;   
const double c=3e8;                                                                    //Speed of light  
 //const double eps =8.85e-12;                                                         //Free space permittivity  


    
mytime time_count;
time_count.start();
std::cout.precision(10);                                                               //change the default prcision std   
    
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

double sigma_ch=(1.5/3.0);
int n_gaussian=1;                                                                      //super gaussian order
double sigma_in=sigma_ch*3e-3;                                                         // Input beam waist
//10/0.02a 3b 5c 7d 11e 13f, where e+f is either 0 or 1
double lambda_0    = 1030.4e-9;
double tau_fwhm    = 250e-12;                                                          //Input pulse fwhm
double tau = tau_fwhm/(sqrt(2*log(2)));                                                //Converting from FWHM to 1/e^2 value. e^(-t^2/(2*tau^2))
double f_max       = 45e12;    //4stages 120                                           //2*f_max =frequency range ,120
double df          = 0.0003e12;//0.0002e12 500ps 0.0004e12 150ps;                      //frequency step size
double omega_0 =  2*pi*c/lambda_0;
double f_c_thz=0.53e12;                                                                //define thz center frequency 


//peak fluence is fixed to damage and the energy is scaled with respect to peak fluence.

double energy=(energy_scale/1.1180e+04)*2*pi*pow(sigma_in,2)*std::tgamma(1.0*(n_gaussian+1.0)/n_gaussian)*1e5*sqrt(tau_fwhm/2/1e-8)/pow(2,(1.0/n_gaussian));   //pulse energy J
double F_peak=pow(2,(1.0/n_gaussian))*energy/(2*pi*pow(sigma_in,2)*std::tgamma(1.0*(n_gaussian+1.0)/n_gaussian)); //peak fluence
double energy_comb=100e-3;
double sigma_comb=sigma_in;
double F_peak_comb=pow(2,(1.0/n_gaussian))*energy_comb/(2*pi*pow(sigma_comb,2)*std::tgamma(1.0*(n_gaussian+1.0)/n_gaussian)); //peak fluence

//-----------------------------------------------------------------------
//spacial discritization
//---------------------------------------------------------------------
double dr =sigma_ch*8e-5;  //1.0e-4 2e-5 for phase study
double crystal_size=sigma_ch*10e-3; //5e-3 for flattop
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
//---------------------------------------------------------------------
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

int Nppln=400; //145 for 150ps 190for 500ps stages                                   // Total number of poling periods                                                     
int nz_ppln=75; //75                                                                 //number of steps in one period
double dppln=c/(2*LiNb.f_c_thz*(LiNb.n_thz_c-LiNb.n_g));                             // Thickness of each poling period
int   Num_stage=1;
MESH::spacial_z my_z(dppln,Nppln,nz_ppln,chi_2);


std::vector<double> l_sta={1};                                                        //smaller than 1, crystal length 
double stage_phase=0;
int   N_zz;

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
U_ir*=sqrt(1-pow(1-LiNb.n_ir_0,2)/pow(1+LiNb.n_ir_0,2));                             //fresnel loss

//---------------------------------------------------
//define calculation constants
//-----------------------------------------------------
P_const myconst(tf,LiNb, n2);

myconst.df=df;
myconst.n_nested=n_nested;

Eigen::MatrixXcd P_ir=Eigen::MatrixXcd::Zero(N_f,N_r);
Eigen::MatrixXcd P_thz=Eigen::MatrixXcd::Zero(N_f_thz,N_r);

std::string save_path="/gpfs/cfel/ufox/scratch/user/luwang/2D_250ps_2lines/bc_dieletric/";
std::string extra_name=std::to_string(energy_scale)+"mjpercm2_";
//std::string extra_name+="err_stepz" +std::to_string(nz_ppln)+"_dz_"+std::to_string(my_z.dz*1e6);
save_path+=extra_name;


Output mydata(save_path,tf,my_r,my_z,U_ir);



//--------------------------------------------------------------
//numerical method
//---------------------------------------------------------
//**************************************************
//RK method initialization//lowest storage rk method
//**************************************************
lowest_storage_RK RK( N_f_thz,N_f, my_r);
RK.assign();


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
        mydata.Energy_thz=0;
        U_thz.setZero();
        U_thz.array()+=delta_safe;
        N_zz=my_z.N_z*l_sta[N_sta];

        for(int j=0;j<N_zz-1;j++){  
            myconst.chi_2=my_z.chi_2_z(j);
            myconst.phase(my_z.z_0(j),stage_phase,my_z.dz,LiNb);
            //-----------------------------------------------------
            //numerical method
            //------------------------------------------------------
            RK.method(world_size,n_thread, Fmy, myconst,  U_ir, P_ir,  U_thz, P_thz,  LiNb ,my_z.dz,iterations );
           
            mydata.spectrum_energy(U_ir,U_thz,j,my_r,my_z,LiNb);
            
            print(j); 
            
            if(world_rank==0){  
                if((j%int(10e-3/my_z.dz))==0){
                  
                       save_array(my_z.N_z-1,mydata.Eff.data(),save_path+"eff.txt" );
                       save_array(tf.N_f,mydata.ir_spectrum.data(),save_path+"ir_spectrum_"+std::to_string(int(j*my_z.dz*1000))+".txt" );
                       
                }
            }
            

        }

    
        if(world_rank==0){   
            mydata.save(N_sta,my_z,Fmy,U_ir,U_thz,tf,myconst.phase_thz,my_r,LiNb,int(my_z.z_0[N_zz-1]*1000));
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

