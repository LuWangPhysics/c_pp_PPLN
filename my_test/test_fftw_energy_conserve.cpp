
#include <iostream>
#include <sstream>
#include <fstream>//read and write file
#include <math.h> // included to use the "fabs" function
#include <fftw3.h>

int main(void)
{



int size =8;
double dt=1;
double df=1/(size*dt);
fftw_complex *out_cpx,*out_compare,*array_tc;
double  *out,*array_t;
fftw_plan fft;
fftw_plan fft_compare;
fftw_plan ifft;
out_cpx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(size/2+1));
out_compare = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(size));
array_tc = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(size));
array_t = (double*)fftw_malloc(sizeof(double) *size);
out=(double*)fftw_malloc(sizeof(double) * size);


for(int i=0;i<size;i++)
{
array_t[i]=i*1.0;
array_tc[i][0]=i*1.0;
array_tc[i][1]=0.0;
}

fft = fftw_plan_dft_r2c_1d(size, array_t, out_cpx, FFTW_ESTIMATE);  //Setup fftw plan for fft
fftw_execute(fft);


ifft = fftw_plan_dft_c2r_1d(size, out_cpx, out, FFTW_ESTIMATE);   //Setup fftw plan for ifft
fftw_execute(ifft);
fft_compare=fftw_plan_dft_1d(size, array_tc, out_compare, 1,FFTW_ESTIMATE);
fftw_execute(fft_compare);






double energy=0;
for(int i=0;i<size;i++)
{
std::cout<<array_t[i]<<std::endl;
energy+=(pow(array_t[i],2))*dt;
}
std::cout<<energy<<std::endl;
std::cout<<"******************************************"<<std::endl;
energy=0;

for(int i=0;i<(size/2+1);i++)
{
std::cout<<out_cpx[i][0]<<","<<out_cpx[i][1]<<std::endl;
energy+=2*(pow(out_cpx[i][0],2)+pow(out_cpx[i][1],2))*pow(dt,2)*df;
}
energy-=((pow(out_cpx[0][0],2)+pow(out_cpx[0][1],2))+(pow(out_cpx[size/2][0],2)+pow(out_cpx[size/2][1],2)))*pow(dt,2)*df;
std::cout<<energy<<std::endl;
std::cout<<"******************************************"<<std::endl;


energy=0;
for(int i=0;i<(size);i++)
{
std::cout<<out_compare[i][0]<<","<<out_compare[i][1]<<std::endl;
energy+=(pow(out_compare[i][0],2)+pow(out_compare[i][1],2))*pow(dt,2)*df;
}
std::cout<<energy<<std::endl;
std::cout<<"******************************************"<<std::endl;
energy=0;
for(int i=0;i<size;i++)
{
std::cout<<out[i]/size<<std::endl;
energy+=(pow(out[i],2))*dt/size/size;
}
std::cout<<energy<<std::endl;
std::cout<<"******************************************"<<std::endl;


fftw_destroy_plan(fft);
fftw_destroy_plan(ifft);
fftw_destroy_plan(fft_compare);
fftw_free(out_cpx);
fftw_free(array_t);
fftw_free(array_tc);
fftw_free(out_compare);
fftw_free(out);
return 0;
}


// //test of E_ir to E_t_ir. correct
// if(world_rank==0){              
// double test_1=0;
// for(int i=0;i<N_f;i++)
// test_1+=(pow(U_ir[i].real(),2)+pow(U_ir[i].imag(),2));
// std::cout<<test_1*myconst.df<<std::endl;
// 
//  test_1=0;
// for(int i=0;i<N_f;i++)
// test_1+=(pow(Fmy.E_t_ir[i][0],2)+pow(Fmy.E_t_ir[i][1],2))*pow(myconst.df,2);
// std::cout<<test_1/(N_f*myconst.df)<<std::endl;}


// //test of E_thz to E_t_thz. correct c2r;
// if(world_rank==0){              
// double test_1=0;
// for(int i=0;i<N_f_thz;i++)
// test_1+=(pow(U_thz[i].real(),2)+pow(U_thz[i].imag(),2));
// std::cout<<2*test_1*myconst.df<<std::endl;
// 
//  test_1=0;
// for(int i=0;i<N_f;i++)
// test_1+=(pow(Fmy.E_t_thz[i],2))*pow(myconst.df,2);
// std::cout<<test_1/(N_f*myconst.df)<<std::endl;}
// }
//test of r2c correct if remove the first and last data point once.
// double test_1=0;
// for(int i=0;i<N_f;i++)
// test_1+=(pow(Fmy.I_t_ir[i],2));
// std::cout<<test_1*dt<<std::endl;
// fftw_execute_dft_r2c(Fmy.r2c,&Fmy.I_t_ir[0],&Fmy.out[0]);
// 
//  test_1=0;
// for(int i=0;i<N_f_thz+1;i++)
// test_1+=(pow(Fmy.out[i][0],2)+pow(Fmy.out[i][1],2))*pow(1/(N_f*df),2)*2;
// 
// test_1-=(pow(Fmy.out[0][0],2)+pow(Fmy.out[0][1],2))*pow(1/(N_f*df),2)+(pow(Fmy.out[N_f_thz][0],2)+pow(Fmy.out[N_f_thz][1],2))*pow(1/(N_f*df),2);
// std::cout<<test_1*df<<std::endl;
//  
// save_array(N_f_thz,1,1,Fmy.out,"I_f_test.txt");
// 
// fftw_execute_dft_c2r(Fmy.c2r,&Fmy.out[0],&Fmy.I_t_ir[0]);
// test_1=0;
// for(int i=0;i<N_f;i++)
// test_1+=(pow(Fmy.I_t_ir[i]/N_f,2));
// std::cout<<test_1*dt<<std::endl;






/*
double test_1=0;
test_1=0;
for(int i=0;i<N_f_thz;i++)
test_1+=(pow(U_thz(i,0).real(),2)+pow(U_thz(i,0).imag(),2));
std::cout<<test_1*df*2<<std::endl;
         
                for(int i=0;i<N_f_thz;i++){
                    Fmy.E_thz[1+i][0]=(U_thz(i,0)*conj(myconst.phase_thz(i))).real();
                    Fmy.E_thz[1+i][1]=(U_thz(i,0)*conj(myconst.phase_thz(i))).imag();
                }
                Fmy.E_thz[0][0]=0;
                Fmy.E_thz[0][1]=0;
                save_array(N_f_thz,1,1,Fmy.E_thz,"E_thz.txt");
                test_1=0;
for(int i=0;i<N_f_thz;i++)
test_1+=(pow(Fmy.E_thz[i][0],2)+pow(Fmy.E_thz[i][1],2));
std::cout<<test_1*df*2<<std::endl;
fftw_execute_dft_c2r(Fmy.c2r,&Fmy.E_thz[0],&Fmy.E_t_thz[0]);


test_1=0;
for(int i=0;i<N_f;i++)
test_1+=(pow(Fmy.E_t_thz[i],2))*df*df;
std::cout<<test_1*dt<<std::endl;
fftw_execute_dft_r2c(Fmy.r2c,&Fmy.E_t_thz[0],&Fmy.out[0]);

 test_1=0;
for(int i=0;i<N_f_thz+1;i++)
test_1+=(pow(Fmy.out[i][0],2)+pow(Fmy.out[i][1],2))*dt*dt*2*df*df;

std::cout<<test_1*df<<std::endl;





save_array(N_f_thz,1,2,U_thz,"u_thz.txt");
U_thz.col(0)=conj(myconst.phase_thz.array());
save_array(N_f_thz,1,2,U_thz,"p_thz.txt");*/
