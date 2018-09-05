#ifndef MY_FFTW_H
#define MY_FFTW_H
//**************************************************
//fftw initialization
// %matlab fft from time to frequency is int exp(-iwt) dt
//************************************************
//fftw normalisation , overall operation(fft and ifft back) should be normalized by 1/N_f   (df*dt);
//fftw spectrum range [0 fs/2 -fs/2 0]
//ifft exp(iwt) from frequeny to time normalized by *df
//fft  exp(-iwt) from time to frequncy normalized by *dt
//*******************************************************************

struct myfft{

fftw_plan fft=nullptr;
fftw_plan ifft=nullptr;
fftw_plan r2c=nullptr;
fftw_plan c2r=nullptr;
int N_f;
int N_f_thz;
int n_thread;
int n_nested;
fftw_complex *in,*out,*E_t_ir,*E_ir,*E_thz;
double *E_t_thz, *I_t_ir;

myfft(const int& N_f_,const int&N_f_thz_,const int& n_thread_,const int& n_nested_): N_f(N_f_),N_f_thz(N_f_thz_),n_thread(n_thread_),n_nested(n_nested_){}
~myfft(){};
void my_fftw_clean();
void assign();
}; 


void myfft::assign(){
E_t_ir=fftw_alloc_complex(N_f*n_thread);
out= fftw_alloc_complex(N_f*n_thread);
in=fftw_alloc_complex(N_f*n_thread);
E_ir=fftw_alloc_complex(N_f*n_thread);
E_thz=fftw_alloc_complex((N_f_thz+1)*n_thread);
//define real array; E_thz negative frequency part is complex conjugate of itself, so E_t_thz always real;
E_t_thz=(double*)fftw_malloc(sizeof(double) * N_f*n_thread);
I_t_ir=(double*)fftw_malloc(sizeof(double) * N_f*n_thread);

ifft=fftw_plan_dft_1d(N_f,fftw_alloc_complex(N_f),fftw_alloc_complex(N_f),1, FFTW_MEASURE);
c2r= fftw_plan_dft_c2r_1d(N_f,fftw_alloc_complex(N_f_thz+1),fftw_alloc_real(N_f),FFTW_MEASURE);
fft=fftw_plan_dft_1d(N_f,fftw_alloc_complex(N_f),fftw_alloc_complex(N_f),-1,FFTW_MEASURE);
r2c=fftw_plan_dft_r2c_1d(N_f,fftw_alloc_real(N_f),fftw_alloc_complex(N_f_thz+1),FFTW_MEASURE);
}


void myfft::my_fftw_clean()
{
fftw_cleanup_threads();
//fftw_mpi_cleanup();                                                           
fftw_destroy_plan(fft);
fftw_destroy_plan(ifft);
fftw_destroy_plan(r2c);
fftw_destroy_plan(c2r);
fftw_free(E_t_thz);
fftw_free(out);
fftw_free(E_ir);
fftw_free(E_thz);
fftw_free(in);
fftw_free(E_t_ir);
fftw_free(I_t_ir);
fftw_cleanup();
}




#endif
