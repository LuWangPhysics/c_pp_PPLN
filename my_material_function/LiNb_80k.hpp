//this code only valid for 80k
 void LiNb_80k_thz(const Eigen::ArrayXd& omega_thz,Eigen::ArrayXd &n_thz,Eigen::ArrayXd &alpha_thz){
double T=80;
Eigen::ArrayXd f=omega_thz/(2*M_PI)/3e10; //Converting to /cm units

double a_c=2.706*exp(-8.344e-4*T)+1.993*exp(0.001173*T);
double b_c=-1.613e-13*pow(T,3)+9.679e-11*pow(T,2)+2.254e-9*T+3.197e-5;
double c_c=3e-10;
int N=omega_thz.size();
n_thz=1.035*(a_c+b_c*f.pow(2) + c_c*f.pow(4));

//absorption coefficient fitting for 80k only valid for 1.5THz
double   p1 =-413.8 ;  
double   p2 = 260.2;
double   p3 = -414.2;
double   p4=-25.02;
Eigen::ArrayXd f_thz=omega_thz/(2*M_PI*1e12);
//alpha is negative value in unite /m;
alpha_thz=p1*f_thz.pow(3)+p2*f_thz.pow(2)+p3*f_thz.pow(1)+p4; 
int N_1=(f_thz<6).count();
for(int j=N_1;j<N;j++){
    alpha_thz(j)=alpha_thz(N_1);
}
}
 
 

void LiNb_80k_ir(const Eigen::ArrayXd& omega_ir,Eigen::ArrayXd &n_ir_b){
int N=omega_ir.size();
double lambda;

double T=80;
double k = (T - 273 - 24.5) / ( T - 273 + 570.82);
for(int i=0;i<N;i++){
   lambda =(2*M_PI)*3e8/1e-6/omega_ir[i];
 n_ir_b[i] = sqrt( 4.5820+( 0.09921 + 5.2716e-8 *k)/(pow(lambda,2) -pow((0.21090-4.9143e-8*k ),2))+2.2971e-8 * k- 0.021940*pow(lambda,2));   
}
}