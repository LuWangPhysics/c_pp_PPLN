double Gro_v(const Eigen::ArrayXd& omega_ir,const Eigen::ArrayXd& n_p){
// //calculate group refractvie index;
double n_g;

int N=omega_ir.size()/2;
n_g=(omega_ir[N]*n_p[N]-omega_ir[N-1]*n_p[N-1])/(omega_ir[N]-omega_ir[N-1]);

return n_g;

}