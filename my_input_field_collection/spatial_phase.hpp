double spatial_phase(const MESH::spacial_r& my_r,const int i){
const double pi = M_PI;
//-------------------------------------- 
//random phase
//--------------------------------------
    
   //double my_p=2.0*pi*(rand() % 10 + 1)/10.0;  
   
//-------------------------------------- 
//sin phase
//--------------------------------------
double my_p=pi*cos(2.0*pi*4.0*i*my_r.dr/my_r.r_0(my_r.N_r));
//-------------------------------------- 
//linear phase
//--------------------------------------
//double my_p=2.0*pi*4.0*i*my_r.dr/my_r.r_0(my_r.N_r);

return my_p;
}
