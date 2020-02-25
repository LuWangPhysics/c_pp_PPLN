//all my print functions
#include<string>
#include <typeinfo>
#include <sstream>

#ifndef MY_PRINT_H
#define MY_PRINT_H


void print(const std::vector<int>& iterations)
{
    int index = 0;
    for (const auto i : iterations)
    {
        std::cout << "Thread " << index++ << " -> " << i << " iterations" << std::endl;
    }
}
template<typename T>
void print(const T& x)
{std::cout<<x<<std::endl;}

void test_line(int n){
     std::cout<<"****************************************************************"<<n<<std::endl;
}




//---------------------------------------------------------------------------------------------------------
//my save functions
//---------------------------------------------------------------------------------------------------------
        
      


        void save_array(int N_f,const double* a,std::string b )
        {
            
            std::ofstream myfile(b);
            myfile.precision(15);
            myfile.setf(std::ios::showpoint);
	    myfile.setf(std::ios::scientific);
            for (int i=0;i<N_f;i++)myfile <<a[i]<<"\n "; 
            myfile.close(); 
      

        }

        void save_array(int N_f,fftw_complex* a,std::string b )
        {
       
          
            std::ofstream myfile(b);
            myfile.precision(15);
            myfile.setf(std::ios::scientific);
            myfile.setf(std::ios::showpoint);
            for (int i=0;i<N_f;i++)myfile <<a[i][0]<<","<<a[i][1]<<"\n "; 

            myfile.close(); 
            
        }

        void save_array(int N_f,int N_x,const Eigen::MatrixXcd& a,std::string b )
        { 
        
            std::ofstream myfile(b);
            myfile.precision(15);
            myfile.setf(std::ios::scientific);
            myfile.setf(std::ios::showpoint);

                for (int i=0;i<N_f;i++)
                {
                    for (int j=0;j<N_x;j++)
                    {
                        myfile<< a(i,j).real()<<","<<a(i,j).imag() << ","; 
                    }
                  myfile << "\n";
                }
            myfile.close();
            
        }

#endif
