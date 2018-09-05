.SUFFIXES: .c .cpp .o .ex

# the compiler: gcc for C program, define as g++ for C++
    # compiler flags:
  #  -g    adds debugging information to the executable file
  #  -Wall turns on most, but not all, compiler warnings
  CC =mpic++ -ggdb -Wall  -std=c++11 -fopenmp  -Wfatal-errors
  #-fopenmp 
  #for compile without mpi, change mpic++ to g++, mpirun
  # g++ -g -Wall 


  # to include external libraries to your code


INCLUDE= -I /home/luwang/cpp_ppln/c++_ppln_3d/
INCLUDE+=  -I /home/luwang/cpp_ppln/c++_ppln_3d/boost_1_65_0/
INCLUDE+=-I /home/luwang/cpp_ppln/c++_ppln_3d/eigen/ 

LIBS = -L /home/luwang/cpp_ppln/c++_ppln_3d/fftw-3.3.4/ -lfftw3_threads -lfftw3 -lm  
LIBS+=  -L/usr/lib64/mpich-3.2/bin/mpicc
LIBS+=-L /home/luwang/cpp_ppln/c++_ppln_3d/boost_1_65_0/stage/lib


#
 






.cpp.o:
	$(CC) $(INCLUDE) -c $< 

.c.o:
	gcc $(INCLUDE) -c $< 

.o.ex:
	@echo g++ ... -o $@ $< ... $(OBJS) ... $(LIBS) 
	@$(CC) -o $@ $< $(OBJS) $(LIBS)

######################################################################
# targets
######################################################################
OBJS  = 

objs: $(OBJS)

clean:
	rm -f *.o *.ex *~
