
void my_mpi_initial(int argc,char **argv,int& world_size,int& world_rank)
{    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    // Get the number of processes

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
 
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);


    printf("Hello world from processor %s, rank %d"
           " out of %d processors\n",
           processor_name, world_rank, world_size);
}

//detact wether the colum can be distributed to all processors
//int my_mpi_test(int& world_size,int& world_rank,int& n_thread,const int& N_x)
//{
//	if(N_x % (world_size*n_thread) != 0)
//	    {
//	        MPI_Finalize();
//	        if(world_rank == 0)
//	            std::cout << "**** N_x="<<N_x<<" is not divisible by " << world_size*n_thread << 				" ...quitting..."<< std::endl;
//        	return 1;
// 	   }
//return 0;
//}

void adjust_mesh_r(int& N_r, int& world_size, int& n_thread, double& crystal_size, double& dr)
{

if((world_size*n_thread)>int (crystal_size/dr))
	{
   	 dr=crystal_size/(world_size*n_thread);
   	 N_r=world_size*n_thread;
	}
else
	{   
	 N_r=(crystal_size)/dr-(int) ((crystal_size)/dr)%(world_size*n_thread);
	 dr=crystal_size/N_r;
	}

}
