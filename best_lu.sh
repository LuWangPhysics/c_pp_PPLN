#!/bin/bash
#SBATCH --partition=cfel         					## this is the default partition.
#SBATCH -t 90:10:00              					## default is 1h. The maximum is partition dependent, have a look at sview or scontrol for details.
#SBATCH --nodes=1               					 ## number of nodes 
##SBATCH -n 50                   					 ## Number of threads. 
#SBATCH --output    my_batch_out/%j-%N.out                                # File to which STDOUT will be written
#SBATCH --error     my_batch_out/%j-%N.err                                # File to which STDERR will be written
#SBATCH --mail-type END                                                   # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user lu.wang@desy.de                                       # Email to which notifications will be sennet 

. /etc/profile.d/modules.sh
module load mpi/mpich-3.2-x86_64

##echo "test cascaded line phase front"

mpirun -np 1 ./two_line.ex
