#!/bin/bash
#
# submission-script for MPI-jobs, example
#
# To run this example, please copy it into your scratch directory (/scratch/<your userID>)
# then run 'sbatch mpi_hello_world.sbatch' -- you can rename the file if you'd like, the
# file extension (eg: .sbatch) is meaningless.
# 
# the SBATCH lines are not comments, they are processed by sbatch as arguements 
#
# -p <queuename>; what queue are you submitting the job to? You must have submit-access to the queue
# to identify what queues that you may submit to, please run 'queues'. To see what other parameters you
# may use with the 'sbatch' command, please run 'man sbatch'
#
#SBATCH -p standard
#
# -n <cores>; number of cores required by your job
#
#SBATCH -A dubayhamblin
#
#SBATCH -N 1
#
# --ntasks-per-node; the number of cores per-node. You should pick a value that will get you the best performance,
# consider the possiblity that your code may run on a compute-node already running different codes started by 
# different users.
#
#SBATCH --ntasks-per-node=20
#
# --mem-per-cpu=<MB>; amount of memory your job will need, per-core. Your job will be killed, if it uses
# more memory per-core then you request.
#
#SBATCH --mem-per-cpu=500
#
#SBATCH --export=ALL
#
#SBATCH --export=variables
#
#SBATCH -t 7-00:00:00
#
# This part, below the SBATCH lines, is the "job", and what will be run on the compute-node
#

# Please run 'module avail' to see what software we have provided for you.
module load gcc/7.1.0 python/3.6.8 ffmpeg intel/18.0 intelmpi/18.0 cuda pgi openmpi

cd $SLURM_SUBMIT_DIR

echo `pwd`

mpiexec -np 20 /home/rh3bf/mylammps/build/lmp -i run.lmp
