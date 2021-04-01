#!/bin/sh
#SBATCH -p hernquist
#SBATCH -J b5.256
#SBATCH -n 256
#SBATCH --ntasks-per-node=64
#SBATCH -o ../output/OUTPUT.lsf
#SBATCH -e ../output/ERROR.lsf
#SBATCH --exclusive
#SBATCH -t 8640 # 6d #5760 # 4 day in min
#SBATCH --mail-user=philip.mocz@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem-per-cpu=3900

module purge
module load intel/17.0.4-fasrc01 openmpi/1.10.4-fasrc01
module load fftw/3.3.7-fasrc02
module load zlib/1.2.8-fasrc07  
module load szip/2.1-fasrc02  
module load hdf5/1.10.1-fasrc03 
module load gsl/2.4-fasrc01  
module load gmp/6.1.2-fasrc01


#srun -n $SLURM_NTASKS --mpi=pmi2 ../arepo/./Arepo param256_B0.25_M10.txt 1 > ../output/OUTPUT.$SLURM_JOB_ID 2> ERROR
#srun -n $SLURM_NTASKS --mpi=pmi2 ../arepo/./Arepo param256_B0.5_M10.txt 1 > ../output/OUTPUT.$SLURM_JOB_ID 2> ERROR
srun -n $SLURM_NTASKS --mpi=pmi2 ../arepo/./Arepo param256_B1_M10.txt 1 > ../output/OUTPUT.$SLURM_JOB_ID 2> ERROR
#srun -n $SLURM_NTASKS --mpi=pmi2 ../arepo/./Arepo param256_B2_M10.txt 1 > ../output/OUTPUT.$SLURM_JOB_ID 2> ERROR
#srun -n $SLURM_NTASKS --mpi=pmi2 ../arepo/./Arepo param256_B4_M10.txt 1 > ../output/OUTPUT.$SLURM_JOB_ID 2> ERROR
