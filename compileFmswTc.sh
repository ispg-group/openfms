#!/bin/bash -l
#SBATCH --job-name="CompileTC"
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --output=tc.%j.o
#SBATCH --error=tc.%j.e
#SBATCH --mail-type=REQUEUE
#SBATCH --partition=rooster
#SBATCH --exclude=gpu202,gpu203,gpu204,gpu205
#SBATCH --mem=200G   

nproc
curl -V

module purge
module load cuda/11.8.0-gcc-8.5.0
module load intel-oneapi-compilers/2023.1.0
module load intel-oneapi-mkl/2023.1.0
module load openmm/7.3.0_cpu
module load intel-oneapi-mpi/2021.9.0

#export INCLUDE=${INCLUDE}:/usr/include
#export LIBRARY_PATH=/gpfsnyu/spack/opt/spack/linux-rhel8-icelake/gcc-8.5.0/cuda-11.8.0-jkl6337z3egvvf6pnbsraqt4l6qq4j5p/lib64/stubs:${LIBRARY_PATH}
#export LD_LIBRARY_PATH=/gpfsnyu/spack/opt/spack/linux-rhel8-icelake/gcc-8.5.0/cuda-11.8.0-jkl6337z3egvvf6pnbsraqt4l6qq4j5p/lib64/stubs:${LD_LIBRARY_PATH}


#export LD_LIBRARY_PATH=/scratch/wjg3/openmm_dev/build/7.3.0/lib/plugins:${LD_LIBRARY_PATH}
#export LD_LIBRARY_PATH=/scratch/wjg3/openmm_dev/build/7.3.0/lib:${LD_LIBRARY_PATH}
#export LIBRARY_PATH=/scratch/wjg3/openmm_dev/build/7.3.0/lib/plugins:${LIBRARY_PATH}
#export LIBRARY_PATH=/scratch/wjg3/openmm_dev/build/7.3.0/lib:${LIBRARY_PATH}

#export LD_LIBRARY_PATH=/gpfsnyu/spack/opt/spack/linux-rhel8-icelake/gcc-8.5.0/intel-oneapi-compilers-2023.1.0-zmeb7pybskx6cbhuv6rdjvzjcoawgedx/compiler/2023.1.0/linux/bin/intel64:${LD_LIBRARY_PATH}

export LD_LIBRARY_PATH=/gpfsnyu/packages/glovergroup/lib/simple-dftd3/1.2.1/lib64:${LD_LIBRARY_PATH}
export LIBRARY_PATH=/gpfsnyu/packages/glovergroup/lib/simple-dftd3/1.2.1/lib64:${LIBRARY_PATH}
export INCLUDE=/gpfsnyu/packages/glovergroup/lib/simple-dftd3/1.2.1/include:${INCLUDE}

echo "compiling on $HOSTNAME"

./configure --prefix build_2026_3_17_restraints_xyz_mpi
make -j 32
