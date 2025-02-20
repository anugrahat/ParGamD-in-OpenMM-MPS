#!/bin/bash
#SBATCH --job-name="chignolin_run"
#SBATCH --output="job.out"
#SBATCH -p GPU-shared
#SBATCH -N 1
#SBATCH --gpus=h100-80:1
#SBATCH -t 04:00:00



module load gcc/13.3.1-p20240614
module load openmpi/5.0.3-gcc13.2.1

module load cuda/12.4.0
module load anaconda3/2024.10-1

#module load amber/22
conda activate westpa


# Navigate to the directory containing the input files
#cd /home/anugraha/pargamd_openmm/ParGaMD/common_files/       # Run your simulation

python gamdRunner -p CUDA  -r  xml upper_dual_temp.xml 
