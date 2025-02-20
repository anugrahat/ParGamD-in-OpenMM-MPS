#!/bin/bash
#SBATCH --job-name=gamd_run
#SBATCH --account=ahnlab
#SBATCH --partition=gpu-ahn
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --gres=gpu:2
#SBATCH --time=120:00:00

# GPU and CUDA settings (No MPS)
export WORKERS_PER_GPU=8   # Number of workers per GPU (adjust based on memory/performance)
export MPS_PERCENTAGE=100  # Unused when MPS is disabled, but left here to keep consistency
export CUDA_DEVICE_MAX_CONNECTIONS=32
export OPENMM_CPU_THREADS=32

set -x
cd "$SLURM_SUBMIT_DIR"
source ~/.bashrc
module load conda3/4.X cuda/11.8.0 amber/22
source activate openmm_env
export PATH=/home/anugraha/.conda/envs/openmm_env/bin:$PATH
export WEST_SIM_ROOT=$SLURM_SUBMIT_DIR

# Just check available GPUs for clarity
echo "Initially, CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
nvidia-smi -L

# Use the comma-separated list of GPUs from SLURM
IFS=',' read -ra DEVICES <<< "$CUDA_VISIBLE_DEVICES"

# Confirm which devices are being used
nvidia-smi

# WESTPA Initialization
./init.sh || { echo "Error: init.sh failed"; exit 1; }
source env.sh || exit 1
SERVER_INFO=$WEST_SIM_ROOT/west_zmq_info.json

# Add before w_run
export WM_ZMQ_MASTER_HEARTBEAT=100
export WM_ZMQ_WORKER_HEARTBEAT=100
export WM_ZMQ_TIMEOUT_FACTOR=300

w_run --work-manager=zmq --n-workers=0 --zmq-mode=master \
     --zmq-write-host-info="$SERVER_INFO" --zmq-comm-mode=tcp \
     &> "west-$SLURM_JOBID-local.log" &

for ((n=0; n<60; n++)); do
   if [ -e "$SERVER_INFO" ]; then
       cat "$SERVER_INFO"
       break
   fi
   sleep 1
done
[ ! -e "$SERVER_INFO" ] && { echo "Error: Server failed to start"; exit 1; }

# Worker Launch
# Calculate total workers based on the number of GPUs * workers per GPU
total_workers=$(( ${#DEVICES[@]} * WORKERS_PER_GPU ))

# Debug statements
echo "Debug: Total workers = $total_workers"
echo "Debug: DEVICES array = ${DEVICES[@]}"
echo "Debug: Node count = $SLURM_NNODES"

# First node, single task
srun --nodes=1 --ntasks=1 bash node.sh \
    "$SLURM_SUBMIT_DIR" \
    "$SLURM_JOBID" \
    "$(hostname)" \
    "$CUDA_VISIBLE_DEVICES" \
    --work-manager=zmq \
    --n-workers="$total_workers" \
    --zmq-mode=client \
    --zmq-read-host-info="$SERVER_INFO" \
    --zmq-comm-mode=tcp &

# Other nodes, if any
srun -N$SLURM_NNODES -n$SLURM_NNODES bash node.sh \
    "$SLURM_SUBMIT_DIR" \
    "$SLURM_JOBID" \
    "$(hostname)" \
    "$CUDA_VISIBLE_DEVICES" \
    --work-manager=zmq \
    --n-workers="$total_workers" \
    --zmq-mode=client \
    --zmq-read-host-info="$SERVER_INFO" \
    --zmq-comm-mode=tcp &
wait

# No MPS cleanup since MPS was never enabled
echo "Run complete. No MPS was used."
