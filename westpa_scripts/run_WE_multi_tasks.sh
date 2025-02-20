#!/bin/bash
#SBATCH --job-name=multiGPU_mps_test
#SBATCH --account=ahnlab
#SBATCH --partition=gpu-ahn
#SBATCH --nodes=1
#SBATCH --ntasks=16             # e.g. 16 total tasks
#SBATCH --cpus-per-task=4       # 4 CPU cores per task
#SBATCH --gres=gpu:2            # request 2 GPUs on a single node
#SBATCH --time=48:00:00

set -x  # Debugging

#############################################
# 0. Load Modules & Activate Environment
#############################################
cd "$SLURM_SUBMIT_DIR"
source ~/.bashrc

module load conda3/4.X cuda/11.8.0 amber/22  # Adjust as needed
source activate openmm_env

echo "Running on host: $(hostname)"
echo "SLURM_JOB_GPUS: $SLURM_JOB_GPUS"  # e.g. "0,1"
export CUDA_VISIBLE_DEVICES=$(echo "$SLURM_JOB_GPUS" | tr ',' ' ')
echo "CUDA_VISIBLE_DEVICES after assignment: $CUDA_VISIBLE_DEVICES"

IFS=' ' read -ra DEVICES <<< "$CUDA_VISIBLE_DEVICES"
NUM_GPUS=${#DEVICES[@]}
echo "Parsed GPU IDs: ${DEVICES[*]}"
echo "Number of GPUs: $NUM_GPUS"

# Optionally, confirm OpenMM sees the CUDA platform:
python -c "import openmm; print('OpenMM version:', openmm.__version__)"

#############################################
# 1. Set Up MPS on Each GPU
#############################################
# This will let multiple tasks share each GPU concurrently
echo quit | nvidia-cuda-mps-control 2>/dev/null || true
sleep 2

for gpuid in "${DEVICES[@]}"; do
    echo "Starting MPS on GPU $gpuid"

    export CUDA_VISIBLE_DEVICES=$gpuid
    export CUDA_MPS_PIPE_DIRECTORY=/tmp/nvidia-mps-$SLURM_JOB_ID-$gpuid
    export CUDA_MPS_LOG_DIRECTORY=/tmp/nvidia-log-$SLURM_JOB_ID-$gpuid

    rm -rf "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
    mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"

    # Start MPS daemon
    nvidia-cuda-mps-control -d &

    # Verify GPU
    sleep 5
    python -c "import openmm; p=openmm.Platform.getPlatformByName('CUDA'); print('MPS test on GPU=$gpuid =>', p.getName())" || {
        echo "Error: GPU=$gpuid not accessible."
        exit 1
    }
done

# Restore visibility to all GPUs (e.g. "0 1")
export CUDA_VISIBLE_DEVICES=$(IFS=' '; echo "${DEVICES[*]}")
echo "Final CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"

#############################################
# 2. Initialize the Master Process (ZMQ or MPI)
#############################################
# Example with ZeroMQ
./init.sh
source env.sh
SERVER_INFO="west_zmq_info.json"

# Launch master with zero local workers
w_run --work-manager=zmq --n-workers=0 --zmq-mode=master \
      --zmq-write-host-info="$SERVER_INFO" \
      --zmq-comm-mode=tcp \
      &> "master_${SLURM_JOB_ID}.log" &

sleep 5
[[ ! -f "$SERVER_INFO" ]] && { echo "Error: Master server info not found!"; exit 1; }

#############################################
# 3. Launch Multiple Tasks in Parallel
#############################################
# We'll do 16 tasks, each one picks GPU by:
#    gpu_assigned = (SLURM_PROCID % NUM_GPUS)
# Use one srun for full concurrency.
#############################################
srun --ntasks=8 --cpus-per-task=4 bash node.sh \
     --work-manager=zmq \
     --n-workers=1 \
     --zmq-mode=client \
     --zmq-read-host-info="$SERVER_INFO" \
     --zmq-comm-mode=tcp

wait

#############################################
# 4. Cleanup MPS
#############################################
for gpuid in "${DEVICES[@]}"; do
    echo "Stopping MPS on GPU $gpuid"
    export CUDA_VISIBLE_DEVICES=$gpuid
    echo quit | nvidia-cuda-mps-control
    rm -rf "/tmp/nvidia-mps-$SLURM_JOB_ID-$gpuid" "/tmp/nvidia-log-$SLURM_JOB_ID-$gpuid"
done

echo "Job completed successfully."
