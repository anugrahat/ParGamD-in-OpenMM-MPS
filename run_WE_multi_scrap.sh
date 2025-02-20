#!/bin/bash
#SBATCH --job-name=gamd_run
#SBATCH --account=ahnlab
#SBATCH --partition=gpu-ahn
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --gres=gpu:2
#SBATCH --time=120:00:00

set -x  
cd "$SLURM_SUBMIT_DIR"

source ~/.bashrc
module load conda3/4.X cuda/11.8.0 amber/22
source activate openmm_env

export WEST_SIM_ROOT="$SLURM_SUBMIT_DIR"
export DESIRED_WORKERS_PER_GPU=4
export PATH="/home/$USER/.conda/envs/openmm_env/bin:$PATH"

# (Optional) start MPS
IFS=',' read -ra DEVICES <<< "$CUDA_VISIBLE_DEVICES"
for gpuid in "${DEVICES[@]}"; do
    export CUDA_VISIBLE_DEVICES=$gpuid
    export CUDA_MPS_PIPE_DIRECTORY="/tmp/nvidia-mps-$SLURM_JOB_ID-$gpuid"
    export CUDA_MPS_LOG_DIRECTORY="/tmp/nvidia-log-$SLURM_JOB_ID-$gpuid"
    mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
    nvidia-cuda-mps-control -d
done

# Restore the full GPU list after MPS setup
export CUDA_VISIBLE_DEVICES="$(IFS=','; echo "${DEVICES[*]}")"

# (Optional) GPU monitor
nvidia-smi --query-gpu=timestamp,index,utilization.gpu,utilization.memory,memory.used,memory.total \
           --format=csv,nounits -l 10 > "gpu_usage_${SLURM_JOB_ID}.csv" &
GPUMON_PID=$!

# Initialize
./init.sh || { echo "Error: init.sh failed"; exit 1; }
source env.sh || exit 1

export SERVER_INFO="${WEST_SIM_ROOT}/west_zmq_info.json"

# Start w_run (master) in the background
w_run --work-manager=zmq \
      --n-workers=0 \
      --zmq-mode=master \
      --zmq-write-host-info="$SERVER_INFO" \
      --zmq-comm-mode=tcp \
      &> "west-${SLURM_JOB_ID}-master.log" &

# Wait for server info
for ((n=0; n<60; n++)); do
    [ -e "$SERVER_INFO" ] && cat "$SERVER_INFO" && break
    sleep 1
done
if [ ! -e "$SERVER_INFO" ]; then
    echo "ERROR: ZMQ master did not start properly."
    exit 1
fi

# ------------------------------
# **Here** is where you embed the worker loop:
# ------------------------------
IFS=',' read -ra GPUS <<< "$CUDA_VISIBLE_DEVICES"
worker_pids=()
for gpuid in "${GPUS[@]}"; do
    for (( i=0; i<DESIRED_WORKERS_PER_GPU; i++ )); do
        CUDA_VISIBLE_DEVICES="$gpuid" \
        w_launch_worker \
            --n-workers=4 \
            --worker-name="${SLURM_NODENAME}-GPU${gpuid}-WORKER${i}" \
            --zmq-mode=client \
            --zmq-read-host-info="${SERVER_INFO}" \
            --zmq-comm-mode=tcp \
            &> "west-${SLURM_NODENAME}-gpu${gpuid}-worker${i}.log" &
        worker_pids+=($!)
    done
done

wait "${worker_pids[@]}"
wait  # wait for w_run as well

# Cleanup
kill $GPUMON_PID 2>/dev/null
for gpuid in "${DEVICES[@]}"; do
    export CUDA_VISIBLE_DEVICES=$gpuid
    export CUDA_MPS_PIPE_DIRECTORY="/tmp/nvidia-mps-$SLURM_JOB_ID-$gpuid"
    echo "Stopping MPS for GPU $gpuid"
    echo quit | nvidia-cuda-mps-control
done

echo "Run complete."
exit 0
