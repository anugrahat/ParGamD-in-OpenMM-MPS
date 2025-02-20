#!/bin/bash
#SBATCH --job-name="chignolin_run"
#SBATCH --output="job.out"
#SBATCH -p GPU-shared
#SBATCH -N 1
#SBATCH --cpus-per-task=12
#SBATCH --gpus=h100-80:1
#SBATCH -t 04:00:00




##############################################################################
# 0) Basic environment and module loads
##############################################################################
set -x
cd "$SLURM_SUBMIT_DIR"

module load gcc/13.3.1-p20240614
module load openmpi/5.0.3-gcc13.2.1

module load cuda/12.4.0
module load anaconda3/2024.10-1

#module load amber/22
conda activate westpa

#export PATH="/home/$USER/.conda/envs/openmm_env/bin:$PATH"

export WEST_SIM_ROOT="$SLURM_SUBMIT_DIR"

cd $WEST_SIM_ROOT

##############################################################################
# 1) (Optional) Configure MPS and set workers
##############################################################################
# How many worker processes per GPU do you want?
export WORKERS_PER_GPU=8

# If you want to limit fraction of GPU SM usage, set MPS_PERCENTAGE
export MPS_PERCENTAGE=100
export CUDA_MPS_ACTIVE_THREAD_PERCENTAGE=60

# Tweak concurrency environment vars
#export CUDA_DEVICE_MAX_CONNECTIONS=2
export OPENMM_CPU_THREADS=2
export OMP_NUM_THREADS=1

# Show which GPUs Slurm gave us
echo "Initially, CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
nvidia-smi -L

# Split out the GPU IDs
IFS=',' read -ra DEVICES <<< "$CUDA_VISIBLE_DEVICES"

# (Optional) Start MPS for each GPU
for gpuid in "${DEVICES[@]}"; do
    export CUDA_MPS_PIPE_DIRECTORY="/tmp/nvidia-mps-$SLURM_JOB_ID-$gpuid"
    export CUDA_MPS_LOG_DIRECTORY="/tmp/nvidia-log-$SLURM_JOB_ID-$gpuid"
    mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
done
# Restore the full list of GPU IDs
export CUDA_VISIBLE_DEVICES="$(IFS=','; echo "${DEVICES[*]}")"
#nvidia-smi

nvidia-cuda-mps-control -d
sleep 5
echo set_active_thread_percentage 60 | nvidia-cuda-mps-control






##############################################################################
# 2) Initialize the WESTPA simulation
##############################################################################
# ./init.sh || { echo "Error: init.sh failed"; exit 1; }
source env.sh || exit 1

export SERVER_INFO="${WEST_SIM_ROOT}/west_zmq_info.json"

# These environment variables tweak the ZMQ heartbeat & timeouts
export WM_ZMQ_MASTER_HEARTBEAT=50
export WM_ZMQ_WORKER_HEARTBEAT=50
export WM_ZMQ_TIMEOUT_FACTOR=100

##############################################################################
# 3) Start w_run (ZMQ master) in the background
##############################################################################
w_run --work-manager=zmq \
      --n-workers=0 \
      --zmq-mode=master \
      --zmq-write-host-info="$SERVER_INFO" \
      --zmq-comm-mode=tcp \
      &> "west-${SLURM_JOB_ID}-master.log" &

# Wait up to 60s for the server info to appear
for ((n=0; n<60; n++)); do
   if [ -e "$SERVER_INFO" ]; then
       cat "$SERVER_INFO"
       break
   fi
   sleep 1
done
if [ ! -e "$SERVER_INFO" ]; then
   echo "ERROR: ZMQ master did not start properly."
   exit 1
fi

##############################################################################
# 4) Launch GPU workers
##############################################################################

# Start monitoring GPU usage in the background
#LOG_FILE="gpu_utilization_${SLURM_JOB_ID}.log"
#echo "Logging GPU utilization to $LOG_FILE"
#nvidia-smi --query-gpu=timestamp,index,name,utilization.gpu,utilization.memory,memory.total,memory.used --format=csv -l 10 > "$LOG_FILE" &
#NVIDIA_SMI_PID=$!


# We'll spawn (WORKERS_PER_GPU) processes for each GPU. This ensures you have
# multiple segments running concurrently on the same GPU via MPS.

total_workers=$(( ${#DEVICES[@]} * WORKERS_PER_GPU ))
echo "Debug: DEVICES = ${DEVICES[@]}"
echo "Debug: total workers = $total_workers"

# Because we are on a single node (nodes=1), we can simply background-run
# each worker and let them connect to the master. Each worker will set
# CUDA_VISIBLE_DEVICES to a single GPU. In total, we’ll have $total_workers
# processes. If your cluster requires `srun` for each process, you can adapt
# it accordingly.

for gpuid in "${DEVICES[@]}"; do
    for ((w=1; w<=WORKERS_PER_GPU; w++)); do

        # Launch in background
        (
          # Pin this worker to one GPU
          export CUDA_VISIBLE_DEVICES="$gpuid"

          # Optional: you might want to isolate CPU cores here as well,
          # e.g., via taskset or cgroups if you have many CPUs. 
          
          # Launch the worker in ZMQ client mode.
          w_run --work-manager=zmq \
                --n-workers=1\
                --zmq-mode=client \
                --zmq-read-host-info="$SERVER_INFO" \
                --zmq-comm-mode=tcp \
                &> "west-${SLURM_JOB_ID}-worker-gpu${gpuid}-w${w}.log"
        ) &
    done
done

# Wait for all background workers plus the master
wait

##############################################################################
# 5) Clean up MPS
##############################################################################
# Stop MPS only if no other jobs are using it
if [[ -z "$(pgrep w_run)" ]]; then
    for gpuid in "${DEVICES[@]}"; do
        export CUDA_VISIBLE_DEVICES=$gpuid
        export CUDA_MPS_PIPE_DIRECTORY="/tmp/nvidia-mps-$SLURM_JOB_ID-$gpuid"
        echo "Stopping MPS for GPU $gpuid"
    done
    echo quit | nvidia-cuda-mps-control
else
    echo "MPS is still in use by other jobs, not stopping."
fi

echo "Run complete."
exit 0
