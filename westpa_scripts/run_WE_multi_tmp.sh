#!/bin/bash
#SBATCH --job-name=gamd_run
#SBATCH --account=ahnlab
#SBATCH --partition=gpu-ahn
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --gres=gpu:2
#SBATCH --time=120:00:00

##############################################################################
# 0) Setup Local Storage & Environment
##############################################################################
set -x
cd "$SLURM_SUBMIT_DIR"

source ~/.bashrc
module load conda3/4.X cuda/11.8.0 amber/22
source activate openmm_env

# Set up TMPDIR for job
TMPDIR="/tmp/${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
cp -r "$SLURM_SUBMIT_DIR"/* "$TMPDIR"
cd "$TMPDIR"

export WEST_SIM_ROOT="$TMPDIR"
export SERVER_INFO="$TMPDIR/west_zmq_info.json"

# Store logs and checkpoints in RAM
export TMPDIR_RAM="/dev/shm/gamd_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR_RAM"

##############################################################################
# 1) Configure MPS and set workers
##############################################################################
export WORKERS_PER_GPU=16
export MPS_PERCENTAGE=100
export CUDA_MPS_ACTIVE_THREAD_PERCENTAGE=60

IFS=',' read -ra DEVICES <<< "$CUDA_VISIBLE_DEVICES"

for gpuid in "${DEVICES[@]}"; do
   export CUDA_VISIBLE_DEVICES=$gpuid
   export CUDA_MPS_PIPE_DIRECTORY="/tmp/nvidia-mps-$SLURM_JOB_ID-$gpuid"
   export CUDA_MPS_LOG_DIRECTORY="/tmp/nvidia-log-$SLURM_JOB_ID-$gpuid"
   mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
   nvidia-cuda-mps-control -d
done

export CUDA_VISIBLE_DEVICES="$(IFS=','; echo "${DEVICES[*]}")"

##############################################################################
# 2) Initialize the WESTPA simulation
##############################################################################
./init.sh || { echo "Error: init.sh failed"; exit 1; }
source env.sh || exit 1

export WM_ZMQ_MASTER_HEARTBEAT=100
export WM_ZMQ_WORKER_HEARTBEAT=100
export WM_ZMQ_TIMEOUT_FACTOR=300

##############################################################################
# 3) Start w_run (ZMQ master) in the background
##############################################################################
w_run --work-manager=zmq \
      --n-workers=0 \
      --zmq-mode=master \
      --zmq-write-host-info="$SERVER_INFO" \
      --zmq-comm-mode=tcp \
      &> "$TMPDIR/west-${SLURM_JOB_ID}-master.log" &

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
# 4) Launch GPU Workers from Local Storage
##############################################################################
total_workers=$(( ${#DEVICES[@]} * WORKERS_PER_GPU ))

for gpuid in "${DEVICES[@]}"; do
    for ((w=1; w<=WORKERS_PER_GPU; w++)); do
        (
          export CUDA_VISIBLE_DEVICES="$gpuid"
          
          w_run --work-manager=zmq \
                --n-workers=1 \
                --zmq-mode=client \
                --zmq-read-host-info="$SERVER_INFO" \
                --zmq-comm-mode=tcp \
                &> "$TMPDIR/west-${SLURM_JOB_ID}-worker-gpu${gpuid}-w${w}.log"
        ) &
    done
done

wait

##############################################################################
# 5) Clean up MPS and Temporary Files
##############################################################################
for gpuid in "${DEVICES[@]}"; do
   export CUDA_VISIBLE_DEVICES=$gpuid
   export CUDA_MPS_PIPE_DIRECTORY="/tmp/nvidia-mps-$SLURM_JOB_ID-$gpuid"
   echo "Stopping MPS for GPU $gpuid"
   echo quit | nvidia-cuda-mps-control
done

rm -rf "$TMPDIR"

echo "Run complete."
exit 0
