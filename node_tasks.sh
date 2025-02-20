#!/bin/bash
set -x

# node.sh
# Usage: node.sh [any extra arguments to w_run?]

# 1) Identify Slurm environment
echo "node.sh: SLURM_PROCID=$SLURM_PROCID"
echo "node.sh: HOSTNAME=$(hostname)"

# Let's read the total number of GPUs from our environment or just assume 2
NUM_GPUS=2

# Simple round-robin: GPU = (task rank) mod (number of GPUs)
GPU_ID=$(( SLURM_PROCID % NUM_GPUS ))
echo "node.sh: Using GPU $GPU_ID"

# 2) Set CUDA_VISIBLE_DEVICES
export CUDA_VISIBLE_DEVICES=$GPU_ID
nvidia-smi

# 3) Actually run the worker
# Our main script passed ZMQ client arguments in $@
w_run "$@" &> "w_run_${SLURM_PROCID}.log"

# That might do a single segment or an entire iteration chunk,
# depending on how your WESTPA workflow is set up.

exit 0
