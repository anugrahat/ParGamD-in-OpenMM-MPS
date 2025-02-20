#!/bin/bash
set -x  # Enable debugging output

module load cuda/11.8.0

# Optional debugging if SEG_DEBUG is set
if [ -n "$SEG_DEBUG" ]; then
    env | sort
    nvidia-smi || echo "Warning: nvidia-smi failed"
    python -c "import openmm; print('OpenMM version:', openmm.__version__)" || echo "Warning: OpenMM import failed"
fi

##############################################################################
# STEP 1: Dynamically Assign Workers to GPUs Across Multiple Nodes
##############################################################################
WORKER_INDEX=${WM_WORKER_INDEX:-0}   # Worker index (varies dynamically)
TOTAL_WORKERS=${WM_N_WORKERS:-$(scontrol show job $SLURM_JOBID | grep NumTasks | awk '{print $3}' | cut -d= -f2)}
NUM_NODES=${SLURM_NNODES:-1}        # Number of nodes allocated
GPUS_PER_NODE=$(nvidia-smi -L | wc -l)  # Detect GPUs per node
TOTAL_GPUS=$(( NUM_NODES * GPUS_PER_NODE ))  # Total GPUs across all nodes

# Ensure we don't divide by zero
if [[ "$TOTAL_GPUS" -eq 0 ]]; then
    echo "Error: No GPUs detected!"
    exit 1
fi

# Assign worker to a GPU based on total workers & GPUs
NODE_RANK=$(( WORKER_INDEX / (TOTAL_WORKERS / NUM_NODES) ))  # Node index
GPU_ID=$(( WORKER_INDEX % GPUS_PER_NODE ))  # Assign GPU within the node

echo "Worker $WORKER_INDEX => Assigned GPU $GPU_ID on Node $NODE_RANK"
export CUDA_VISIBLE_DEVICES=$GPU_ID  # Assign correct GPU

# Verify CUDA access
nvidia-smi
python -c "import openmm; print('OpenMM CUDA platform:', openmm.Platform.getPlatformByName('CUDA'))" || echo "Warning: OpenMM CUDA not accessible"

##############################################################################
# STEP 2: Create and Move into the Simulation Directory
##############################################################################
mkdir -pv "$WEST_CURRENT_SEG_DATA_REF"
cd "$WEST_CURRENT_SEG_DATA_REF" || exit 1

##############################################################################
# STEP 3: Link Necessary Files (Topology, Coordinates, XML)
##############################################################################
ln -sfv "$WEST_SIM_ROOT/common_files/chignolin.parm7" .
ln -sfv "$WEST_SIM_ROOT/common_files/gamd-restart.dat" .
ln -sfv "$WEST_SIM_ROOT/common_files/upper_dual_temp.xml" .
ln -sfv "$WEST_SIM_ROOT/common_files/chignolin.rst7" .

# Fix XML output directory to the current location
sed -i 's|<directory>.*</directory>|<directory>.</directory>|' upper_dual_temp.xml

##############################################################################
# STEP 4: Handle GaMD Checkpoints (First vs. Subsequent Iterations)
##############################################################################
if [ "$WEST_CURRENT_ITER" -eq 1 ]; then
    echo "First iteration - copying initial checkpoint file"
    cp -v "$WEST_SIM_ROOT/common_files/gamd_restart.checkpoint" ./gamd_restart.checkpoint
else
    echo "Subsequent iteration - using parent checkpoint"
    cp -v "$WEST_PARENT_DATA_REF/gamd_restart.checkpoint" ./gamd_restart.checkpoint
fi

##############################################################################
# STEP 5: Run the GaMD Simulation (OpenMM)
##############################################################################
echo "Starting GaMD simulation..."
python "$WEST_SIM_ROOT/common_files/gamdRunner" -r -p CUDA xml upper_dual_temp.xml 

##############################################################################
# STEP 6: Verify That GaMD Output Was Generated
##############################################################################
if [ ! -f "gamd_restart.checkpoint" ] || [ ! -f "output_restart.dcd" ]; then
    echo "Error: GaMD simulation failed to generate output files"
    ls -l
    exit 1
fi

##############################################################################
# STEP 7: Post-Processing with cpptraj (RMSD + Radius of Gyration)
##############################################################################
RMSD_FILE="rmsd_ca.xvg"
RG_FILE="rg_ca.xvg"
CPPTRAJ_LOG="cpptraj.log"

cat <<EOF | cpptraj > "$CPPTRAJ_LOG" 2>&1
parm chignolin.parm7
trajin output_restart.dcd
reference $WEST_SIM_ROOT/common_files/chignolin.pdb
rms rmsd_ca @CA reference out $RMSD_FILE mass
radgyr rg_ca @CA out $RG_FILE
go
EOF

if [ $? -ne 0 ]; then
    echo "Error: cpptraj failed. Check $CPPTRAJ_LOG for details."
    exit 1
fi

##############################################################################
# STEP 8: Write Final RMSD and Rg Data to $WEST_PCOORD_RETURN
##############################################################################
if [ -f "$RMSD_FILE" ] && [ -f "$RG_FILE" ]; then
    > "$WEST_PCOORD_RETURN"
    paste <(awk 'NR>1 {print $2}' "$RMSD_FILE") <(awk 'NR>1 {print $2}' "$RG_FILE") >> "$WEST_PCOORD_RETURN"
else
    echo "Error: Missing $RMSD_FILE or $RG_FILE"
    exit 1
fi

##############################################################################
# STEP 9: Optional Debugging
##############################################################################
if [ -n "$SEG_DEBUG" ]; then
    echo "Preview of $WEST_PCOORD_RETURN:"
    head -v "$WEST_PCOORD_RETURN"
fi

echo "Segment run completed successfully."
