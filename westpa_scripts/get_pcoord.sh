#!/bin/bash
##############################################################################
# Example driver script to compute RMSD and Rg for each frame using MDAnalysis
# for a single WESTPA segment. It mimics the cpptraj-based script but replaces
# it with Python+MDAnalysis for analysis.
#
# Expected WESTPA environment variables:
#   - $WEST_SIM_ROOT        : Path to the simulation root.
#   - $WEST_STRUCT_DATA_REF : Trajectory file (DCD, etc.) for this segment.
#   - $WEST_PCOORD_RETURN   : File to which pcoord data must be written.
#   - $SEG_DEBUG (optional) : If set, script runs in debug mode (prints env, etc.).
#
# This script assumes your trajectory contains ~101 frames (or more),
# matching the shape (101, 2) that WESTPA expects for pcoords:
# i.e., 101 lines with [RMSD, Rg].
##############################################################################

# 1. Load Python (or activate a conda environment with MDAnalysis)
#module load python/3.9
# Or: conda activate my_mdanalysis_env

# 2. If debugging is enabled, turn on shell debugging and print environment
if [ -n "$SEG_DEBUG" ]; then
  set -x
  env | sort
fi

# 3. Move to the simulation root directory
cd "$WEST_SIM_ROOT" || {
  echo "Error: Could not cd to \$WEST_SIM_ROOT=$WEST_SIM_ROOT" >&2
  exit 1
}

# 4. Create temporary files for storing RMSD and Rg data
RMSD_FILE=$(mktemp --tmpdir rmsd_XXXX.xvg)
RG_FILE=$(mktemp --tmpdir rg_XXXX.xvg)

# 5. Create a small Python script to run MDAnalysis
cat << EOF > mdanalysis_rmsd_rg.py
#!/usr/bin/env python

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import rms

# Input files
topology = "/ocean/projects/chm240001p/thyagatu/pargamd/ParGaMD/common_files/chignolin.parm7"
trajectory = "${WEST_STRUCT_DATA_REF}"
ref_pdb = "/ocean/projects/chm240001p/thyagatu/pargamd/ParGaMD/common_files/chignolin.pdb"

# Load the reference and the simulation trajectory
ref = mda.Universe(ref_pdb)
u = mda.Universe(topology, trajectory)

# Select CA atoms
ref_ca = ref.select_atoms("name CA")
mobile_ca = u.select_atoms("name CA")

# Prepare lists for RMSD and Rg
rmsd_values = []
rg_values = []

# Loop through frames in the trajectory
for ts in u.trajectory:
    # -- Mass-weighted best-fit RMSD to the reference CA --
    #    center=True, superposition=True => best-fit alignment
    #    weights=mobile_ca.masses => mass-weighted
    this_rmsd = rms.rmsd(
        mobile_ca.positions,
        ref_ca.positions,
        center=True,
        superposition=True,
        weights=mobile_ca.masses
    )

    # -- Mass-weighted radius of gyration for CA atoms --
    #    By default, radius_of_gyration() is geometric (not mass-weighted).
    #    We'll compute a mass-weighted Rg manually:
    coords = mobile_ca.positions
    masses = mobile_ca.masses
    com = np.average(coords, axis=0, weights=masses)
    sq_dist = np.sum((coords - com)**2, axis=1)
    rg = np.sqrt(np.average(sq_dist, weights=masses))

    rmsd_values.append(this_rmsd)
    rg_values.append(rg)

# Write out data to match the cpptraj style .xvg (frame vs. value).
# We will keep a simple two-column format: "frame value".
with open("${RMSD_FILE}", "w") as f_rms:
    f_rms.write("# frame RMSD(Angstrom)\n")
    for i, val in enumerate(rmsd_values):
        f_rms.write(f"{i+1} {val}\n")

with open("${RG_FILE}", "w") as f_rg:
    f_rg.write("# frame Rg(Angstrom)\n")
    for i, val in enumerate(rg_values):
        f_rg.write(f"{i+1} {val}\n")
EOF

chmod +x mdanalysis_rmsd_rg.py

# 6. Run the Python script
if ! ./mdanalysis_rmsd_rg.py; then
  echo "Error: MDAnalysis script execution failed!" >&2
  rm -f "$RMSD_FILE" "$RG_FILE" mdanalysis_rmsd_rg.py
  exit 1
fi

# 7. Extract the RMSD and Rg columns from all frames and paste them
#    We skip the header line (#) in each .xvg file with 'NR>1'.
#    $2 is the second column (the actual RMSD or Rg value).
paste <(awk 'NR>1 {print $2}' "$RMSD_FILE") \
      <(awk 'NR>1 {print $2}' "$RG_FILE") \
      > "$WEST_PCOORD_RETURN"

# 8. If debug is enabled, show a quick preview of the resulting pcoord file
if [ -n "$SEG_DEBUG" ]; then
  echo "Preview of \$WEST_PCOORD_RETURN:"
  head -v "$WEST_PCOORD_RETURN"
  echo "Number of lines in \$WEST_PCOORD_RETURN: \$(wc -l < "$WEST_PCOORD_RETURN")"
fi

# 9. Clean up temporary files and exit
rm -f "$RMSD_FILE" "$RG_FILE" mdanalysis_rmsd_rg.py
exit 0
