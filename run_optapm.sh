#!/bin/bash
#
# run_optapm.sh - one command to run the full optAPM uncertainty envelope for a plate model.
#
# Runs both optimisation variants and their diagnostics:
#   1. 'best'   - net rotation bounded to 0.08-0.20 deg/Myr           (the preferred frame)
#   2. 'nr_max' - net rotation 0.08-0.30 deg/Myr, half NR weight      (high-NR end-member)
# The no-net-rotation (zero-NR) frame is produced automatically by each run, and each run also
# writes per-timestep diagnostics (net rotation + trench migration) to model_output/.
#
# USAGE:
#   ./run_optapm.sh <model> [cores] [start_age] [end_age]
#
#   <model>     EITHER the name of a model already under data/ (a "data_model", e.g.
#               Zahirovic_etal_2022_GDJ), OR a path to a plate-model directory anywhere on disk.
#               A directory that is NOT already under data/ is symlinked into data/ automatically,
#               and its rotation (.rot) and topology (.gpml/.gpmlz) files are auto-discovered
#               (searched recursively) - no editing of Optimised_config.py is needed.
#   cores       number of CPU cores for mpiexec                       (default 8)
#   start_age   oldest age to optimise, in Ma  (required for a model with no built-in default age)
#   end_age     youngest age, in Ma                                   (default 0)
#
# EXAMPLES:
#   ./run_optapm.sh Zahirovic_etal_2022_GDJ 8
#   ./run_optapm.sh /data/plate_models/MyModel 16 600 0
#
# OUTPUTS (in the model's optimisation/ sub-directory):
#   optimised_rotation_model_<name>.rot, optimised_rotation_model_<name>_nr_max.rot,
#   no_net_rotation_model_<name>.rot
# LOGS + DIAGNOSTICS (in model_output/):
#   <name>_run.log, <name>_nr_max_run.log, <name>[_nr_max]_diagnostics.csv + *.png
#
# REFERENCE FRAME: auto-discovered models optimise Africa (plate 701), seeded from the previous
#       optimised interval. To change this without editing Optimised_config.py, set these environment
#       variables before calling the script (they are inherited by the run):
#         OPTAPM_REF_PLATE_ID=<id>                          reference plate, e.g. 101
#         OPTAPM_REF_ROTATION_FILE=<optimised|nnr|path>     search seed (path is relative to data/)
#       e.g.:  OPTAPM_REF_PLATE_ID=101 ./run_optapm.sh /data/plate_models/MyModel 8 1000 0
#       A model that is not already in its paleomagnetic frame needs a pmag seed file pointed to by
#       OPTAPM_REF_ROTATION_FILE. If the model carries its pmag path on a plate you would anchor
#       (eg, 701701), generate that seed file first with:
#         python make_pmag_seed.py --pmag-plate 701701 --ref-plate 701 --output data/<model>/pmag/seed.rot
#
set -e
cd "$(dirname "$0")"

if [ -z "$1" ]; then
    sed -n '7,33p' "$0"   # print the usage block above
    exit 1
fi

MODEL="$1"
CORES="${2:-8}"
START_AGE="${3:-}"
END_AGE="${4:-}"

# Resolve <model> to a data_model name. If it is a directory path, symlink it into data/ (keeping
# the shared data/ directory as the root, so the bundled hotspot data still resolves) and use its
# basename as the data_model.
if [ -d "$MODEL" ]; then
    ABS=$(cd "$MODEL" && pwd -P)
    NAME=$(basename "$ABS")
    if [ ! -e "data/$NAME" ]; then
        ln -s "$ABS" "data/$NAME"
        echo "Linked data/$NAME -> $ABS"
    fi
    export OPTAPM_DATA_MODEL="$NAME"
else
    export OPTAPM_DATA_MODEL="$MODEL"
fi
[ -n "$START_AGE" ] && export OPTAPM_START_AGE="$START_AGE"
[ -n "$END_AGE" ]   && export OPTAPM_END_AGE="$END_AGE"

mkdir -p model_output

echo "=========================================================="
echo " optAPM run  |  model: $OPTAPM_DATA_MODEL  |  cores: $CORES"
[ -n "$START_AGE" ] && echo " age range: ${START_AGE}-${END_AGE:-0} Ma"
echo " $(date)"
echo "=========================================================="

# --- Run 1/2: 'best' ------------------------------------------------------------------------
unset OPTAPM_VARIANT
echo ""
echo "---- Run 1/2: best (NR 0.08-0.20 deg/Myr) ----"
mpiexec -n "$CORES" python Optimised_APM.py
echo "---- diagnostics: best ----"
python model_diagnostics.py

# --- Run 2/2: 'nr_max' ----------------------------------------------------------------------
export OPTAPM_VARIANT=nr_max
echo ""
echo "---- Run 2/2: nr_max (NR 0.08-0.30 deg/Myr, half NR weight) ----"
mpiexec -n "$CORES" python Optimised_APM.py
echo "---- diagnostics: nr_max ----"
python model_diagnostics.py
unset OPTAPM_VARIANT

echo ""
echo "=========================================================="
echo " Done. Rotation models in data/$OPTAPM_DATA_MODEL/optimisation/ ;"
echo " logs and diagnostics in model_output/."
echo " To assemble all reference frames into one file, next run:"
echo "     python add_ref_frames_to_unoptimised.py"
echo "=========================================================="
