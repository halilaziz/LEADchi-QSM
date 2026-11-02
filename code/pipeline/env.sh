# env.sh - LEAD-chi environment, sourced by EVERY job on server (task P5b, 2026-09-25; revised 2026-09-26; registration Annex B.1).
#
# Rule: nothing may be written on the server outside /path/to/Penn_LEAD. Several tools write to
# $HOME by default (MATLAB ~/.matlab, apptainer ~/.apptainer and $HOME binds, TemplateFlow ~/.cache/templateflow,
# matplotlib ~/.config/matplotlib, QSMxT ~/.qsmxt and ~/.cache/qsm-rs), so every home/config/cache/tmp path is
# redirected into $LX/tools/. The real home is remembered once (REAL_HOME) only for the read-only home check.
#
# Sourcing this file creates NOTHING. Directories are created only by `lx_init_tree` (README step 3, after
# approval). Usage:  source /path/to/Penn_LEAD/LEADchi/scripts/pipeline/env.sh
# Dry-run convention for the shell scripts: LX_DRY=1 (or --dry-run) prints every command instead of running it.

export PL=/path/to/Penn_LEAD
if [ -n "${LX_TEST_PL:-}" ]; then                      # synthetic tests off-server only; refused on server
  case "$(hostname)" in server*) echo "LX_TEST_PL is not allowed on server" >&2; return 9 2>/dev/null || exit 9 ;; esac
  export PL=$LX_TEST_PL
fi
export LX=$PL/LEADchi
export LX_SCRIPTS=${LX_SCRIPTS:-$LX/scripts/pipeline}

# ---- remember the real home once, then redirect every home-like location into the tree
if [ -z "${REAL_HOME:-}" ]; then
  if [ "${HOME:-}" != "$LX/tools/home" ]; then export REAL_HOME=$HOME
  else REAL_HOME=$(getent passwd "$(id -un)" | cut -d: -f6); export REAL_HOME; fi
fi
export HOME=$LX/tools/home
export XDG_CONFIG_HOME=$LX/tools/xdg_config
export XDG_CACHE_HOME=$LX/tools/xdg_cache
export XDG_DATA_HOME=$LX/tools/xdg_data
export XDG_STATE_HOME=$LX/tools/xdg_state
export MPLCONFIGDIR=$LX/tools/mpl
export MATLAB_PREFDIR=$LX/tools/matlab_prefs
export MATLAB_LOG_DIR=$LX/tools/matlab_logs
export TMPDIR=$LX/tools/tmp TMP=$LX/tools/tmp TEMP=$LX/tools/tmp
export APPTAINER_CACHEDIR=$LX/tools/apptainer_cache
export APPTAINER_TMPDIR=$LX/tools/tmp
export APPTAINER_CONFIGDIR=$LX/tools/apptainer_config
export SINGULARITY_CACHEDIR=$APPTAINER_CACHEDIR SINGULARITY_TMPDIR=$APPTAINER_TMPDIR
export TEMPLATEFLOW_HOME=$LX/tools/templateflow
export APPTAINERENV_TEMPLATEFLOW_HOME=$LX/tools/templateflow
export APPTAINERENV_TMPDIR=$LX/tools/tmp
export APPTAINERENV_MPLCONFIGDIR=$LX/tools/mpl
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1

# ---- QSMxT v9: never self-update, never write ~/.qsmxt or ~/.cache/qsm-rs (E3 row 1.4)
export QSMXT_NO_UPDATE_CHECK=1
export QSMXT_UPDATE_CACHE=$LX/tools/qsmxt/update-check.json
export QSM_MODEL_CACHE=$LX/tools/qsm_models
export QSMXT_EXAMPLE_CACHE=$LX/tools/qsmxt/examples
export QSMXT_TAG=v9.22.0
export QSMXT_BIN=$LX/tools/qsmxt/qsmxt

# ---- CPU etiquette: nice and thread caps (<= 16 concurrent jobs; protocol sec. 12)
export LX_NICE=${LX_NICE:-10}
export LX_MAXJOBS=${LX_MAXJOBS:-16}
export LX_THREADS=${LX_THREADS:-2}                   # threads per job
export OMP_NUM_THREADS=$LX_THREADS OPENBLAS_NUM_THREADS=$LX_THREADS MKL_NUM_THREADS=$LX_THREADS
export NUMEXPR_NUM_THREADS=$LX_THREADS VECLIB_MAXIMUM_THREADS=$LX_THREADS
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$LX_THREADS
export RAYON_NUM_THREADS=$LX_THREADS                 # Rust (QSMxT/QSM.rs) worker pools [VERIFY at A1]
export ANTS_RANDOM_SEED=20260924
export APPTAINERENV_OMP_NUM_THREADS=$LX_THREADS APPTAINERENV_ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$LX_THREADS
umask 027

# ---- pinned tools (read/execute only outside the tree)
export ANTS_ROOT=/path/to/opt/ants-2.5.3
if [ -x "$ANTS_ROOT/bin/antsRegistration" ]; then export ANTSPATH=$ANTS_ROOT/bin/; else export ANTSPATH=$ANTS_ROOT/; fi
export FSLDIR=/path/to/opt/fsl-6.0.7.11-centos8_64      # default PATH resolves to FSL 5.0.6 (protocol sec. 12)
export FSLOUTPUTTYPE=NIFTI_GZ
export APPT=/path/to/apptainers           # A6: read/execute only
export QSIPREP_SIF=$APPT/qsiprep-26.0.0.sif          # SynthStrip on PATH; SynthSeg at /opt/freesurfer/bin/mri_synthseg
export SYNTHSEG_BIN=/opt/freesurfer/bin/mri_synthseg
export PY=${PY:-/path/to/apptainers/envs/AFQ_env/bin/python}   # numpy, scipy, nibabel (A6: execute only)
export SYSPY=${SYSPY:-/bin/python3}                                         # stdlib-only scripts
export OPENSSL=${OPENSSL:-$(command -v openssl || true)}
export PATH=$LX/tools/qsmxt:${ANTSPATH%/}:$FSLDIR/bin:$PATH

# ---- helpers ------------------------------------------------------------------------------------------------

# lx_init_tree: create the project tree (inside Penn_LEAD only). Run once, after approval (README step 3).
lx_init_tree() {
  local d
  for d in tools/home tools/xdg_config tools/xdg_cache tools/xdg_data tools/xdg_state tools/mpl tools/matlab_prefs \
           tools/matlab_logs tools/tmp tools/apptainer_cache tools/apptainer_config tools/templateflow tools/qsmxt \
           tools/qsm_models scripts/pipeline prereg/inputs locked blind reg qsm rois/templates rois/CIT168 rois/SUIT \
           rois/JHU rois/xfm rois/afq_templates tables_coded flags_coded qc freeze logs work gates dwi dwiqc; do
    lx_run mkdir -p "$LX/$d"
  done
  lx_run chmod 700 "$LX/locked" "$LX/blind" "$LX/tables_coded" "$LX/dwiqc"
  lx_run touch "$LX/logs/tree.stamp"
}

# lx_run CMD...: run a command under nice, or print it when LX_DRY=1.
lx_run() {
  if [ "${LX_DRY:-0}" = "1" ]; then printf '[DRY-RUN]'; printf ' %q' "$@"; printf '\n'; return 0; fi
  nice -n "$LX_NICE" "$@"
}

# lx_apx IMG CMD...: run CMD in a container with no home, a contained /tmp and only $LX bound (A6).
lx_apx() {
  local img=$1; shift
  lx_run apptainer exec --containall --no-home -B "$LX:$LX" --pwd "$LX" "$img" "$@"
}

# lx_assert_inside PATH: refuse any write target outside $LX.
lx_assert_inside() {
  local p r
  p=$("$SYSPY" -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$1")
  r=$("$SYSPY" -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$LX")
  case "$p/" in "$r"/*) return 0 ;; *) echo "REFUSED: $1 is outside $LX" >&2; exit 3 ;; esac
}

# lx_stamp NAME / lx_home_check NAME: read-only check that nothing new appeared in the real home (protocol 12).
lx_stamp() { [ "${LX_DRY:-0}" = "1" ] || touch "$LX/logs/$1.stamp"; }
lx_home_check() {
  [ "${LX_DRY:-0}" = "1" ] && { echo "[DRY-RUN] home check after $1"; return 0; }
  local hits; hits=$(find "$REAL_HOME" -maxdepth 3 -newer "$LX/logs/$1.stamp" 2>/dev/null | head -20)
  if [ -n "$hits" ]; then
    echo "STOP: new entries in $REAL_HOME after step $1 (nothing deleted; report to the user):" >&2
    echo "$hits" >&2
    return 4
  fi
  echo "home check after $1: clean"
}

# lx_require FILE...: abort unless every file exists (used for approval-gated inputs).
lx_require() {
  local f
  for f in "$@"; do [ -e "$f" ] || { echo "missing: $f" >&2; [ "${LX_DRY:-0}" = "1" ] || exit 5; }; done
}
