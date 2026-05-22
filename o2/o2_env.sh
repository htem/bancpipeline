#!/bin/bash
# Common O2 environment for bancpipeline SBATCH jobs

# sbatch runs a non-interactive, non-login shell where LMOD's `module`
# function is not defined. Source it explicitly so `module load` works
# regardless of how this script is invoked.
# On O2 the LMOD init lives at /etc/profile.d/modules.sh (sometimes
# /etc/profile.d/lmod.sh on other clusters). The cron-launched updateids
# job kept failing silently for ~3 weeks because only the lmod.sh path was
# being checked — see banc_updateids_3*.err for evidence.
if ! command -v module >/dev/null 2>&1; then
  for f in /etc/profile.d/modules.sh /etc/profile.d/lmod.sh /etc/profile.d/00-modulepath.sh /etc/profile.d/z00_lmod.sh; do
    [ -r "$f" ] && . "$f"
  done
fi

module purge
module load gcc/14.2.0
module load R/4.4.2
module load cmake/3.31.2
module load java/jdk-23.0.1

export UDUNITS2_INCLUDE=/n/app/udunits/2.2.28-gcc-9.2.0/include
export UDUNITS2_LIBS=/n/app/udunits/2.2.28-gcc-9.2.0/lib

# Disable core dumps. Some R prep scripts hit `free(): invalid pointer` in
# libc teardown *after* all output files are written (benign, tolerated by
# our sbatch scripts via `|| true` + require_files). Without this ulimit the
# kernel dumps multi-GB core.<pid> files into the repo on every such exit.
ulimit -c 0
