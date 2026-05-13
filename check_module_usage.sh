#!/bin/bash
# ParticleModule, TrajectoryModule and BundleModule
# override the assignment operator (`=`) for t_Particle, t_Trajectory and t_TrajectoryBundle
# derived types. However, this overriding doesn't work if the module is used
# with the `only: ` statement, and that statement does not include `assignment(=)` in its list.
# Because this is so error prone, we simply forbid using "only" with these modules,
# (while for all other modules we enforce the opposite via a fortitude linter).
#
# In the future, it might be better to not use the assignment like this at all,
# but that's not where we are at the moment.

set -euo pipefail

if [[ ${#@} -eq 0 ]]; then
    echo "ERROR: pass files to check as parameters to this script"
    exit 1
fi

if grep -Ein 'use (trajectory|particle|bundle)Module.*only' $@; then
    echo "ERROR: Found improper use of Particle|Trajectory|Bundle modules!"
    echo "Never use the 'only' keyword with these modules"
    exit 1
fi
