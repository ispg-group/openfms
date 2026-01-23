#!/bin/bash

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
