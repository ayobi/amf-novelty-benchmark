#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
python3 -u run_pilot.py --threads "${AMF_THREADS:-8}" 2>&1 | tee -a pilot_run.log
