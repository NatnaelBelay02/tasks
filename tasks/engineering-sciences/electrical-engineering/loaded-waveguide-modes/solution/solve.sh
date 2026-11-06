#!/bin/bash
set -euo pipefail
cd /root
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 /solution/solve.py /root/data/guides.json /root/results/modes.json
