#!/bin/bash
set -euo pipefail
mkdir -p /root/results
cp /solution/synthesize.py /root/results/synthesize.py
# check the program on the examples
for f in /root/data/examples/[A-Z]*.json; do
    python3 /root/results/synthesize.py "$f" "/tmp/$(basename "$f" .json).network.json"
done
