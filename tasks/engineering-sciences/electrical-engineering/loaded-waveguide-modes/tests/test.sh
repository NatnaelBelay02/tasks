#!/bin/bash
set -euo pipefail

mkdir -p /logs/verifier
echo 0 > /logs/verifier/reward.txt

if pytest --ctrf /logs/verifier/ctrf.json /tests/test_modes.py -rA; then
    echo 1 > /logs/verifier/reward.txt
else
    echo 0 > /logs/verifier/reward.txt
fi
