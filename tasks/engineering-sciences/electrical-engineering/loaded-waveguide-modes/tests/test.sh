#!/bin/bash
set -euo pipefail
mkdir -p /logs/verifier
# Fail safe: reward is 0 unless the full test suite passes.
echo 0 > /logs/verifier/reward.txt
chmod 700 /logs/verifier
chmod -R go-rwx /tests
if timeout 4200 python3 -m pytest -p no:cacheprovider --ctrf /logs/verifier/ctrf.json /tests/test_modes.py -rA -s; then
    echo 1 > /logs/verifier/reward.txt
else
    echo 0 > /logs/verifier/reward.txt
fi
