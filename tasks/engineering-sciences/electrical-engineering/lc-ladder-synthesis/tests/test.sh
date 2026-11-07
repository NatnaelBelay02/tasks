#!/bin/bash
set -euo pipefail

mkdir -p /logs/verifier
echo 0 > /logs/verifier/reward.txt

# the held-out specifications and the reference must not be readable by the program under test
chmod -R go-rwx /tests

if pytest --ctrf /logs/verifier/ctrf.json /tests/test_synthesis.py -rA -p no:cacheprovider; then
    echo 1 > /logs/verifier/reward.txt
else
    echo 0 > /logs/verifier/reward.txt
fi
