#!/usr/bin/env bash
#doesn't work rn probably due to WSL issue.
#conda activate my-rdkit-env
cd ../LogP || exit 1
python3 MakePolymer.py -j ../testing/test.json
