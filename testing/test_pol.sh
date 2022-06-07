#!/usr/bin/env bash
#doesn't work rn probably due to WSL issue.
#conda activate my-rdkit-env
cd ../LogP
python3 MakePolymer.py -j ../testing/test.json
