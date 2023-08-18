#!/bin/bash
#run from the directory containing setup.py

version_number=$(awk -F "\"" '/version/ {print $2}' setup.py)

python3 -m build

find dist -name "*${version_number}*" -exec python3 -m twine upload --repository pypi "{}" \;

#enter login info
