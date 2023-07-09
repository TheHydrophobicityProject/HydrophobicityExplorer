#!/bin/bash
#run from the directory containing setup.py

python -m build

python3 -m twine upload --repository pypi dist/*

#enter login info

