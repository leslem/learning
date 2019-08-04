#!/usr/bin/env bash

# Create a virtualenv with the desired packages
# and set the project directory for `workon`.
mkvirtualenv -r requirements.txt -a . --python=python3.7 ds-from-scratch

# Compile a pinned requirements file.
pip-compile requirements.in

# Install packages.
pip-sync requirements.txt
