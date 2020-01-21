#!/usr/bin/env zsh

# Trying to run this as a zsh script didn't work because I don't think I have my .zshrc working properly
# It couldn't find mkvirtualenv

# Create a virtualenv with the desired packages
# and set the project directory for `workon`.
mkvirtualenv -r requirements.txt -a . pandas-learning
