#!/usr/bin/env zsh

# Have to source zshrc to get access to virtualenvwrapper commands
source ~/.zshrc

# Create a virtualenv with the desired packages
# and set the project directory for `workon`.
mkvirtualenv -r ~/devel/learning/python/tutorials/requirements.txt -a ~/devel/learning/python/tutorials/ py-tutorials
