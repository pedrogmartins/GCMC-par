#!/bin/bash

#Script to push all added scripts to this directory

git add .
git commit -m "push-update"
git branch -M main
git push -u origin main
