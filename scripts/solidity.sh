#!/bin/sh

docker run --pull always --rm -it -u "$(id -u)":"$(id -g)" \
       -v "$(pwd)":/home/software/solidity/run/ \
       ghcr.io/imperialcollegelondon/solidity_project_dev:Solidity_RCS \
       DEMPlus "$@"