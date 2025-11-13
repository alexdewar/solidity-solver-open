$ErrorActionPreference = "Stop"
docker run --pull always --rm -it -v ${pwd}:/home/software/solidity/run/ ghcr.io/imperialcollegelondon/solidity_project_dev:Solidity_RCS DEMPlus ${args}