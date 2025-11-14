# Solidity Solver (Open Source)

This repository contains the **public, open-source** version of the Solidity solver.  

## Getting Started

### Prerequisities

In order to run Solidity you will need docker installed.

* [Windows](https://docs.docker.com/windows/started)
* [OS X](https://docs.docker.com/mac/started/)
* [Linux](https://docs.docker.com/linux/started/)

### Usage

On Mac/Linux run:

```shell
./scripts/solidity.sh path/to/input/files
```

On Windows (Powershell) run:

```shell
.\scripts\solidity.ps1 path\to\input\files
```

### Running examples

To run Solidity on the example data, on Mac/Linux run:

```shell
cd examples
../scripts/solidity.sh BST.Y3D layout.txt pellet.txt mat.txt
```

On Windows (Powershell) run:

```shell
cd examples
..\scripts\solidity.ps1 BST.Y3D layout.txt pellet.txt mat.txt
```

## Development Workflow

### 1. Branching Strategy
```
main # Stable version
dev # Active development
```

### 2. Contributing
1. Fork the repository.
2. Create a branch from `dev`.
3. Submit a PR with:
   - tests
   - documentation updates
   - clear description of the change