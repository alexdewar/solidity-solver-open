FROM ubuntu:20.04 AS solidity-build

# Install dependencies
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y --no-install-recommends --fix-missing build-essential libvtk7-dev

COPY src/* /home/software/solidity/src/
WORKDIR /home/software/solidity/src

# Build using all cores
RUN make -j `nproc`

# Target image
FROM ubuntu:20.04

# Install dependencies
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y --no-install-recommends --fix-missing libvtk7.1p

# Copy over executable and examples
COPY --from=solidity-build /home/software/solidity/src/bin/Release/DEMPlus /usr/local/bin
COPY examples/* /home/software/solidity/examples/

WORKDIR /home/software/solidity/run/
