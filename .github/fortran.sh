#!/usr/bin/sh

gfortran-9 --version
gfortran-10 --version
sudo apt-get install -y gcc-11
sudo apt-get install -y gfortran-11
gfortran-11 --version
clang --version
sudo apt-get install -y flang libgmp-dev libomp-11-dev
flang --version
