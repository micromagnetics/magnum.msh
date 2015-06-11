magnum.msh
==========

magnum.msh is a meshing helper for magnum.fe (http://micromagnetics.org/magnum.fe). It uses Gmsh to import a variety of mesh file formats and offers functionality to create cuboid shells as required by the shell transformation method for open boundary problems.

Installation
------------
### Prerequisites
magnum.msh requires the following software/libraries to be installed:

* FEniCS >= 1.5
* Gmsh >= 2.8.0 (headers and library)
* CMake >= 2.8
* SWIG >= 2.0
* G++ >= 4.0

#### Install dependencies in Ubuntu 14.04
Install FEniCS

    $ sudo add-apt-repository ppa:fenics-packages/fenics
    $ sudo apt-get update
    $ sudo apt-get install fenics

Install Gmsh with headers and build dependencies

    $ sudo apt-get install libgmsh-dev g++ swig cmake

### Build and install
To build magnum.msh with cmake do

    $ cd /path/to/magnum.msh
    $ mkdir build
    $ cd build
    $ cmake ..
    $ make
    $ sudo make install

License and Disclaimer
----------------------
magnum.msh is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

magnum.msh is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.
