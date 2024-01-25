# XrayPhysics
This is a C/C++ library and Python bindings with the following capabilities:
1) x-ray cross sections (energies from 1 keV to 20 MeV and elements 1-100)
2) x-ray source spectra modeling (for any source voltage and take-off angle, and the following anode types: Cu, Mo, W, Au)
3) single-material beam hardening correction algorithms (theoretically exact and polynomial-based)
4) Calculation of the LLNL-defined effective atomic number; For more information see [here](https://ieeexplore.ieee.org/document/8638824)
5) dual energy decomposition algorithms

And here are the features we are working on for the next relase:
1) multi-material beam hardening correction algorithms
2) incoherent and coherent scattering angle distributions

![The x-ray source models are at least as accurate as TASMICS and much more flexible.](https://github.com/kylechampley/XrayPhysics/blob/main/comparisonWithTASMICS.png)

The x-ray source models are at least as accurate as TASMICS and much more flexible.


## Installation and Usage

Installation and usage information is posted on the wiki page here: https://github.com/kylechampley/XrayPhysics/wiki

Designed to be used in conjuction with (but is not required) [LEAP](https://github.com/LLNL/LEAP) 

## Author
Kyle Champley (champley@gmail.com)


## License
XrayPhysics is distributed under the terms of the MIT license. All new contributions must be made under this license. See LICENSE in this directory for the terms of the license.
See [LICENSE](LICENSE) for more details.  
SPDX-License-Identifier: MIT  

Please cite our work by referencing this github page

