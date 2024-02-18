# XrayPhysics
This is a C/C++ library with Python bindings with the following capabilities:
1) x-ray cross sections (energies from 1 keV to 20 MeV and elements 1-100) specified by chemical formula or element/compound mass fractions
2) incoherent and coherent x-ray scattering angle distributions specified by chemical formula or element/compound mass fractions
3) x-ray source spectra modeling (for any source voltage and take-off angle, and the following anode types: Cu, Mo, W, Au)
4) one- and two-material beam hardening correction algorithms (theoretically exact and polynomial-based)
5) Calculation of the LLNL-defined effective atomic number; For more information see [here](https://ieeexplore.ieee.org/document/8638824)
6) dual energy decomposition algorithms (and SIRZ if used in conjuction with [LEAP](https://github.com/LLNL/LEAP))

There are very few dependencies and we don't use any specialized data structures, making it very easy to incorporate with other software packages.  The cross section tables are hard-coded into C++ arrays, so queries of the database are instant.

![The x-ray source models are at least as accurate as TASMICS and much more flexible.](https://github.com/kylechampley/XrayPhysics/blob/main/comparisonWithTASMICS.png)

The x-ray source models are at least as accurate as TASMICS and much more flexible.  The cross section tables are based on [EPDL97.](https://www-nds.iaea.org/epdl97/)


## Installation and Usage

Installation and usage information is posted on the [wiki page](https://github.com/kylechampley/XrayPhysics/wiki)

Designed to be used in conjuction with (but is not required) [LEAP](https://github.com/LLNL/LEAP) 

## Author
Kyle Champley (champley@gmail.com)


## License
XrayPhysics is distributed under the terms of the MIT license. All new contributions must be made under this license. See LICENSE in this directory for the terms of the license.
See [LICENSE](LICENSE) for more details.  
SPDX-License-Identifier: MIT  

Please cite our work by referencing this github page and the following [paper](https://ieeexplore.ieee.org/document/8638824):

Champley, Kyle M., Stephen G. Azevedo, Isaac M. Seetho, Steven M. Glenn, Larry D. McMichael, Jerel A. Smith, Jeffrey S. Kallman, William D. Brown, and Harry E. Martz. "Method to extract system-independent material properties from dual-energy X-ray CT." IEEE Transactions on Nuclear Science 66, no. 3 (2019): 674-686.


