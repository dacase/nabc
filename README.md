
# Jargon-filled summary
This is conversion of AmberTool's NAB/SFF functionality, to keep the "SFF" part,
but convert the "NAB" part to C.  Optionally, one can also install the
original NAB compiler.

# More complete summary
The AmberTools package (https://ambermd.org) contains a molecular
maniupulation language, NAB (Nucleic Acid Builder, written primarily by Tom
Macke) that interfaces to SFF (Simple Force Field, created originally by
Dave Case).

The NAB language is expressive for some things, especially for short
programs -- its original design was to serve as a "molecular AWK".  But it
doesn't scale well to larger projects, requires users to learn a new,
C-like, language that isn't quite C, and doesn't do a good job of catching
syntax errors at compile time.

The companion SFF library is written in C, and has a lot of nice features.
It provides a lightweight code for carrying out simulations of non-periodic
systems, and is file-compatible with the Amber package of programs.  
All of the Amber generalized Born options are availabe, along with
some unique (or at least, unusual) features, including:

* generalized-Born second derivatives
* normal modes and Langevin modes
* "low-mode" conformational search and docking routines
* an optimized nonbonded list generator
* CPU parallelization via both MPI and OpenMP.

This is an attempt to create driver programs in C that interface with SFF
API.  I have also included a slightly updated version of the NAB compiler
itself.  Right now it is a one-person (DAC) effort, but let me know if you
would like to help out.  

# Primary Contributors

*  Tom Macke designed the language and wrote the basic code
*  Dave Case provided the force field routines
*  Andreas Svrcek-Seiler contributed in general, and in particular to the GB code
*  Russ Brown wrote the second-derivative code, and the list generator
*  István Kolossváry contributed the "low-mode" code
*  Yannick Bomble worked on normal modes and Langevin modes
*  Jason Swails added support for variable 1-4 scaling and GB neck models
*  Ramu Anandakrishnan and Alexey V. Onufriev added the hierarchical charge partitioning approximation. 

# License
Mostly GPLv3; see the LICENSE file for full details.

# How to install

*  ./configure --help   # then run configure with the options you want; no options is the most common choice
*  make install   # creates lib/lib{arpack,nabc,sff,lapack,blas}.a
*  cd test && make test  

*  make nab  # optional: only if you want or need the NAB compiler
*  cd test/nab && make test  # only makes sense if you installed NAB

# Documentation

* See the doc/nabc.pdf file -- mostly OK, but there are some out-of-date things there!
* Look at the example driver files in the test folder; see if you can modify those to meet your needs.

# Reporting problems

* Best is to create an issue at github.com/dacase/nabc
* I've run tests on Ubuntu 22.04 and OSX Ventura (with HomeBrew gcc, intel hardware)
