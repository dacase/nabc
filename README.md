*Work in progress:*

* Jargon-filled summary
This is conversion of AmberTool's NAB/SFF functionality, to keep the "SFF" part,
but convert the "NAB" part to C.

* More complete summary
The AmberTools package (https://ambermd.org) contains a molecular
maniupulation language, NAB (Nucleic Acid Builder, written primarily by Tom
Macke) that interfaces to SFF (Simple Force Field, created originally by
Dave Case).

The NAB language is expressive for some things, especially for short
programs -- its original design was to serve as a "molecular AWK".  But it
doesn't scale well to larger projects, requires users to learn a new,
C-like, language that isn't quite C, and doesn't do a good job of catching
syntax errors.

The companion SFF library is written in C, and has a lot of nice features.
It provides a very lightweight code for carrying out simulations of
non-periodic systems, and is file-compatible with the Amber package of
programs.

This is an attempt to create driver programs in C that interface with SFF
API.  Right now it is a one-person (DAC) effort, but let me know if you
would like to help out.  

* Primary Contributors

**  Tom Macke designed the language and wrote the basic code
**  Dave Case provided the force field routines
**  Andreas Svrcek-Seiler contributed in general, and in particular to the GB code
**  Russ Brown wrote the second-derivative code
**  István Kolossváry contributed the "low-mode" code
**  Yannick Bomble worked on normal modes and Langevin modes
**  Jason Swails added support for variable 1-4 scaling and GB neck models
**  Ramu Anandakrishnan and Alexey V. Onufriev added the hierarchical charge partitioning approximation. 

* License
Mostly GPLv3; see the LICENSE file for full details.
