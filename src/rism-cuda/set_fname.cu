#include <iostream>
#include <fstream>
#include "rism3d.h"

void RISM3D :: set_fname (string control, string structure) {
  if (structure.empty()) {
    fname.append(control);
  } else {
    fname.append(structure);
  }
  fname = fname.substr(0, fname.rfind("."));
}
