#include <iostream>
#include "rism3d.h"

void RISM3D :: initialize(string control, string structure, bool centering) {
  read_input(control, structure, centering);
  set_cuda();
  set_fname(control, structure);
  initialize_g();
  set_solvent();
  // DAC: since cal_potential depends on the coordinates, do this inside
  //     mme_rism()
  // cal_potential();
} 
