#include <iostream>
#include <fstream>
#include "rism3d.h"

void RISM3D :: set_solvent () {
  sv -> read(fsolvent);
  sv -> spline(ga, indga, nga, ce -> ngrid);
}
