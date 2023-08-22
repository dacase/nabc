#include <iostream>
#include "rism3d.h"

// calculate and return P/kBT as pressure

double RISM3D :: cal_pressure () {

  double pressure = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
    pressure +=  sv -> rhov[iv];
  }

  pressure *= 1.0e24 / avogadoro;

  pressure = 0.5 * (pressure + 1.0/(sv -> xikt));
// P/kBT (pressure ) is in [mol/cc]

  return pressure;
}
