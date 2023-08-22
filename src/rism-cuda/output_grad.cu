#include <iostream>
#include <fstream>
#include "rism3d.h"
#include "extension.h"

void RISM3D :: output_grad(double * & du) {

  cout << "outputting grad to file:  " << fname + extgra << "  ..." << endl;

  ofstream out_file;
  out_file.open ((fname + extgra).c_str());

  double dv = ce -> dv / kcal2J;
  for (int iu = 0; iu < su -> num; ++iu) {
    int num = iu * 3;
    out_file << du[num] * dv << " "
	     << du[num + 1] * dv << " "
	     << du[num + 2] * dv << endl;
  }

  cout << "done." << endl;

  out_file.close();
} 
