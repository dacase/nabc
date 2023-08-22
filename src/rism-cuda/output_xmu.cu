#include <iostream>
#include <fstream>
#include <iomanip>
#include "rism3d.h"
#include "extension.h"

void RISM3D :: output_xmu(double * & xmu, double pmv, double pressure) {
    
  ofstream out_file;
  out_file.open((fname + extxmu).c_str());

  double ibeta = avogadoro * boltzmann * sv -> temper;

  double xmua = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
    xmua += xmu[iv];
  }

  double gf = 0.0;
  for (int iv = 0; iv < sv -> natv; ++iv) {
    gf += xmu[sv -> natv + iv];
  }

  double pcterm = - pressure * pmv * ibeta;

  out_file << "SFE(SC)= " << fixed << setprecision(5) 
  	   << ibeta * xmua << " !(J/mol)" << endl;

  for (int iv = 0; iv < sv -> natv; ++iv) {
    out_file << "  " << iv << "= " << fixed << setprecision(5) 
             << ibeta * xmu[iv] << endl;
  }
  out_file << endl;

  out_file << "SFE(GF)= " << fixed << setprecision(5) 
	   << ibeta * gf << " !(J/mol)" << endl;

  for (int iv = 0; iv < sv -> natv; ++iv) {
    out_file << "  " << iv << "= " << ibeta * xmu[sv -> natv + iv] 
	     << endl;
  }
  out_file << endl;

  out_file << "PMV= " << fixed << setprecision(5)
           << pmv << " !(cc/mol)" << endl;

  out_file << "Pressure= " << fixed << setprecision(5)
           << pressure * kcal2J << " !(J/cc)" << endl;

  out_file << "Correction_Term= " << fixed << setprecision(5) 
             << pcterm << " !(J/mol)" << endl;

  out_file.close();
} 
