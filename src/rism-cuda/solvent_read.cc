#include <iostream>
#include <fstream>
#include <string>

#include "solvent.h"
#include "version.h"

void Solvent :: read(string fsolvent) {
  void alloc3D (vector <vector <double * > > &, int, int, int);
  
  ifstream in_file;
  in_file.open(fsolvent.c_str());

  string check;
  in_file >> check;
  if (check != version) {
    cout << "This xvva file is for old version." << endl;
    exit (1);
  }

  double dielr;
  in_file >> temper >> dielr >> xikt;

  in_file >> natv;

  multv = new int[natv];
  rhov = new double[natv];
  qv = new double[natv];
  sigv = new double[natv];
  epsv = new double[natv];

  for (int iv = 0; iv < natv; ++iv) {
    in_file >> multv[iv] >> rhov[iv] >> qv[iv] >> sigv[iv] >> epsv[iv];
  }

  double wlmvv;

  for (int iv1 = 0; iv1 < natv; ++iv1) {
    for (int imv = 0; imv < multv[iv1]; ++imv) {
      for (int iv2 = 0; iv2 < natv; ++iv2) {
	in_file >> wlmvv;
      }
    }
  }

  in_file >> ntab;

  ttab = new double[ntab];
    
  alloc3D (xvv, natv, natv, ntab);

  for (int n = 0; n < ntab; ++n) {
    in_file >> ttab[n];
    for (int iv2 = 0; iv2 < natv; ++iv2) {
      for (int iv1 = 0; iv1 < natv; ++iv1) {
	in_file >> xvv[iv2][iv1][n];
      }
    }
  }

  in_file.close();
}
