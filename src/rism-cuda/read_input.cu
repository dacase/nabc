#include <iostream>
#include <fstream>
#include <string>

#include "rism3d.h"
#include "version.h"

bool isNaturalNumber(const std::string& str) {
  try {
    size_t pos;
    int num = std::stoi(str, &pos);
    if (pos != str.length()) {
      return false;
    }
    return num >= 0;
  } catch (const std::exception& e) {
    return false;
  }
}

void RISM3D :: read_input (string control, string structure, bool centering) {

  ifstream in_file;
  in_file.open (control.c_str());

  cout << "reading input data file  : " << control << endl;

  string check;
  in_file >> outlist >> co -> ksave >> check;
  if (check != version) {
    cout << "Error: This input file is for old version." << endl;
    exit (1);
  }

  string closure;
  in_file >> closure;
  if (closure == "KH") {
    clos = 0;
  } else if (closure == "HNC") {
    clos = 1;
  } else {
    cout << "Error: Unexpected closure switch " << endl;
    exit(1);
  }
  in_file >> fsolvent;
  in_file >> co -> convergence >> co -> maxstep;
  in_file >> ma -> count >> ma -> m >> ma -> mp;
  in_file >> ce -> box[0] >> ce -> box[1] >> ce -> box[2];
  in_file >> ce -> grid[0] >> ce -> grid[1] >> ce -> grid[2];

  if (!structure.empty()) {
    in_file.close ();
    in_file.open (structure.c_str());
    cout << "reading solute data file : " << structure << endl;
  }

  int num;

  string tmp;
  in_file >> tmp;

  if (isNaturalNumber(tmp)) {
    num = std::stoi(tmp);
  } else {
    cout << "Error: Number of solute atom is not a natural number." << endl;
    exit(1);
  }

  ce -> setup();
  su -> init(num);

  for (int iu = 0; iu < su -> num; ++iu) {
    int n = iu * 4;
    in_file >> su -> q[iu] >> su -> sig[iu] >> su -> eps[iu]
	    >> su -> r[n] >> su -> r[n + 1]
	    >> su -> r[n + 2];
  }

  if (centering) {
    ce -> shift = su -> centering();
  }

  in_file.close ();

  su -> setup_cuda();
}
