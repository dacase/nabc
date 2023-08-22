#include <iostream>
#include <fstream>
#include <stdio.h>
#include "rism3d.h"
#include "extension.h"

void RISM3D :: output_huv() {

  cout << "outputting Huv to file:  " << fname + exthuv << "  ..." << endl;

  ofstream out_file;
  out_file.open((fname + exthuv).c_str());

  out_file << ce -> box[0] << " " 
	   << ce -> box[1] << " " 
	   << ce -> box[2] << " " 
	   << endl
	   << ce -> grid[0] << " " 
	   << ce -> grid[1] << " " 
	   << ce -> grid[2] << " " 
	   << endl 
	   << su -> num << endl;

  for (int iu = 0; iu < su -> num; ++iu) {
    int num = iu * 3;
    out_file << su -> q[iu] << " "
	     << su -> sig[iu] << " " 
	     << su -> eps[iu] << " "
	     << su -> r[num] << " "
	     << su -> r[num + 1] << " "
	     << su -> r[num + 2] << endl;
  }
  out_file << sv -> natv << endl;

  for (int ig = 0; ig < ce -> ngrid; ++ig) {
    for (int iv = 0; iv < sv -> natv; ++iv) {
      out_file << huv[iv][ig].real() << " ";
    }
    out_file << endl;
  }

  cout << "done." << endl;

  out_file.close();
} 
