#include <iostream>
#include <fstream>
#include <string>
#include "rism3d.h"
#include "extension.h"

void RISM3D :: write_tuv () {

  cout << "saving Tuv to file:  " << fname + exttuv << "  ..." << endl;

  for (int iv = 0; iv < sv -> natv; ++iv) {
    cudaMemcpyAsync(tuv[iv], dt + (iv * ce -> ngrid), 
		    ce -> ngrid * sizeof(double), cudaMemcpyDefault);
  }

  ofstream out_file;
  out_file.open((fname + exttuv).c_str());

  out_file << ce -> grid[0] << " " 
	   << ce -> grid[1] << " " 
	   << ce -> grid[2] << " " 
	   << sv -> natv << endl;

  for (int iv = 0; iv < sv -> natv; ++iv) {
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
      out_file << tuv[iv][ig] << endl;
    }
  }

  cout << "done." << endl;

  out_file.close();
}
