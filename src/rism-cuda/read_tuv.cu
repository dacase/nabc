#include <iostream>
#include <fstream>
#include "rism3d.h"
#include "extension.h"

void RISM3D :: read_tuv () {

  ifstream in_file;
  in_file.open((fname + exttuv).c_str());

  int ng0, ng1, ng2;
  int natvo;
  in_file >> ng0 >> ng1 >> ng2 >> natvo;
  if (ng0 != ce -> grid[0] || ng1 != ce -> grid[1] || 
      ng2 != ce -> grid[2] || natvo != sv -> natv) {
    cout << " RXRISM:  actual   NGr = " << ce -> ngrid
	 << "   NatV = " << sv -> natv << endl 
	 << "           saved   NGr = " << ng0 * ng1 * ng2 
	 << "   NatV = " << natvo << endl; 
  }

  for (int iv = 0; iv < sv -> natv; ++iv) {
    for (int ig = 0; ig < ce -> ngrid; ++ig) {
      in_file >> tuv[iv][ig];
    }
  }

  in_file.close();

  for (int iv = 0; iv < sv -> natv; ++iv) {
    cudaMemcpyAsync(dt + (iv * ce -> ngrid), tuv[iv],
		    ce -> ngrid * sizeof(double), cudaMemcpyDefault);
  }
}
