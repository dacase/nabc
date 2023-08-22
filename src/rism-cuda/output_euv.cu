#include <iostream>
#include <fstream>
#include "rism3d.h"
#include "extension.h"
  
void RISM3D :: output_euv (double euv) {
     
  ofstream out_file;
  out_file.open ((fname + exteuv).c_str());
  
  double ibeta = avogadoro * boltzmann * sv -> temper;
  
  out_file << "Soulte_Solvent_Interaction_Energy= " << ibeta * euv
	   << " (J/mol)" << endl ;
  
  out_file.close () ;
}
