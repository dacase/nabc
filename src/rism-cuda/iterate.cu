#include <iostream>
#include <fstream>
#include "rism3d.h"
#include "extension.h"

void RISM3D :: iterate(int cu) {
  void alloc2D (vector <double *> &, int, int);
  void calloc2D (vector <complex <double> *> &, int, int);

  double cf, cuf;

  calloc2D (guv, sv -> natv, ce -> ngrid);
  calloc2D (huv, sv -> natv, ce -> ngrid);
  alloc2D (tuv, sv -> natv, ce -> ngrid);

  cudaMalloc(&dguv, ce -> ngrid * sv -> natv * sizeof(double2));
  cudaMalloc(&dhuv, ce -> ngrid * sv -> natv * sizeof(double2));
  cudaMalloc(&dt, ce -> ngrid * sv -> natv * sizeof(double));
  cudaMalloc(&dtr, ce -> ngrid * sv -> natv * sizeof(double));
  cudaMalloc(&ds, ce -> grid[1] * ce -> grid[2] * sizeof(double));

  ma -> initialize (ce, sv);
  fft -> initialize (ce);

  ifstream in_file ;
  in_file.open((fname + exttuv).c_str());
  bool saved = in_file.is_open();
  in_file.close();

  if (saved) {
    read_tuv();
    cu = 0;
    cf = 1.0;
  } else if (cu == 0) {
    initialize_tuv();
    cf = 1.0;
  } else {
    cudaMemset(dt, 0.0, ce -> ngrid * sv -> natv * sizeof(double));
    cuf = 1.0 /cu;
    cf = 0.0;
  }
  for (int c = 0; c <= cu; ++c) {
    if (c >0) {
      cf += cuf;
      if (cf > 1.0) cf = 1.0;
      add_tuv(cuf);
    }
    cout << "relaxing 3D UV RISM: Charge Up Factor = " << cf << endl;
    bool conver = false;
    bool diverge = false;
    for (int istep = 1; istep <= co -> maxstep; ++istep) {
      calculate(cf);
      double rms = cal_rms ();
      diverge = !isfinite(rms);
      if (diverge) {
	break;
      }
      if (rms <= co -> convergence) {
	conver = true;
      } else {
	ma -> calculate (dt, dtr);
      }
      if( istep%50 == 0 )
         cout << " Step = " << istep << " Reside = " << rms << endl;
      if (co -> ksave > 0 && istep % co -> ksave == 0) {
	write_tuv();
      }
      if (conver) {
	if (co -> ksave != 0 && c == cu) {
	  write_tuv();
	}
	break;
      }
    }
    if (diverge) {
      cout << "Calculation diverged." << endl;
    } else if (!conver) {
      cout << "3D UV RISM: reached limit # of relaxation steps: "
	   << co -> maxstep << endl;
      break;
    }
  }
  for (int iv = 0; iv < sv -> natv; ++iv) {
    cudaMemcpyAsync(huv[iv], dhuv + (iv * ce -> ngrid), 
	       ce -> ngrid * sizeof(double2), cudaMemcpyDefault);
    cudaMemcpyAsync(guv[iv], dguv + (iv * ce -> ngrid), 
	       ce -> ngrid * sizeof(double2), cudaMemcpyDefault);
  }
  delete ma;
  delete fft;
} 
