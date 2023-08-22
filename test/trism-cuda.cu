#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "sff.h"
#include "/home/case/rism/3D-RISM-CUDA/src-cuda/rism3d.h"
FILE* nabout;

//  Very primitive short rism-cuda minimizer, no argument-
//     checking, etc., etc.
//
//  Usage:  trism-cuda  <parm-file>  <restart-file>
//       output to stdout

void init_rism() {

   RISM3D *system;
   int dn=0;
   cudaSetDevice(dn);
   fprintf( stderr, "back from cudaSetDevice\n" );
   system = new RISM3D;
   system -> initialize( "trpcage_c", "trpcage_s", 0 );
   fprintf( stderr, "back from system->initialize\n" );

}

int main( int argc, char *argv[] )
{
   RISM3D *system;      //  main struct for rism
   PARMSTRUCT_T *prm;   //  struct to hold info from a prmtop file
   XMIN_OPT_T xo;       //  options for the minimizer
   double *xyz,  *grad;
   double energy, grms;
   double start_time = 0.0;   // dummy, since this is minimization, not md

   nabout = stdout;    // change to redirect output (historical kludge)

//   options for the minimizer:

   xmin_opt_init( &xo );  // sets the default parameters;

   xo.maxiter = 5;                  // non-default minimization options:
   xo.grms_tol = 0.0005;
   xo.ls_maxatmov = 0.15;
   xo.print_level = 1;
   xo.method = 2;

//   read in the prmtop file and the coordinates:

   prm = rdparm( argv[1] );    // reads the prmtop file
   int natm = prm->Natom;
   xyz = (double *) malloc( 3 * natm * (sizeof(double)) );
   grad = (double *) malloc( 3 * natm * (sizeof(double)) );
   getxv( argv[2], natm, start_time, xyz, grad );  // reads a restart file

//   setup the force field parameters, and get an initial energy:

   mm_options( "ntpr=1, cut=99.0, diel=C " );
   mm_options( "ntpr=1, gb=3, cut=9999.0" );
   // mm_options( "xvvfile=../rism1d/spc-kh/spc.xvv.save" );
   // mm_options( "verbose=0" );
   // mm_options( "buffer=-1, ng=30,30,30, solvbox=15,15,15, solvcut=999" );
   // mm_options( "tolerance=1e-7" );
   // mm_options( "apply_rism_force=1, centering=2" );
   // mm_options( "ntpr_rism=1, ntwrism=0" );
   // mm_options( "uccoeff=-0.149818,-0.1136266,-0.00053163,0.0103954" );

   // nothing frozen or constrained for now:
   int* frozen = parseMaskString( "@ZZZ", prm, xyz, 2 );
   int* constrained = parseMaskString( "@ZZZ", prm, xyz, 2 );

   mme_init_sff( prm, frozen, constrained, NULL, NULL );
   fprintf( stderr, "ready for init_rism\n");
   init_rism();
   fprintf( stderr, "back from init_rism\n");
   exit(0);   // just for initial testing

   int verbose = -1;   // historical flag to give more verbose output
   energy = mme( xyz, grad, &verbose );

//   run the minimization:

   char title[] = " non-rattle minimization";
   energy = xmin( mme,  &natm, xyz, grad,  &energy,  &grms,  &xo );
// putxv( argv[3], title, natm, start_time, xyz, xyz );
   energy = mme( xyz, grad, &verbose );

//  optional: get solvent distribution at the end:
//   mm_options( "guvfile=g.xmin, apply_rism_force=0, ntwrism=1" );
//   mme_init_sff( prm, frozen, constrained, NULL, NULL );
//   mme( xyz, grad, &verbose );

}

