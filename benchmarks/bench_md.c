#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "sff.h"

#ifdef MPI
   int   mpierror( int );
   int   mpifinalize( void );
   int   mpiinit( int*, char**, int*, int* );
   int mytaskid;
   int numtasks;
#else
   int mytaskid = 0;
   int numtasks = 1;
#endif

FILE* nabout;

int main( int argc, char *argv[] )
{

#ifdef MPI
   if ( mpiinit(&argc, argv, &mytaskid, &numtasks) != 0 ) {
      printf("Error in mpiinit!");
      fflush(stdout);
      exit(1);
   }
#endif

//   read in the prmtop file and the coordinates:

   PARMSTRUCT_T* prm = rdparm( argv[1] );    // reads the prmtop file
   int natm = prm->Natom;
   int natm3 = 3 * natm;
   double* x = malloc( natm3 * (sizeof(double)) );
   double* f = malloc( natm3 * (sizeof(double)) );
   double* v = malloc( natm3 * (sizeof(double)) );

   //   option: read a formatted restart file:
   float start_time;
   getxv( argv[2], natm, start_time, x, v );

#if 0
   //  (next lines could be bundled into a netcdfReadRestart() function)
   struct AmberNetcdf ain;
   int ier = netcdfLoad( &ain, argv[2] );
   // netcdfDebug( &ain );
   // netcdfInfo( &ain );
   netcdfGetFrame(  &ain, 0, x, NULL );
   netcdfGetVelocity( &ain, 0, v );
#endif

// Initialize molecular mechanics..

mm_options("cut=20.0, rgbmax=20.0, ntpr=100, nsnb=10, gb=1, diel=C, rattle=1");

   // nothing frozen or constrained for now:
   int frozen[natm], constrained[natm];
   memset( frozen, 0, natm * sizeof(int) );
   memset( constrained, 0, natm * sizeof(int) );

   mme_init_sff( prm, frozen, constrained, NULL, NULL );


// Do some molecular dynamics.

   mm_options("tautp=0.4, temp0=100.0, ntpr_md=10, tempi=50.0, ntpr=10000");
   md(natm3, 100, x, f, v, mme);
   if( mytaskid == 0 )
      netcdfWriteRestart( "xfin.md2.x", natm, x, v, NULL, 1.0, 0.0 );

#ifdef MPI
   mpifinalize();
#endif

}

