#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabcode.h"
extern char NAB_rsbuf[];
static int mytaskid, numtasks;





























 struct xmin_opt{INT_T mol_struct_opt;INT_T maxiter;REAL_T grms_tol;INT_T method;INT_T numdiff;INT_T m_lbfgs;INT_T iter;REAL_T xmin_time;INT_T ls_method;INT_T ls_maxiter;REAL_T ls_maxatmov;REAL_T beta_armijo;REAL_T c_armijo;REAL_T mu_armijo;REAL_T ftol_wolfe;REAL_T gtol_wolfe;INT_T ls_iter;INT_T print_level;INT_T error_flag;};
































 struct lmod_opt{INT_T niter;INT_T nmod;INT_T kmod;INT_T nrotran_dof;INT_T nconf;REAL_T minim_grms;REAL_T energy_window;REAL_T conf_separation_rms;INT_T eig_recalc;INT_T ndim_arnoldi;INT_T lmod_restart;INT_T n_best_struct;INT_T mc_option;REAL_T rtemp;REAL_T lmod_step_size_min;REAL_T lmod_step_size_max;INT_T nof_lmod_steps;REAL_T lmod_relax_grms;INT_T nlig;INT_T apply_rigdock;INT_T nof_poses_to_try;INT_T random_seed;INT_T print_level;REAL_T lmod_time;REAL_T aux_time;INT_T error_flag;};




static  struct xmin_opt xo;

static  struct lmod_opt lo;


static INT_T i, j, k, natm, all_frames, n_frames;

static MOLECULE_T *mol;

static ATOM_T *ai;

static REAL_T ene, grms, glob_min_energy;

static INT_T *lig_start,  *lig_end,  *lig_cent;


static REAL_T *x,  *g,  *x_ref,  *conflib,  *lmod_traj;


static REAL_T *tr_min,  *tr_max,  *rot_min,  *rot_max;

static STRING_T *fname = NULL,  *sys_command = NULL;

static REAL_T dx, dy, dz;

static POINT_T dummy;


int main( argc, argv )
	int	argc;
	char	*argv[];
{
	nabout = stdout; /*default*/

	if ( mpiinit(&argc, argv, &mytaskid, &numtasks) != 0 ) {
		printf("Error in mpiinit!");
		fflush(stdout);
		exit(1);
	}
static INT_T __gdab0001__;
static INT_T __gdab0002__;
static INT_T __gdab0003__;
static INT_T __gdab0004__;
static INT_T __gdab0005__;
static INT_T __gdab0006__;
static INT_T __gdab0007__;
static INT_T __gdab0008__;
static INT_T __gdab0009__;
static INT_T __gdab0010__;
static INT_T __gdab0011__;
static INT_T __gdab0012__;
static INT_T __it0001__;
static STRING_T *__st0001__ = NULL;
if( argc == 1 ){

fprintf( nabout, "Usage:  orterun -np n_cores  %s  input_pdb_file  input_prmtop_file  output_file_base_name  ", argv[1 - 1] );
fprintf( nabout, "number_of_vib_modes (max 200)  number_of_frames_per_half_swing (max 50) > log_file 2>&1 &\n" );
exit( 0 );
}

if( ( atoi( argv[5 - 1] ) ) > 200 || ( atoi( argv[6 - 1] ) ) > 50 ){fprintf( nabout, "Input error!\n" );exit( 1 );}

mol = getpdb( argv[2 - 1], NULL );
readparm( mol, argv[3 - 1] );
natm =  *( NAB_mri( mol, "natoms" ) );

__gdab0004__ = 3 * natm;DA_ALLOC( x = ( REAL_T * )malloc( __gdab0004__ * ( sizeof( REAL_T ) ) ), "main", "x" );__gdab0005__ = 3 * natm;DA_ALLOC( g = ( REAL_T * )malloc( __gdab0005__ * ( sizeof( REAL_T ) ) ), "main", "g" );__gdab0006__ = 3 * natm;DA_ALLOC( x_ref = ( REAL_T * )malloc( __gdab0006__ * ( sizeof( REAL_T ) ) ), "main", "x_ref" );

setxyz_from_mol(  &mol, NULL, x );



xmin_opt_init(  &xo );
lmod_opt_init(  &lo,  &xo );

lo . niter = 1;
lo . nconf = atoi( argv[5 - 1] );
lo . minim_grms = 1.000000E+00;
lo . nmod = atoi( argv[5 - 1] );
lo . kmod = atoi( argv[5 - 1] );
lo . nrotran_dof = 6;
lo . eig_recalc = 1;

lo . ndim_arnoldi = 50;
lo . lmod_step_size_min = 3.000000E-01;
lo . lmod_step_size_max = 3.000000E-01;
lo . nof_lmod_steps = atoi( argv[6 - 1] );
lo . lmod_relax_grms = 1.000000E+01;
lo . print_level = 5;

xo . maxiter = 5000;
xo . grms_tol = 1.000000E-12;
xo . method = 3;
xo . m_lbfgs = 5;
xo . print_level = 1;
lo . random_seed = 0;

all_frames = lo . niter * ( lo . kmod + 1 ) * ( 4 * lo . nof_lmod_steps + 3 );

__gdab0007__ = ( lo . kmod + 1 ) * 3 * natm;DA_ALLOC( conflib = ( REAL_T * )malloc( __gdab0007__ * ( sizeof( REAL_T ) ) ), "main", "conflib" );
__gdab0008__ = all_frames * 3 * natm;DA_ALLOC( lmod_traj = ( REAL_T * )malloc( __gdab0008__ * ( sizeof( REAL_T ) ) ), "main", "lmod_traj" );

setmol_from_xyz(  &mol, NULL, x );
setxyz_from_mol(  &mol, NULL, x_ref );

mm_options( "ntpr=1,gb=0,kappa=0.0,cut=999.9,diel=R" );
mme_init( mol, NULL, "::ZZZ", dummy, NULL );
mme( x, g, ITEMP( __it0001__, 1 ) );

if( mytaskid == 0 ){

NAB_strcpy(  &sys_command, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "/bin/rm -f conflib.dat" ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
system( sys_command );
}



glob_min_energy = lmod(  &natm, x, g,  &ene, conflib, lmod_traj, lig_start, lig_end, lig_cent, tr_min, tr_max, rot_min, rot_max,  &xo,  &lo );



for( i = 1;i <= all_frames;i = i + 1 )
{
setmol_from_xyz(  &mol, NULL,  &lmod_traj[( i - 1 ) * 3 * natm + 1 - 1] );
for( ai = NULL;ai = NAB_mnext( mol, ai ); ){
k = ( i - 1 ) / ( 4 * lo . nof_lmod_steps + 3 );
j =  *( NAB_ari( ai, "tatomnum" ) ) - 1;
dx = conflib[k * 3 * natm + 3 * j + 1 - 1];
dy = conflib[k * 3 * natm + 3 * j + 2 - 1];
dz = conflib[k * 3 * natm + 3 * j + 3 - 1];
ai->a_bfact = 1000 * ( sqrt( dx * dx + dy * dy + dz * dz ) );

}
NAB_strcpy(  &fname, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "%s-frame_%04d.pdb", argv[4 - 1], i ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
putpdb( fname, mol, NULL );
}

NAB_strcpy(  &sys_command, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "sleep 60" ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );




if( mytaskid == 0 ){


n_frames = 4 * lo . nof_lmod_steps + 3;
k = 1;

for( j = 1;j <= lo . kmod + 1;j = j + 1 )
{
NAB_strcpy(  &sys_command, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "/bin/rm -f %s-vib-mode_%04d.pdb; touch %s-vib-mode_%04d.pdb", argv[4 - 1], j, argv[4 - 1], j ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
system( sys_command );

for( i = 1;i <= n_frames;i = i + 1 )
{
NAB_strcpy(  &sys_command, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "sed -i '1s/^/MODEL%8d\\n/' %s-frame_%04d.pdb", i, argv[4 - 1], k ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
system( sys_command );
NAB_strcpy(  &sys_command, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "cat '%s-frame_%04d.pdb' >> %s-vib-mode_%04d.pdb", argv[4 - 1], k, argv[4 - 1], j ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
system( sys_command );
NAB_strcpy(  &sys_command, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "echo 'ENDMDL' >> %s-vib-mode_%04d.pdb", argv[4 - 1], j ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
system( sys_command );

k = k + 1;
}
}

NAB_strcpy(  &sys_command, SPRINTF( assert( ( snprintf( NAB_rsbuf, 10000, "/bin/rm -f %s-frame_????.pdb conflib.dat", argv[4 - 1] ) ) < 10000 ), NAB_strcpy(  &__st0001__, NAB_rsbuf ) ) );
system( sys_command );
}


	mpifinalize();

	exit( 0 );
}
