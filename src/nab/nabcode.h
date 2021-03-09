#ifndef	NABCODE_H
#define	NABCODE_H

	/* "stub" types:	*/

#include "defreal.h"
#include "nabc.h"

typedef	char		HASH_T;
typedef	struct	curhash_t	{
	int	index;
	char	*pointer;
} CURHASH_T;

	/* nab builtins (but no libc or libm calls):	*/

#define		I2R(i)	((REAL_T)(i))
INT_T		dumpatom( FILE_T*, RESIDUE_T*, INT_T, INT_T );
INT_T		dumpmolecule( FILE_T*, MOLECULE_T*, INT_T, INT_T, INT_T);
INT_T		dumpresidue( FILE_T*, RESIDUE_T*, INT_T, INT_T );

/*  AmberNetcdf routines: (no argument checking) */
INT_T netcdfDebug();
INT_T netcdfLoad();
INT_T netcdfClose();
INT_T netcdfWriteRestart();
INT_T netcdfCreate();
INT_T netcdfGetVelocity();
INT_T netcdfGetFrame();
INT_T netcdfGetNextFrame();
INT_T netcdfWriteFrame();
INT_T netcdfWriteNextFrame();
INT_T netcdfInfo();

	/* output from molecular mechanics routines goes to nabout  */
FILE_T*  nabout;

	/* defines for hash table references:	*/

typedef	struct	{
	int	v_type;
	union	{
		int	v_ival;
		SIZE_T	v_szval;
		REAL_T	v_fval;
		char	*v_cval;
		REAL_T	v_ptval[ 3 ];
		char	*v_matval;
		FILE	*v_fpval;
		char	*v_atomval;
		char	*v_molval;
		char	*v_resval;
		char	*v_bval;
		char	*v_uval;
	} v_value;
} VALUE_T;

VALUE_T		*NAB_href();
STRING_T	*NAB_hin();
STRING_T	*NAB_hfirst();
STRING_T	*NAB_hnext();

#define	HRI(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_ival
#define	HRSZ(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_szval
#define	HRF(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_fval
#define	HRC(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_cval
#define	HRPT(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_ptval
#define	HRMAT(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_matval
#define	HRFP(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_fpval
#define	HRATOM(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_atomval
#define	HRRES(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_resval
#define	HRMOL(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_molval
#define	HRB(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_bval
#define	HRU(h,k,t,s,c)	(*((c)NAB_href((h),(k),(t),(s))->v_value.v_uval))

	/* defines & declares for points	*/

#define	RECIP(f)	( 1./(f) )

#define	PTX(p)	((p)[0])
#define	PTY(p)	((p)[1])
#define	PTZ(p)	((p)[2])

#define	PTEQ(p,q)	\
	((p)[0]==(q)[0]&&(p)[1]==(q)[1]&&(p)[2]==(q)[2])
#define	PTNE(p,q)	\
	((p)[0]!=(q)[0]||(p)[1]!=(q)[1]||(p)[2]!=(q)[2])
#define	PT_ISTRUE(p)	((p)[0]!=0.0||(p)[1]!=0.0||(p)[2]!=0.0)

POINT_T	*NAB_ptcpy();
POINT_T	*NAB_ptadd();
POINT_T	*NAB_ptsub();
POINT_T	*NAB_ptscl();
POINT_T	*NAB_ptcrs();
REAL_T	NAB_ptdot();

	/* defines for string compares:	*/

int	NAB_strcmp( char *, char * );
#define	LT(a,b)	(NAB_strcmp((a),(b))<0)
#define	LE(a,b)	(NAB_strcmp((a),(b))<=0)
#define	EQ(a,b)	(NAB_strcmp((a),(b))==0)
#define	NE(a,b)	(NAB_strcmp((a),(b))!=0)
#define	GE(a,b)	(NAB_strcmp((a),(b))>=0)
#define	GT(a,b)	(NAB_strcmp((a),(b))>0)

	/* String stuff	*/

#define	length(s)	strlen(s)

char	*NAB_readstring();
char	*NAB_strcpy();
char	*NAB_strcat();
int	NAB_index( char [], char [] );
int	NAB_rematch( char [], char [] );
int NAB_gsub(int, char **, char **, char **);

	/*  Other NAB declarations:  */

int	unlink();
int	NAB_aematch( ATOM_T *ap, char aex[] );

	/* defines for assigning then using temp vars in exprs	*/

#define	ITEMP(t,e)	((t)=(e),&(t))
#define	SZTEMP(t,e)	((t)=(e),&(t))
#define	FTEMP(t,e)	((t)=(e),&(t))
#define	STEMP(t,e)	(NAB_strcpy(&(t),(e)),&(t))
#define	FPTEMP(t,e)	((t)=(e),&(t))
#define	SPRINTF(s,t)	((s),(t))

	/* trig functions in degrees:	*/

#define	R2D	57.29577951308232090712
#define	D2R	 0.01745329251994329576
#define	ACOS(c)	(R2D*acos(c))
#define	ASIN(s)	(R2D*asin(s))
#define	ATAN(t)	(R2D*atan(t))
#define	ATAN2(y,x)	(R2D*atan2((y),(x)))
#define	COS(a)	cos(D2R*(a))
#define	SIN(a)	sin(D2R*(a))
#define	TAN(a)	tan(D2R*(a))

	/* dynamic array allocation macro:	*/

#define DA_ALLOC(a,fn,an)	\
	if(!(a)){fprintf( stderr,	\
	"%s: can't allocate space for dynamic array \"%s\"\n",	\
		(fn),(an));exit(1);}

	/* local string array zero macro:	*/

#define	NAB_AINIT(a)	memset((a),0,sizeof(a))
#define	NAB_SINIT(s)	memset(&(s),0,sizeof(s))

#endif
