/// ptrajmask: The enhanced atom selection mask parser from ptraj.
/// Originally written by Viktor Hornak, Stony Brook University.
/// Adapted as standalone code by Dan Roe, NIST.
#ifdef __cplusplus
extern "C" {
#endif
// The NAME parameters control how the mask parser expects strings to look.
// Originally, the parameter file assumes the atom, residue, symbol, etc. 
// names to be four characters. When stored as a string, the NULL character
// is required, requiring a size of 5. It has been increased to 6 in cpptraj
// to accomodate the slightly larger names that can be found in MOL2 files.
// NAMESIZE: Default size for atom and residue names, 5 + NULL.
// Amber atom/residue names are 4, but some mol2 atom types are larger.
#define NAMESIZE 6
#define NAME_DEFAULT "     "
// Default type for atom names, res names, atom types etc
typedef char NAME[NAMESIZE];
// More defines
#define  MAXSELE 1000
#define  ALL      0
#define  NUMLIST  1
#define  NAMELIST 2
#define  TYPELIST 3
#define  ELEMLIST 4
/* parseMaskString()
 * The main interface to the mask parser. Takes a mask expression and some
 * information from a parameter file (# atoms, # residues, atom names, residue
 * names, an array containing the first atom # of each residue, atomic coords
 * in X0Y0Z0X1Y1Z1... format, atom types, and a debug value controlling how
 * much debug information is printed (internally the global int prnlev).
 * It returns a character mask array mask[i]='T'|'F', i=0,atoms-1
 * which contains the resulting atom selection
 */
int *parseMaskString(char*,int,int,NAME*,NAME*,int*,double*,NAME*,int);
#ifdef __cplusplus
}
#endif
