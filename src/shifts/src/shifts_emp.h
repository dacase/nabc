//   general indices:
atom		a;
int			i,j; 
point		p;

#define MAXWAT 2000
#define MAXSTRAND 5
#define MAXP 6
#define MAXATBIN 1

int		verbose;           // control amount of printing
int		rinfo;             // determine whether rinfo file is read

string		obsfile, bmrbfile, rdbfile, watfile, outfile;
string		anisfile, banisfile;
string		neighfile;
string		aname, rname;  // for nmr-star output processing
string		map3to1[ hashed ]; //        ditto
string		map1to3[ hashed ]; //        ditto

int			do_CSA;     //  true if CSA tensors are to be computed
int			use_rmd;    //  true if using ring magnetic dipole approx.
int			use_lp;     //  use Sitkoff/Case lone pair charge model
int			version2;   //  always true for now: turns on Osapay/Case
			            //   protein parameters

point		wocoor[ MAXWAT ], wmin, wmax;
int			nwo, woh[ MAXWAT ];

int			obs_found;
float		observed[ hashed ];

int			no_coil[ hashed ];


int			nswap;                  // for swapping pro-chiral pairs
//float		helhn_sw, helha_sw, hfly_sw, hvdw_sw, hvdwc_sw, hvdwh_sw;
float		calc_swap[ hashed ];

file		rdbf, bmrbf, outf;
float		calc_plus_rc[ hashed ];
string		idp;

float		hel_hn, hel_hc, helst, nneigh, aneigh;
int			nbonds;
atom		neighbors[ 8 ];

float		hFly, hFly_aa[ dynamic ], hFlyconv;
int			ipep, npep;
point		pc[ dynamic ], pn[ dynamic ], pbis[ dynamic ]; 
float		pb[ dynamic ];
int			pepres[ dynamic ];

//float		h1p[ dynamic ], h2p1[ dynamic ], h2p2[ dynamic ], 
//			h3p[ dynamic ], h4p[ dynamic ];

point		aa_v[ dynamic ], aa_c[ dynamic ];
int			do_banis, naa, nb, do_aa[ dynamic ], ni, btype[ dynamic ];
int			aa1[ dynamic ], aa2[ dynamic ];
atom		ai[ 10 ];

#define MAXCHI 500
point		chi_v[ MAXCHI ], chi_c[ MAXCHI ];
float		chi_iso[ MAXCHI ], chi_anis[ MAXCHI ];
int			nchi;

int			nh;
point		wato[ MAXP, MAXATBIN ];

float		shhm, she;  // partial and total shifts
float		shp, shb, shb_r, shv, shvc, shvh, sconst;
float		shdist, shhnca, shhnco;
float		sigma[ 3,3 ];    // total shift tensor
float		sigma_p[ 3,3 ];  // result for each peptide bond

float		consthm, constha, consthapg, consthn; 	// constants for proteins
float		consth1p_pyrimidine, consth1p_purine, consth8_guanine;
                                                    // constants for rna,dna

					// total shifts for equivalent protons
float		shhm_t, she_t, shp_t, shb_t, shv_t, shvc_t, shvh_t;
float		shex2_t, shep2_t; 
float		shdist_t, shhnca_t, shhnco_t;
float		hm_t[ hashed ];

					// total  shifts for printing
float		hm[ hashed ];
float		shhm_p[ hashed ], shehn_p[ hashed ], sheha_p[ hashed ];
float       shhm_G6[ hashed ], shhm_G[ hashed ], 
            shhm_A6[ hashed ], shhm_A[ hashed ], 
            shhm_C[ hashed ], shhm_U[ hashed ]; 
float		shp_p[ hashed ];
float		shb_p1[ hashed ];
float		sconst_p[ hashed ];
float		shv_p[ hashed ], shvc_p[ hashed ], shvh_p[ hashed ];
float		shex2_p[ hashed ], shep2_p[ hashed ];
float		shdist_p[ hashed ], shhnca_p[ hashed ], shhnco_p[ hashed ];
string		fields_p[ 4 ];
float		hm_pA[ hashed ], hm_pA6[ hashed ], hm_pG[ hashed ];
float		hm_pG6[ hashed ], hm_pC[ hashed ], hm_pT[ hashed ];
float		hm_pF[ hashed ], hm_pY[ hashed ], hm_pH[ hashed ];
float		hm_pW[ hashed ], hm_pW6[ hashed ], hm_pF1[ hashed ];
float		hm_pF2[ hashed ], hm_pF3[ hashed ];

float		shhm_save[ hashed ], she_save[ hashed ]; // avging aromatics
float		shp_save[ hashed ];
float		shb_save1[ hashed ];
float		shv_save[ hashed ], shvc_save[ hashed ], shvh_save[ hashed ];
float		shdist_save[ hashed ], shhnca_save[ hashed ];
float		shhnco_save[ hashed ];
float		hmF_save[ hashed ], hmY_save[ hashed ], hmH_save[ hashed ];
float		hmW_save[ hashed ], hmW6_save[ hashed ];
float		hmF1_save[ hashed ],hmF2_save[ hashed ],hmF3_save[ hashed ];
string		harom;
int			nprot;    // no. of protons in current equivalent set.

string      ringskip[ hashed ];
