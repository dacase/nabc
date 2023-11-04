//
//============================================================================
//
//    Values in this file provide the parameterization used in the
//      empirical calculations.  (Parameters for aromatic rings are
//      given in ringinfo.h)
//
//    References: Osapay & Case, JACS 113: 9436-9444 (1991).
//                Sitkoff & Case, JACS 119: 12262-12272 (1997).
//                Dejaegere, Bryce & Case, An empirical analysis of
//                   proton chemicalshifts in nucleic acids.  In "Modeling
//                   NMR Chemical Shifts", A.C. de Dios and J. Facelli, eds.
//                   (Washington, American Chemical Society, 1999), in press.

//============================================================================
//
//    Parameters for the peptide group:

//   (1) anisotropy parameters, units are A^3

hFly = -7.8828 * 1.66058;  // suscept anis for peptide grps (from Osap & Case)

//    next two lines are for a non-axial peptide anisotropy model, discussed
//    in Sitkoff & Case.
//	hFlychi1 = 3.653;	// suscept anis prm1 assuming nonsymm pept grp (Flygare)
//	hFlychi2 = 13.285;	// suscept anis prm2 assuming nonsymm pept grp (Flygare)

hFlyconv = 1.66058;	// converts Flygare's units to A**3 (for bond anis)

//    (2) electrostatic parameters, units are 10**-12 esu**-1

hel_hn = -1.2 * 4.803;	// value from Osapay & Case paper
hel_hc = -1.2 * 4.803;	// value from Osapay & Case paper

//helst    = -14.93;	// Buckingham A value for HM fit to nucleic acid
                    	//    rings, from Table 3 of JBNMR paper

//============================================================================
//
//    Constant contribution to some protein shifts, from Osapay & Case:

consthm = -0.041;	// constant for aliphatic protons, default
constha = -0.754;	// constant for HA
consthapg = -0.51;	// constant for HA from gly, pro
consthn = -0.550;	// constant for HN

//    Constants for nucleic acids, from Dejaegere, Bryce & Case:

consth1p_pyrimidine = 0.32;
consth1p_purine = -0.20;
consth8_guanine = 0.34;

#if 0
//============================================================================
//
//   van der Waals prms, from Sitkoff & Case, units are ppmA**3/eV
//   (N.B.: these are not used for "standard" calculations (yet).  They
//    are only referenced in the code is compiled with -DVDW set.)

//	hvdw_def = 6.79;  	// default van der Waals parameter
//  (The default  vdw prm is actually: 
//	2.81E-18 esu * [ 3/2 * 0.1113E-8 / 183.63 / 0.03779E-28 ppmA**3/eV/esu ] )

//	hvdw_c   = 10.54; 	// for H near C
//	hvdw_h   = 0.0;  	// for H near H
pol["H"] = 0.667; 
pol["O"] = 0.802;
pol["N"] = 1.1;
pol["C"] = 1.76;
pol["S"] = 2.90;
pol["F"] = 0.557;
pol["c"] = 2.18;		// Chlorine
pol["h"] = 0.205;		// Helium
pol["n"] = 0.386;		// Neon
pol["A"]  = 1.64; 		// Argon
ion_en["H"] = 13.6;
ion_en["O"] = 13.6;
ion_en["N"] = 14.5;
ion_en["C"] = 11.3;
ion_en["S"] = 10.3;
ion_en["F"] = 17.4;
ion_en["c"] = 12.97;		// Chlorine
ion_en["h"] = 24.59;		// Helium
ion_en["n"] = 21.56;		// Neon
ion_en["A"]  = 15.76;		// Argon
#endif

//============================================================================
//
// Here are more prms from Sitkoff & Case analyses of gas phase DFT shifts:
//    these are only used if special items are turned on during compilation.

#ifdef PANIS_SC
	hFly = 6.25 * 1.66058;	// peptide anisotropy
#endif
#ifdef EL_SC
	hel_hn = - 3.46 * 4.803;
	hel_hc = - 3.81 * 4.803;
#endif

