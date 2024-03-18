//  Database of information about rings and their ring-currents and
//    magnetic susceptibilities.

if( rinfo ){   //  read the "ringinfo.in" file to get the needed information

	frinfo = fopen( "ringinfo.in", "r" );
	if( frinfo == NULL ) {
		fprintf( stderr, "empty or non-existent ringinfo.in file\n" );
		exit( 1 );
	}
	while( oline = getline( frinfo ) ) {
		rname = oline;
		ringinfo[rname] = getline( frinfo );
		ringskip[rname] = "";
	}

}else{

//  Database of information about rings and their ring-currents and
//    magnetic susceptibilities.

//  N.B.: Susceptibility tensors for most rings have not yet been calibrated,
//  so "default" ring isotropic and anisotropic portions of chi are used;
//  these are values for benzene, in ppm-A**3/molecule
#define DEFAULT_CHI " -99.2 -91.0"

if( version2 ){
//
//============================================================================
//
//   set up hashed variables to define aromatic rings:
//     [this version has ring-current intensities from Osapay and Case,
//      JACS 113: 9436-9444 (1991)]
//

ringinfo["PHE"] = "CG CD1 CE1 CZ CE2 CD2 1.00" + DEFAULT_CHI;
ringskip["PHE"] = "HA HB2 HB3 HD1 HD2 HE1 HE2 HZ";

ringinfo["TYR"] = "CG CD1 CE1 CZ CE2 CD2 0.84" + DEFAULT_CHI;
ringskip["TYR"] = "HA HB2 HB3 HD1 HD2 HE1 HE2 HH";

ringinfo["HIS"] = "CG CD2 NE2 CE1 ND1 0.90" + DEFAULT_CHI;
ringinfo["HIP"] = "CG CD2 NE2 CE1 ND1 0.90" + DEFAULT_CHI;
ringinfo["HID"] = "CG CD2 NE2 CE1 ND1 0.90" + DEFAULT_CHI;
ringinfo["HIE"] = "CG CD2 NE2 CE1 ND1 0.90" + DEFAULT_CHI;
ringskip["HIS"] = "HA HB2 HB3 HD1 HD2 HE1 HE2";
ringskip["HIP"] = "HA HB2 HB3 HD1 HD2 HE1 HE2";
ringskip["HID"] = "HA HB2 HB3 HD1 HD2 HE1 HE2";
ringskip["HIE"] = "HA HB2 HB3 HD1 HD2 HE1 HE2";

ringinfo["TRP"] = "CG CD2 CE2 NE1 CD1 1.04" + DEFAULT_CHI;
ringinfo["TRP6"] = "CD2 CE2 CZ2 CH2 CZ3 CE3 1.02" + DEFAULT_CHI;
ringskip["TRP"] =  "HA HB2 HB3 HD1 HE1 HZ2 HH2 HZ3 HE3";
ringskip["TRP6"] =  "HA HB2 HB3 HD1 HE1 HZ2 HH2 HZ3 HE3";
ringinfo["WRP"]  = "CG CD2 CE2 NE1 CD1 1.04" + DEFAULT_CHI;
ringinfo["WRP6"] = "CD2 CE2 CZ2 CH2 CZ3 CE3 1.02" + DEFAULT_CHI;
ringskip["WRP"]  =  "HD1 HE1 HZ2 HH2 HZ3 HE3";
ringskip["WRP6"] =  "HD1 HE1 HZ2 HH2 HZ3 HE3";

}else{
//
//============================================================================
//
// set up hashed variables to define aromatic rings:
// ring-current intensities from D.A. Case, J. Biomol. NMR 6: 341-346 (1995):
//
ringinfo["PHE"] = "CG CD1 CE1 CZ CE2 CD2 1.46" + DEFAULT_CHI;
ringskip["PHE"] = "HD1 HD2 HE1 HE2 HZ";

ringinfo["TYR"] = "CG CD1 CE1 CZ CE2 CD2 1.24" + DEFAULT_CHI;
ringskip["TYR"] = "HD1 HD2 HE1 HE2 HH";

ringinfo["HIS"] = "CG CD2 NE2 CE1 ND1 1.35" + DEFAULT_CHI;
ringinfo["HIP"] = "CG CD2 NE2 CE1 ND1 1.35" + DEFAULT_CHI;
ringinfo["HID"] = "CG CD2 NE2 CE1 ND1 1.35" + DEFAULT_CHI;
ringinfo["HIE"] = "CG CD2 NE2 CE1 ND1 1.35" + DEFAULT_CHI;
ringskip["HIS"] = "HD1 HD2 HE1 HE2";
ringskip["HIP"] = "HD1 HD2 HE1 HE2";
ringskip["HID"] = "HD1 HD2 HE1 HE2";
ringskip["HIE"] = "HD1 HD2 HE1 HE2";

ringinfo["TRP"]  = "CG CD2 CE2 NE1 CD1 1.32" + DEFAULT_CHI;
ringinfo["TRP6"] = "CD2 CE2 CZ2 CH2 CZ3 CE3 1.24" + DEFAULT_CHI;
ringskip["TRP"]  =  "HD1 HE1 HZ2 HH2 HZ3 HE3";
ringskip["TRP6"] =  "HD1 HE1 HZ2 HH2 HZ3 HE3";
ringinfo["WRP"]  = "CG CD2 CE2 NE1 CD1 1.32" + DEFAULT_CHI;
ringinfo["WRP6"] = "CD2 CE2 CZ2 CH2 CZ3 CE3 1.24" + DEFAULT_CHI;
ringskip["WRP"]  =  "HD1 HE1 HZ2 HH2 HZ3 HE3";
ringskip["WRP6"] =  "HD1 HE1 HZ2 HH2 HZ3 HE3";

}

ringinfo["NMA"] = "CH3 C O N H CA 0.0 -10.5 -79.3"; // 6-31G** dalton results
// ringskip["NMA"] = "HT1 HT2 HT3 HT4 HT5 HT6 H";
ringskip["NMA"] = "HH31 HH32 HH33 H";

	// folate molecule, using intensity from GUA and TYR
ringinfo["FOL1"] = "N1 C2 N3 C4 C4A C8A 1.00" + DEFAULT_CHI;
ringinfo["FOL2"] = "C4A N5 C6 C7 N8 C8A 0.51" + DEFAULT_CHI;
ringinfo["FOL3"] = "C11 C12 C13 C14 C15 C16 1.24" + DEFAULT_CHI;
ringskip["FOL1"] = "HA21 HA22 H3";
ringskip["FOL2"] = "H7";
ringskip["FOL3"] = "H40 H41 H43 H44";

	// porphyrin ring, using intensity from previous shifts version 
ringinfo["POR1"] = "NA C1A C2A C3A C4A 1.08" + DEFAULT_CHI;
ringinfo["POR2"] = "NB C1B C2B C3B C4B 1.08" + DEFAULT_CHI;
ringinfo["POR3"] = "NC C1C C2C C3C C4C 1.08" + DEFAULT_CHI;
ringinfo["POR4"] = "ND C1D C2D C3D C4D 1.08" + DEFAULT_CHI;
		// do not include ring contribs from inner rings to any of
		// protons in inner rings.
ringskip["POR1"] = "H2A H3A H1N H2B H3B H2C H3C H2N H2D H3D";
ringskip["POR2"] = "H2A H3A H1N H2B H3B H2C H3C H2N H2D H3D";
ringskip["POR3"] = "H2A H3A H1N H2B H3B H2C H3C H2N H2D H3D";
ringskip["POR4"] = "H2A H3A H1N H2B H3B H2C H3C H2N H2D H3D";

	// porph macrocycle, using intensity from previous shifts version
ringinfo["PORM"] = "C4D ND C1D CHD C4C NC C1C CHC C4B NB C1B CHB C4A NA C1A CHA 1.68" + DEFAULT_CHI;
ringskip["PORM"] = "H2A H3A H1N H2B H3B H2C H3C H2N H2D H3D";

// Calibrations from D.A. Case, J. Biomol. NMR 6: 341-346 (1995):

ringinfo["DT"] = "N1 C2 N3 C4 C5 C6 0.35" + DEFAULT_CHI;
ringskip["DT"] = "H6 H71 H72 H73 H3 N3";

ringinfo["U"] = "N1 C2 N3 C4 C5 C6 0.30" + DEFAULT_CHI;
ringskip["U"] = "H6 H5 H3 N3";

ringinfo["C"] = "N1 C2 N3 C4 C5 C6 0.37" + DEFAULT_CHI;
ringskip["C"] = "H5 H6 H41 H42";
ringinfo["DC"] = "N1 C2 N3 C4 C5 C6 0.37" + DEFAULT_CHI;
ringskip["DC"] = "H5 H6 H41 H42";

ringinfo["G"]  = "N9 C4 C5 N7 C8 1.00" + DEFAULT_CHI;
ringinfo["G6"] = "C4 N3 C2 N1 C6 C5 0.51" + DEFAULT_CHI;
ringskip["G"] = "H8 H1 H21 H22 N1";
ringskip["G6"] = "H8 H1 H21 H22 N1";
ringinfo["DG"]  = "N9 C4 C5 N7 C8 1.00" + DEFAULT_CHI;
ringinfo["DG6"] = "C4 N3 C2 N1 C6 C5 0.51" + DEFAULT_CHI;
ringskip["DG"] = "H8 H1 H21 H22 N1";
ringskip["DG6"] = "H8 H1 H21 H22 N1";

ringinfo["A"]  = "N9 C4 C5 N7 C8 1.14" + DEFAULT_CHI;
ringinfo["A6"] = "C4 N3 C2 N1 C6 C5 0.90" + DEFAULT_CHI;
ringskip["A"] = "H2 H8 H61 H62";
ringskip["A6"] = "H2 H8 H61 H62";
ringinfo["DA"]  = "N9 C4 C5 N7 C8 1.14" + DEFAULT_CHI;
ringinfo["DA6"] = "C4 N3 C2 N1 C6 C5 0.90" + DEFAULT_CHI;
ringskip["DA"] = "H2 H8 H61 H62";
ringskip["DA6"] = "H2 H8 H61 H62";

//   exp. susc. values here from Flygare, Chem. Rev. 74:653-687 (1974):
//   ring-current intensity from the Case JBNMR paper cited above:
ringinfo["BEN"] = "C1 C2 C3 C4 C5 C6 1.46 -99.2 -91.0";
ringskip["BEN"] = "H1 H2 H3 H4 H5 H6";

}
