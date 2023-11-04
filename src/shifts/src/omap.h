
string		omap[ hashed ];           //  variables for reading observed

//    set up arrays to help read observed shifts files:

omap["GLY"] = "H HA2 HA3";
omap["ALA"] = "H HA HB3";
omap["VAL"] = "H HA HB HG13 HG23";
omap["LEU"] = "H HA HB2 HB3 HG HD13 HD23";
omap["ILE"] = "H HA HB HG12 HG13 HG23 HD13";
omap["SER"] = "H HA HB2 HB3";
omap["THR"] = "H HA HB HG23";
omap["PHE"] = "H HA HB2 HB3 HD2 HE2 HZ";
omap["HIP"] = "H HA HB2 HB3 HE1 HD2";
omap["HIS"] = "H HA HB2 HB3 HE1 HD2";
omap["HID"] = "H HA HB2 HB3 HE1 HD2";
omap["HIE"] = "H HA HB2 HB3 HE1 HD2";
omap["TYR"] = "H HA HB2 HB3 HD2 HE2 HH";
omap["TRP"] = "H HA HB2 HB3 HE1 HD1 HE3 HZ2 HH2 HZ3";
omap["PRO"] = "HA HB2 HB3 HG2 HG3 HD2 HD3";
omap["ASP"] = "H HA HB2 HB3";
omap["GLU"] = "H HA HB2 HB3 HG2 HG3";
omap["ASN"] = "H HA HB2 HB3 HD21 HD22";
omap["GLN"] = "H HA HB2 HB3 HG2 HG3 HE21 HE22";
omap["ARG"] = "H HA HB2 HB3 HG2 HG3 HD2 HD3";
omap["LYS"] = "H HA HB2 HB3 HG2 HG3 HD2 HD3 HE2 HE3";
omap["MET"] = "H HA HB2 HB3 HG2 HG3 HE3";
omap["CYS"] = "H HA HB2 HB3 HG";
omap["CYX"] = "H HA HB2 HB3";
omap["MTH"] = "H1 H2 H3 H4";
omap["PRP"] = "H11 H12 H13 H21 H22 H31 H32 H33";
omap["ETH"] = "H11 H12 H13 H21 H22 H23";
omap["IBU"] = "H11 H12 H13 H2 H31 H32 H33 H41 H42 H43";
omap["IPE"] = "H11 H12 H13 H21 H31 H32 H33 H41 H42 H51 H52 H53";
omap["NEP"] = "H11 H12 H13 H31 H32 H33 H41 H42 H43 H51 H52 H53";
omap["NIA"] = "H";
omap["IAM"] = "HNT";
omap["IAC"] = "HNT";
omap["POR"] = "H2A H3A H2B H3B H2C H3C H2D H3D H1N H2N H2PA H3PA H5PA H6PA H73A H2PB H3PB H5PB H6PB H73B H2PC H3PC H5PC H6PC H73C H2PD H3PD H5PD H6PD H73D ";
omap["FOL"] = "HA21 HA22 H3 H7 H92 H93 H10 H40 H41 H43 H44 H HA HB1 HB2 HG1 HG2";

omap["DC"] = "H1' H2' H2'' H3' H4' H5' H5'' H41 H42 H5 H6";
omap["DT"] = "H1' H2' H2'' H3' H4' H5' H5'' H3 H6 H73";
omap["DA"] = "H1' H2' H2'' H3' H4' H5' H5'' H2 H61 H62 H8";
omap["DG"] = "H1' H2' H2'' H3' H4' H5' H5'' H1 H21 H22 H8";

omap["C"] = "H1' H2' HO'2 H3' H4' H5' H5'' H41 H42 H5 H6";
omap["A"] = "H1' H2' HO'2 H3' H4' H5' H5'' H2 H61 H62 H8";
omap["G"] = "H1' H2' HO'2 H3' H4' H5' H5'' H1 H21 H22 H8";
omap["U"] = "H1' H2' HO'2 H3' H4' H5' H5'' H3 H5 H6";

omap["PSU"] = "H1' H2'1 HO'22 H3' H4' H5'1 H5'2 H1 H3 H6";
omap["H2U"] = "H1' H2'1 HO'2 H3' H4' H5'1 H5'2 H3 H51 H52 H61 H62";
omap["7MG"] = "H1' H2'1 HO'2 H3' H4' H5'1 H5'2 H1 H21 H22 H73 H8";
omap["1MA"] = "H1' H2'1 HO'2 H3' H4' H5'1 H5'2 H11 H12 H13 H2 H6 H8";
omap["2MG"] = "H1' H2'1 HO'2 H3' H4' H5'1 H5'2 H1 H2 HA3 H8";
omap["M2G"] = "H1' H2'1 HO'2 H3' H4' H5'1 H5'2 H1 H8 HB3 HB3";
omap["5MC"] = "H1' H2'1 HO'2 H3' H4' H5'1 H5'2 H41 H42 H6 H73";
omap["5MU"] = "H1' H2'1 HO'2 H3' H4' H5'1 H5'2 H3 H6 H73";
omap["OMG"] = "H1' H2'1 HA'3 H3' H4' H5'1 H5'2 H1 H21 H22 H8";
omap["OMC"] = "H1' H2'1 HA'3 H3' H4' H5'1 H5'2 H41 H42 H5 H6";
omap["YG"]  = "H1' H2'1 HO'2 H3' H4' H5'1 H5'2 H33 H2 H8 H103 H133 H141 H142 H15 H193 H20 H243";
omap["MTH"] = "H1 H2 H3 H4";
omap["BEN"] = "H1 H2 H3 H4 H5 H6";
omap["PRT"] = "HP";
omap["NRG"] = "H HA HB2 HB3 HG2 HG3 HD2 HD3";
omap["CLA"] = "H HA HB3";
omap["NME"] = "H";
omap["NH2"] = "H H1";
