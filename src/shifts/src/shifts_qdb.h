
// quantum chemistry database parameter declarations:

float A_s_N_bb[112,21], A_s_Ca_bb[112,21], A_s_Cb_bb[112,21], A_s_CO_bb[112,21];
float G_s_N_bb[112,21], G_s_Ca_bb[112,21], G_s_Cb_bb[112,21], G_s_CO_bb[112,21];
float Q_s_N_bb[112,21], Q_s_Ca_bb[112,21], Q_s_Cb_bb[112,21], Q_s_CO_bb[112,21];
float S_s_N_bb[112,21], S_s_Ca_bb[112,21], S_s_Cb_bb[112,21], S_s_CO_bb[112,21];
float V_s_N_bb[112,21], V_s_Ca_bb[112,21], V_s_Cb_bb[112,21], V_s_CO_bb[112,21];
float F_s_N_bb[112,21], F_s_Ca_bb[112,21], F_s_Cb_bb[112,21], F_s_CO_bb[112,21];
float L_s_N_bb[112,21], L_s_Ca_bb[112,21], L_s_Cb_bb[112,21], L_s_CO_bb[112,21];

float A_h_N_bb[112,21], A_h_Ca_bb[112,21], A_h_Cb_bb[112,21], A_h_CO_bb[112,21];
float G_h_N_bb[112,21], G_h_Ca_bb[112,21], G_h_Cb_bb[112,21], G_h_CO_bb[112,21];
float Q_h_N_bb[112,21], Q_h_Ca_bb[112,21], Q_h_Cb_bb[112,21], Q_h_CO_bb[112,21];
float S_h_N_bb[112,21], S_h_Ca_bb[112,21], S_h_Cb_bb[112,21], S_h_CO_bb[112,21];
float V_h_N_bb[112,21], V_h_Ca_bb[112,21], V_h_Cb_bb[112,21], V_h_CO_bb[112,21];
float F_h_N_bb[112,21], F_h_Ca_bb[112,21], F_h_Cb_bb[112,21], F_h_CO_bb[112,21];
float L_h_N_bb[112,21], L_h_Ca_bb[112,21], L_h_Cb_bb[112,21], L_h_CO_bb[112,21];

float A_n_N_bb[112,21], A_n_Ca_bb[112,21], A_n_Cb_bb[112,21], A_n_CO_bb[112,21];
float G_n_N_bb[112,21], G_n_Ca_bb[112,21], G_n_Cb_bb[112,21], G_n_CO_bb[112,21];
float Q_n_N_bb[112,21], Q_n_Ca_bb[112,21], Q_n_Cb_bb[112,21], Q_n_CO_bb[112,21];
float S_n_N_bb[112,21], S_n_Ca_bb[112,21], S_n_Cb_bb[112,21], S_n_CO_bb[112,21];
float V_n_N_bb[112,21], V_n_Ca_bb[112,21], V_n_Cb_bb[112,21], V_n_CO_bb[112,21];
float F_n_N_bb[112,21], F_n_Ca_bb[112,21], F_n_Cb_bb[112,21], F_n_CO_bb[112,21];
float L_n_N_bb[112,21], L_n_Ca_bb[112,21], L_n_Cb_bb[112,21], L_n_CO_bb[112,21];

float A_Nchi_p[10,21], G_Nchi_p[10,21], Q_Nchi_p[10,21], S_Nchi_p[10,21], V_Nchi_p[10,21],
      L_Nchi_p[10,21], F_Nchi_p[10,21], T_Nchi_p[10,21],I_Nchi_p[10,21],P_Nchi_p[10,21];
float A_COchi_f[10,21], Q_COchi_f[10,21], G_COchi_f[10,21], V_COchi_f[10,21];
float N_chi12_s[30,21], Ca_chi12_s[30,21], Cb_chi12_s[30,21], CO_chi12_s[30,21];
float N_ref[20,3], Ca_ref[20,3], Cb_ref[20,3], CO_ref[20,3];
float N_rc[20], Ca_rc[20], Cb_rc[20], CO_rc[20]; // exp. ring current refs.

//  declarations for computed geometrical parameters:

#define MAXRES 10000
string R_name[ MAXRES ];
string RES_exp [ MAXRES ];
int first_resid;
string cf[ MAXRES ];
int S_num[ MAXRES ];
int R_num[ MAXRES ];
int R_relnum[ MAXRES ];
float phi[ MAXRES ];
float psi[ MAXRES ];
float chi[ MAXRES ];
float chi2[ MAXRES ];
float chi3[ MAXRES ];
float DHB1[ MAXRES ];
float DHB01[ MAXRES ];
float DHBA1[ MAXRES ];
float DHBT1[ MAXRES ];
float IHB1[ MAXRES ];
float IHBA1[ MAXRES ];
float IHBT1[ MAXRES ];
float DHB2[ MAXRES ];
float DHB02[ MAXRES ];
float DHBA2[ MAXRES ];
float DHBT2[ MAXRES ];
float IHB2[ MAXRES ];
float IHBA2[ MAXRES ];
float IHBT2[ MAXRES ];
float DHB03[ MAXRES ];
float DHBA3[ MAXRES ];
float DHBT3[ MAXRES ];
float IHB3[ MAXRES ];
float IHBA3[ MAXRES ];
float IHBT3[ MAXRES ];
float DHB04[ MAXRES ];
float DHBA4[ MAXRES ];
float DHBT4[ MAXRES ];
float IHB4[ MAXRES ];
float IHBA4[ MAXRES ];
float IHBT4[ MAXRES ];
float DWHB1[ MAXRES ];
float DWHBA1[ MAXRES ];
float DWHBT1[ MAXRES ];
float IWHB1[ MAXRES ];
float IWHBA1[ MAXRES ];
float IWHBT1[ MAXRES ];
float DWHB2[ MAXRES ];
float DWHBA2[ MAXRES ];
float DWHBT2[ MAXRES ];
float IWHB2[ MAXRES ];
float IWHBA2[ MAXRES ];
float IWHBT2[ MAXRES ];
float DWHB3[ MAXRES ];
float DWHBA3[ MAXRES ];
float DWHBT3[ MAXRES ];
float IWHB3[ MAXRES ];
float IWHBA3[ MAXRES ];
float IWHBT3[ MAXRES ];
float DWHB4[ MAXRES ];
float DWHBA4[ MAXRES ];
float DWHBT4[ MAXRES ];
float IWHB4[ MAXRES ];
float IWHBA4[ MAXRES ];
float IWHBT4[ MAXRES ];
float DWHB5[ MAXRES ];
float DWHBA5[ MAXRES ];
float DWHBT5[ MAXRES ];
float IWHB5[ MAXRES ];
float IWHBA5[ MAXRES ];
float IWHBT5[ MAXRES ];
float DWHB6[ MAXRES ];
float DWHBA6[ MAXRES ];
float DWHBT6[ MAXRES ];
float IWHB6[ MAXRES ];
float IWHBA6[ MAXRES ];
float IWHBT6[ MAXRES ];
int  CI1[ MAXRES ];
int  CD1[ MAXRES ];
int  CI2[ MAXRES ];
int  CD2[ MAXRES ];
int  CI3[ MAXRES ];
int  CD3[ MAXRES ];
int  CWI1[ MAXRES ];
int  CWD1[ MAXRES ];
int  CWI2[ MAXRES ];
int  CWD2[ MAXRES ];
int  CWI3[ MAXRES ];
int  CWD3[ MAXRES ];
int  CWI4[ MAXRES ];
int  CWD4[ MAXRES ];
int  CWI5[ MAXRES ];
int  CWD5[ MAXRES ];
int  class_cf[ MAXRES ];
int  range[ MAXRES ];
int  class_res[ MAXRES ];
int  class_chi[ MAXRES ];
int  class_chi2[ MAXRES ];
int  class_chi12[ MAXRES ];
float N_close[MAXRES], O_close[MAXRES];
int   No_close[MAXRES];
int firstwat, lastwat, firstprot, lastprot;

// prediction and observation parameter declarations:

float incf[ MAXRES ];
float dN_bb[MAXRES,3], dCa_bb[MAXRES,3], dCb_bb[MAXRES,3], dCO_bb[MAXRES,3];
float dN_chi_p[MAXRES], dCa_chi_p[MAXRES], dCb_chi_p[MAXRES], dCO_chi_p[MAXRES], dCO_chi_f[MAXRES];
float dN_chi12_s[MAXRES], dCa_chi12_s[MAXRES], dCb_chi12_s[MAXRES],
	dCO_chi12_s[MAXRES];
float dHB[MAXRES], dHB_D[MAXRES], dHB_I[MAXRES], dHB_T[MAXRES], dHB_WD[MAXRES]; 
float dHB_D1[MAXRES], dHB_D2[MAXRES], dHB_D3[MAXRES], dHB_D4[MAXRES];
float dHB_I1[MAXRES], dHB_I2[MAXRES], dHB_I3[MAXRES], dHB_I4[MAXRES];
float dWHB_D1[MAXRES], dWHB_D2[MAXRES], dWHB_D3[MAXRES], dWHB_D4[MAXRES], dWHB_D5[MAXRES], dWHB_D6[MAXRES];
float dWHB_I1[MAXRES], dWHB_I2[MAXRES], dWHB_I3[MAXRES], dWHB_I4[MAXRES], dWHB_I5[MAXRES], dWHB_I6[MAXRES];
float CO_dHB[MAXRES], CO_dHB_D[MAXRES], CO_dHB_I[MAXRES];
float dN_Sum[MAXRES], dCa_Sum[MAXRES], dCb_Sum[MAXRES], dCO_Sum[MAXRES];
float N_REF[MAXRES], Ca_REF[MAXRES], Cb_REF[MAXRES], CO_REF[MAXRES];
float N_RC[MAXRES], Ca_RC[MAXRES], Cb_RC[MAXRES], CO_RC[MAXRES];

float N_pred[MAXRES], Ca_pred[MAXRES], Cb_pred[MAXRES], CO_pred[MAXRES];
float N_exp[MAXRES], Ca_exp[MAXRES], Cb_exp[MAXRES], CO_exp[MAXRES];

// refinement parameter declarations:

int mod_ID[MAXRES], N_ID, Ca_ID, Cb_ID, CO_ID;
float diff_N[MAXRES], diff_Ca[MAXRES], diff_Cb[MAXRES], diff_CO[MAXRES];
float new_ID[30], old_ID, mid_ID;

float  min_D[MAXRES], min_A[MAXRES], min_T[MAXRES];
float min_WD[MAXRES], min_WA[MAXRES], min_WT[MAXRES];


