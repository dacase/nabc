/* divcon CRT */

#include "eprintf.h"


int rdivcrt(char *filename, int *atomnum, ATOM atom[], CONTROLINFO cinfo, MOLINFO minfo)
{
    FILE *fpin;
    int index;
    int numatom;
    int overflow_flag = 0;
    char line[MAXCHAR];
    char tmpchar1[MAXCHAR];
    char tmpchar2[MAXCHAR];
    double tmpfloat1, tmpfloat2, tmpfloat3;

    fpin = efopen(filename, "r");
    initial(cinfo.maxatom, atom, minfo.resname);
    index = 0;
    numatom = 0;
    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL) {
            break;
        }
        index++;
        if (index <= 2)
            continue;
        sscanf(line, "%s%s%lf%lf%lf", tmpchar1, tmpchar2, &tmpfloat1, &tmpfloat2,
               &tmpfloat3);
        if (strcmp(tmpchar1, "END_COORD") == 0)
            break;
        if (overflow_flag == 0) {
            strcpy(atom[numatom].name, tmpchar2);
            atom[numatom].x = tmpfloat1;
            atom[numatom].y = tmpfloat2;
            atom[numatom].z = tmpfloat3;
        }
        numatom++;
        if (numatom >= cinfo.maxatom && overflow_flag == 0) {
            printf("Info: The number of atoms (%d) exceeded MAXATOM.\n", numatom);
            overflow_flag = 1;
        }
    }
    if (cinfo.intstatus == 2)
        printf("Info: Finished reading file (%s); lines read (%d), atoms read (%d).\n",
               filename, index, numatom);
    fclose(fpin);
    *atomnum = numatom;
    return overflow_flag;
}

void wdivcrt(char *filename, int atomnum, ATOM atom[], MOLINFO minfo)
{
    FILE *fpout;
    int i, nelectrons;
    char tmpkeyword[MAXCHAR];

    fpout = efopen(filename, "w");

    for (i = 0; i < strlen(minfo.ekeyword); i++)
        tmpkeyword[i] = toupper(minfo.ekeyword[i]);
    if (strstr(tmpkeyword, "CHAR") == NULL) {
        fprintf(fpout, "%s", minfo.ekeyword);
        fprintf(fpout, " CHARGE=%d", minfo.icharge);
    }
    if (strstr(tmpkeyword, "DOUBLET") == NULL && minfo.multiplicity == 2)
        fprintf(fpout, " DOUBLET");
    if (strstr(tmpkeyword, "TRIPLET") == NULL && minfo.multiplicity == 3)
        fprintf(fpout, " TRIPLET");
/*
#if 0
	fprintf(fpout, "\ncreated by wmopcrt()\n\n" );
#else
	fprintf(fpout, "\ncreated by wmopcrt()\n");
#endif
*/
    fprintf(fpout, "\ncreated by wmopcrt() for divcon\n");
    initialize_elements_in_atom_to_symbols_upto_atomnum(atomnum, atom);
    nelectrons = 0;
    for (i = 0; i < atomnum; i++) {
        fprintf(fpout, "%5d %5s %12.4lf  %12.4lf  %12.4lf\n", i + 1, atom[i].element,
                atom[i].x, atom[i].y, atom[i].z);
/*  check that the number of electrons is even:   */
        nelectrons += atom[i].atomicnum;
    }
    fprintf(fpout, "END_COORD\n\n");
    nelectrons -= minfo.icharge;
    fprintf(stdout, "Total number of electrons: %d; net charge: %d\n", nelectrons,
            minfo.icharge);
    if (nelectrons % 2 != 0) {
        fprintf(stdout, "Number of electrons is odd: %d\n", nelectrons);
        fprintf(stdout, "Please check the total charge and your -nc flag\n");
        exit(1);
    }
    fclose(fpout);
}
