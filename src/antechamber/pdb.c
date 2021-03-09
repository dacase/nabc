/* PDB */
/*
1 - 6 Record name "ATOM "
7 - 11 Integer serial Atom serial number.
13 - 16 Atom name Atom name.
17 Character altLoc Alternate location indicator.
18 - 20 Residue name resName Residue name.
22 Character chainID Chain identifier.
23 - 26 Integer resSeq Residue sequence number.
27 AChar iCode Code for insertion of residues.
31 - 38 Real(8.3) x Orthogonal coordinates for X in Angstroms
39 - 46 Real(8.3) y Orthogonal coordinates for Y in Angstroms
47 - 54 Real(8.3) z Orthogonal coordinates for Z in Angstroms
55 - 60 Real(6.2) occupancy Occupancy.
61 - 66 Real(6.2) tempFactor Temperature factor.
77 - 78 LString(2) element Element symbol, right-justified.
79 - 80 LString(2) charge Charge on the atom.
12345678901234567890123456789012345678901234567890123456789012345678901234567890
ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92           N
ATOM    146  CA  VAL A  25      31.132  16.439  58.160  1.00 11.85           C  
*/

#include <ctype.h>
#include "eprintf.h"


// Stores the trimmed input string into the given output buffer, which must be
// large enough to store the result.  If it is too small, the output is
// truncated.
size_t trimwhitespace(char *out, size_t len, const char *str)
{
    if (len == 0)
        return 0;

    const char *end;
    size_t out_size;

    // Trim leading space
    while (isspace(*str))
        str++;

    if (*str == 0)              // All spaces?
    {
        *out = 0;
        return 1;
    }
    // Trim trailing space
    end = str + strlen(str) - 1;
    while (end > str && isspace(*end))
        end--;
    end++;

    // Set output size to minimum of trimmed string length and buffer size minus 1
    out_size = (end - str) < len - 1 ? (end - str) : len - 1;

    // Copy trimmed string and add null terminator
    memcpy(out, str, out_size);
    out[out_size] = 0;

    return out_size;
}

int rpdb(char *filename, int *atomnum, ATOM * atom, CONTROLINFO cinfo, MOLINFO minfo,
         int pqr)
{
    int numatom;
    int terindex;
    int overflow_flag = 0;
    int warning_flag = 1;
    int resno = -1;             /* used to check for a multiple residue PDB file */
    int resno_counter = 0;      /* used to check for a multiple residue PDB file */
    const int MULTIPLE_RESNO_WARNING_TOLERANCE = 10;
    char tmpchar[MAXCHAR];
    char line[MAXCHAR];
    char elem[MAXCHAR];
    char atomname[MAXCHAR];
    char atomname2[MAXCHAR];
    char resname[MAXCHAR];
    double tmpfloat1, tmpfloat2;
    double x, y, z;
    int id;
    int tmpint, tmpint1, tmpint2, tmpint3, tmpint4, tmpint5, tmpint6, tmpint7;
    int ielem = 1;
    int i;
    int len;
    FILE *fpin;

    fpin = efopen(filename, "r");
    numatom = 0;
    initial(cinfo.maxatom, atom, minfo.resname);
    terindex = -1;
    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL) {
            /*  printf("\nFinished reading %s file.", filename); */
            break;
        }
        if (strncmp("TER", line, 3) == 0) {
            terindex = 1;
            continue;
        }
        if (strncmp("ATOM", line, 4) == 0 || strncmp("HETATM", line, 6) == 0) {
            id = -1;
            if (overflow_flag == 0) {
                line[26] = ' ';
                if (strchr("0123456789", line[5]) != NULL) {
                    sscanf(line + 5, "%d", &id);        /* 'ATOM\s\d+' */
                } else {
                    sscanf(line + 6, "%d", &id);        /* 'HETATM' or 'ATOM\s\s' */
                }
                if (pqr)
                    sscanf(&line[22], "%d%lf%lf%lf%lf%lf%s", &tmpint2, &x, &y, &z,
                           &tmpfloat1, &tmpfloat2, tmpchar);
                else {
                    if (cinfo.verify_pdb_atomname == 1 && ielem == 1) {
                        if (strlen(line) < 1 + 77)
                            ielem = 0;
                        else if (line[76] == ' ' && line[77] >= 'A' && line[77] <= 'Z') {
                            elem[0] = line[77];
                            elem[1] = '\0';
                        } else if (line[76] >= 'A' && line[76] <= 'Z' && line[77] == ' ') {
                            elem[0] = line[76];
                            elem[1] = '\0';
                        } else if (line[76] >= 'A' && line[76] <= 'Z'
                                   && ((line[77] >= 'a' && line[77] <= 'z')
                                       || (line[77] >= 'A' && line[77] <= 'Z'))) {
                            elem[0] = line[76];
                            elem[1] = line[77];
                            elem[2] = '\0';
                        } else
                            ielem = 0;
                    }
                    sscanf(&line[22], "%d%lf%lf%lf", &tmpint2, &x, &y, &z);

                }

/*          --- columns 13-16 have the Brookhaven-formatted name:    */
                atomname[0] = line[12];
                atomname[1] = line[13];
                atomname[2] = line[14];
                atomname[3] = line[15];
                atomname[4] = '\0';
/*          --- now unwrap this to a more understandable convention:    */
                sscanf(atomname, "%s", atomname2);
                if (atomname2[0] >= '0' && atomname2[0] <= '9') {
                    for (i = 1; i < strlen(atomname2); i++)
                        atom[numatom].name[i - 1] = atomname2[i];
                    atom[numatom].name[i - 1] = atomname2[0];
                    atom[numatom].name[i] = '\0';
                } else
                    strcpy(atom[numatom].name, atomname2);

                if (cinfo.rnindex == 0) {
                    resname[0] = line[17];
                    resname[1] = line[18];
                    resname[2] = line[19];
                    resname[3] = '\0';
                    sscanf(resname, "%s", atom[numatom].aa);
                }
                if (line[21] != ' ')
                    atom[numatom].chain[0] = line[21];
                atom[numatom].id = id;
                atom[numatom].ter = terindex;
                atom[numatom].resno = tmpint2;
                if (atom[numatom].resno != resno) {
                    resno = atom[numatom].resno;
                    ++resno_counter;
                    if (warning_flag == 1) {
                        if (resno_counter > MULTIPLE_RESNO_WARNING_TOLERANCE) {
                            fprintf(stdout,
                                    "Warning: Detected more than %d "
                                    "Residue sequence numbers;\n"
                                    "         this may be a large multiple residue PDB file;\n"
                                    "         large multiple residue PDB files are not supported.\n"
                                    "         This warning usually indicates a conceptual"
                                    " misunderstanding.\n"
                                    "         We recommend reviewing the Information flow in"
                                    " Amber documentation\n"
                                    "         and the antechamber tutorials.\n"
                                    "         Continuing, but problems may be encountered.\n",
                                    MULTIPLE_RESNO_WARNING_TOLERANCE);
                            warning_flag = 0;
                        }
                    }
                }
                atom[numatom].x = x;
                atom[numatom].y = y;
                atom[numatom].z = z;
                if (line[16] != ' ' && numatom >= 1) {
                    if (atom[numatom].id == atom[numatom - 1].id
                        && strcmp(atom[numatom].name, atom[numatom - 1].name) == 0
                        && strcmp(atom[numatom].aa, atom[numatom - 1].aa) == 0)
                        continue;

                }
                if (terindex == 1) {
                    atom[numatom].ter = terindex;
                    terindex = -1;
                }
                if (pqr) {
                    atom[numatom].charge = tmpfloat1;
                    atom[numatom].radius = tmpfloat2;
                    strcpy(atom[numatom].ambername, tmpchar);
                }

                if (strcmp(atom[numatom].name, "dumm") == 0)
                    continue;
                if (strcmp(atom[numatom].name, "Du") == 0)
                    continue;
                if (strcmp(atom[numatom].name, "DUMM") == 0)
                    continue;
                if (cinfo.verify_pdb_atomname == 1 && ielem == 1) {
                    len = strlen(elem);
                    if (strncmp(atom[numatom].name, elem, len) != 0) {
                        for (i = 0; i < len; i++)
                            atom[numatom].name[i] = elem[i];
                    }
                }
            }
            numatom++;
            if (numatom >= cinfo.maxatom && overflow_flag == 0) {
                printf("Info: The number of atoms (%d) exceeded MAXATOM.\n", numatom);
                overflow_flag = 1;
            }
        }

    }
    *atomnum = numatom;
    rewind(fpin);
    for (; overflow_flag == 0;) {
        if (fgets(line, MAXCHAR, fpin) == NULL)
            break;
        if (strncmp("CONECT", line, 6) == 0) {
            tmpint1 = -1;
            tmpint2 = -1;
            tmpint3 = -1;
            tmpint4 = -1;
            tmpint5 = -1;
            tmpint6 = -1;
            tmpint7 = -1;
            sscanf(&line[6], "%d%d%d%d%d%d%d", &tmpint1, &tmpint2, &tmpint3, &tmpint4,
                   &tmpint5, &tmpint6, &tmpint7);
            if (tmpint1 == -1)
                continue;
            else {
                tmpint = -1;
                for (i = 0; i < numatom; i++)
                    if (atom[i].id == tmpint1) {
                        tmpint = i;
                        break;
                    }
                if (tmpint == -1) {
                    eprintf("Invalid atom id (%d) in line (%s).", tmpint1, line);
                } else
                    tmpint1 = tmpint;

            }

            if (tmpint2 == -1 || tmpint2 == 0)
                continue;
            else {
                tmpint = -1;
                for (i = 0; i < numatom; i++)
                    if (atom[i].id == tmpint2) {
                        tmpint = i;
                        break;
                    }
                if (tmpint == -1) {
                    eprintf("Invalid atom id (%d) in line (%s).", tmpint2, line);
                } else
                    atom[tmpint1].con[0] = tmpint;

            }

            if (tmpint3 == -1 || tmpint3 == 0) {
                atom[tmpint1].connum = 1;
                continue;
            } else {
                tmpint = -1;
                for (i = 0; i < numatom; i++)
                    if (atom[i].id == tmpint3) {
                        tmpint = i;
                        break;
                    }
                if (tmpint == -1) {
                    eprintf("Invalid atom id (%d) in line (%s).", tmpint3, line);
                } else
                    atom[tmpint1].con[1] = tmpint;

            }

            if (tmpint4 == -1 || tmpint4 == 0) {
                atom[tmpint1].connum = 2;
                continue;
            } else {
                tmpint = -1;
                for (i = 0; i < numatom; i++)
                    if (atom[i].id == tmpint4) {
                        tmpint = i;
                        break;
                    }
                if (tmpint == -1) {
                    eprintf("Invalid atom id (%d) in line (%s).", tmpint4, line);
                } else
                    atom[tmpint1].con[2] = tmpint;

            }

            if (tmpint5 == -1 || tmpint5 == 0) {
                atom[tmpint1].connum = 3;
                continue;
            } else {
                tmpint = -1;
                for (i = 0; i < numatom; i++)
                    if (atom[i].id == tmpint5) {
                        tmpint = i;
                        break;
                    }
                if (tmpint == -1) {
                    eprintf("Invalid atom id (%d) in line (%s).", tmpint5, line);
                } else
                    atom[tmpint1].con[3] = tmpint;

            }

            if (tmpint6 == -1 || tmpint6 == 0) {
                atom[tmpint1].connum = 4;
                continue;
            } else {
                tmpint = -1;
                for (i = 0; i < numatom; i++)
                    if (atom[i].id == tmpint6) {
                        tmpint = i;
                        break;
                    }
                if (tmpint == -1) {
                    eprintf("Invalid atom id (%d) in line (%s).", tmpint6, line);
                } else
                    atom[tmpint1].con[4] = tmpint;

            }

            if (tmpint7 == -1 || tmpint7 == 0) {
                atom[tmpint1].connum = 5;
                continue;
            } else {
                tmpint = -1;
                for (i = 0; i < numatom; i++)
                    if (atom[i].id == tmpint7) {
                        tmpint = i;
                        break;
                    }
                if (tmpint == -1) {
                    eprintf("Invalid atom id (%d) in line (%s).", tmpint7, line);
                } else {
                    atom[tmpint1].con[5] = tmpint;
                    atom[tmpint1].connum = 6;
                }
            }
        }
    }
    fclose(fpin);
    return overflow_flag;

}



void wpdb(char *filename, int atomnum, ATOM atom[])
{
    FILE *fpout;
    int i, j;
    char atmp[5], atmp2[5];
    fpout = efopen(filename, "w");
    for (i = 0; i < atomnum; i++) {
        if (atom[i].ter == 1)
            fprintf(fpout, "TER\n");
        /* use PDB v3 rules to pack atom name into 4 spaces:  */
        trimwhitespace(atmp, 5, atom[i].name);
        if (strlen(atmp) < 4 && strlen(atom[i].element) < 2) {
            atmp2[0] = ' ';
            atmp2[1] = atmp[0];
            atmp2[2] = atmp[1];
            atmp2[3] = atmp[2];
            atmp2[4] = '\0';
        } else {
            strncpy(atmp2, atmp, 4);
        }
        fprintf(fpout, "ATOM%7d %-4s %3s  %4d    %8.3f%8.3f%8.3f%6.2lf%6.2lf%12s\n",
                i + 1, atmp2, atom[i].aa, atom[i].resno, atom[i].x, atom[i].y, atom[i].z,
                1.0, 0.0, atom[i].element);
    }

    for (i = 0; i < atomnum; i++) {
        if (atom[i].connum > 0)
            fprintf(fpout, "CONECT%5d", i + 1);
        else
            continue;
        for (j = 0; j < 6; j++) {
            if (atom[i].con[j] < 0)
                continue;
            else
                fprintf(fpout, "%5d", atom[i].con[j] + 1);
        }
        fprintf(fpout, "\n");
    }

    fclose(fpout);
}

void wmpdb(char *filename, int atomnum, ATOM atom[])
{
    FILE *fpout;
    int i, j;
    fpout = efopen(filename, "w");
    for (i = 0; i < atomnum; i++) {
        if (atom[i].ter == 1)
            fprintf(fpout, "TER\n");
        fprintf(fpout, "ATOM%7d  %-4s%-4s%5d%12.3f%8.3f%8.3f%10.6lf%8.2lf%8s\n", i + 1,
                atom[i].name, atom[i].aa, atom[i].resno, atom[i].x, atom[i].y, atom[i].z,
                atom[i].charge, atom[i].radius, atom[i].ambername);
    }

    for (i = 0; i < atomnum; i++) {
        if (atom[i].connum > 0)
            fprintf(fpout, "CONECT%5d", i + 1);
        else
            continue;
        for (j = 0; j < 6; j++) {
            if (atom[i].con[j] < 0)
                continue;
            else
                fprintf(fpout, "%5d", atom[i].con[j] + 1);
        }
        fprintf(fpout, "\n");
    }
    fclose(fpout);
}
