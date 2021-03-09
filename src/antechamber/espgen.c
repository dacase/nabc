/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    espgen                                                     *
*  Version: version 1.1                                                *
*  Author:  Junmei Wang                                                *
*  GAMESS support added by Robin Betz, 2019                            *
*                                                                      *
*  Department of Pharmaceutical Chemistry                              *
*  School of Pharmacy                                                  *
*  University of California                                            *
*  San Francisco   CA 94143                                            *
*  Octomber, 2001                                                      *
************************************************************************
*/
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "define.h"
# include "atom.h"
# include "eprintf.h"
# define MAXESP 99999
# define MAX_ATOM_CENTER 1000

int i, j, k;
int fail = 0;
int method = 0;                 /* output dipole and quadrupole? */
int maxesp = 0;
int format = 1;                 /* 1: log files of G98, G03, 2: G09 ESP file, specia output triggered by iop(6/50=1) */
char line[MAXCHAR];
char ifilename[MAXCHAR];
char ofilename[MAXCHAR];

double cord[MAX_ATOM_CENTER][3];
ESP *esp;
DM dipole;
QM quadrupole;
FILE *fpin, *fpout;

int Found_Stationary = 0;
int esp_index = 0;
int espvalue_index = 0;
int g09_index = 0;
int opt_index = 0;
int npoint = 0;
int natom = 0;
int nset = 0;

char tmpchar1[MAXCHAR], tmpchar2[MAXCHAR], tmpchar3[MAXCHAR];
char tmpchar4[MAXCHAR], tmpchar5[MAXCHAR], tmpchar6[MAXCHAR];
char tmpchar7[MAXCHAR], tmpchar8[MAXCHAR], tmpchar9[MAXCHAR];
int tmpint1, tmpint2, tmpint3;
double tmpfloat1, tmpfloat2, tmpfloat3, tmpfloat4;

void rgamessdat(void)
{
    int atoms = 0, points = 0;
    int i;

    rewind(fpin);

    // Dipole and quadrupole currently unimplemented in this parser
    if (method == 1)
        eprintf("Dipole and quadrupole support unimplemented in GAMESS parser.");

    // Find line with number of atomic centers
    for (;;) {
        if (!fgets(line, MAXCHAR, fpin))
            eprintf("Premature end of file before electrostatic potential");

        // ELECTROSTATIC POTENTIAL COMPUTED FOR <N> ATOMS, TOTAL CHARGE= <N>
        if (strncmp(line, " ELECTROSTATIC POTENTIAL COMPUTED FOR", 37) == 0) {
            sscanf(line, "%*s%*s%*s%*s%d", &atoms);
            if (atoms > MAX_ATOM_CENTER) {
                eprintf("The number of atomic centers (%d) exceeded "
                        "MAX_ATOM_CENTER.\n Increase MAX_ATOM_CENTER in "
                        "espgen.c and rebuild.", atoms);
            }
            break;
        }
    }

    // Read in each atomic center
    for (i = 0; i < atoms; i++) {
        if (!fgets(line, MAXCHAR, fpin))
            eprintf("Premature end of file at atom %d", i);

        // Line has index, atomicnum, x, y, z coordinates
        sscanf(line, "%*d%*f%lf%lf%lf", &cord[i][0], &cord[i][1], &cord[i][2]);
    }

    // Throw out "ELECTROSTATIC POTENTIAL, IPT,X,Y,Z,ELPOTT" line
    if (!fgets(line, MAXCHAR, fpin))
        eprintf("Premature end of file before grid points definition");

    // Read in number of grid points
    if (!fgets(line, MAXCHAR, fpin))
        eprintf("Premature end of file before grid points number");

    if (strncmp(line, " TOTAL NUMBER OF GRID POINTS NPT=", 33))
        eprintf("Grid points line not found as expected");

    sscanf(line, "%*s%*s%*s%*s%*s%*s%d", &points);

    if (points >= maxesp) {
        printf("\nInfo: The number of ESP exceeded MAXESP of %d; "
               "automatically increasing to %d.\n", maxesp, maxesp + MAXESP);
        maxesp += MAXESP;
        esp = (ESP *) erealloc(esp, maxesp * sizeof(ESP));
    }

    // Read in each grid point
    for (i = 0; i < points; i++) {
        if (!fgets(line, MAXCHAR, fpin))
            eprintf("Premature end of file reading ESP point %d", i);

        sscanf(line, "%*d%lf%lf%lf%lf", &esp[i].x, &esp[i].y, &esp[i].z,
               &esp[i].esp);
    }

    // Now write out to the output file
    // First line is number of atoms, number of points, and 0
    fprintf(fpout, "%5d%5d%5d\n", atoms, points, 0);

    // Then all of the atomic coordinates
    // These are already in units of Bohr
    for (i = 0; i < atoms; i++) {
        fprintf(fpout, "%32.7E%16.7E%16.7E\n", cord[i][0], cord[i][1],
                cord[i][2]);
    }

    // Finally, all of the ESP points, which are already in units of Bohr
    for (i = 0; i < points; i++) {
        fprintf(fpout, "%16.7E%16.7E%16.7E%16.7E\n", esp[i].esp,
                esp[i].x, esp[i].y, esp[i].z);
    }

    printf("Success!");
}

void rglog(void)
{
    int i, j;

    tmpint1 = 0;
    tmpint2 = 0;
    tmpint3 = 0;
    rewind(fpin);
    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL)
            break;
        sscanf(&line[0], "%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4, tmpchar5);
        if (strcmp("--", tmpchar1) == 0 && strcmp("Stationary", tmpchar2) == 0
            && strcmp("point", tmpchar3) == 0 && strcmp("found.", tmpchar4) == 0)
            Found_Stationary = 1;
        if (opt_index == 1 && nset == 2 && Found_Stationary == 0)
            continue;
        if (strcmp("Atomic", tmpchar1) == 0 && strcmp("Center", tmpchar2) == 0) {
            sscanf(&line[32], "%lf%lf%lf", &cord[tmpint1][0], &cord[tmpint1][1],
                   &cord[tmpint1][2]);
            tmpint1++;
            if (tmpint1 > MAX_ATOM_CENTER) {
                eprintf("The number of atomic centers (%d) exceeded MAX_ATOM_CENTER.\n"
                        "Increase MAX_ATOM_CENTER in espgen.c and rebuild.", tmpint1);
            }
        }
        if (strcmp("ESP", tmpchar1) == 0 && strcmp("Fit", tmpchar2) == 0
            && strcmp("Center", tmpchar3) == 0) {
            sscanf(&line[32], "%lf%lf%lf", &esp[tmpint2].x, &esp[tmpint2].y,
                   &esp[tmpint2].z);
            esp[tmpint2].x /= Bohr;
            esp[tmpint2].y /= Bohr;
            esp[tmpint2].z /= Bohr;
            tmpint2++;
            if (tmpint2 >= maxesp) {
                maxesp += MAXESP;
                printf("\nInfo: The number of ESP exceeded MAXESP; automatically increasing to %d.\n",
                       maxesp);
                esp = (ESP *) erealloc(esp, maxesp * sizeof(ESP));
            }
        }
        if (strncmp("Fit", &line[6], 3) == 0) {
            sscanf(&line[12], "%lf", &esp[tmpint3++].esp);
        }
    }
    if (tmpint3 < 10)           // no ESP at all
        fail = 1;
    if (fail == 0) {
        fprintf(fpout, "%5d%5d%5d\n", tmpint1, tmpint2, 0);
        for (i = 0; i < tmpint1; i++)
            fprintf(fpout, "%32.7E%16.7E%16.7E\n", cord[i][0] / Bohr, cord[i][1] / Bohr,
                    cord[i][2] / Bohr);

        for (j = 0; j < tmpint2; j++)
            fprintf(fpout, "%16.7E%16.7E%16.7E%16.7E\n", esp[j].esp, esp[j].x, esp[j].y,
                    esp[j].z);

    }
    if (fail == 0 && method == 1) {
        rewind(fpin);
        Found_Stationary = 0;
        for (;;) {
            if (fgets(line, MAXCHAR, fpin) == NULL)
                break;
            sscanf(&line[0], "%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4,
                   tmpchar5);
            if (strcmp("--", tmpchar1) == 0 && strcmp("Stationary", tmpchar2) == 0
                && strcmp("point", tmpchar3) == 0 && strcmp("found.", tmpchar4) == 0) {
                Found_Stationary = 1;
                continue;
            }
            if (opt_index == 0 || Found_Stationary == 1) {
                if (strcmp("Dipole", tmpchar1) == 0 && strcmp("moment", tmpchar2) == 0) {
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    sscanf(line, "%s%lf%s%lf%s%lf", tmpchar1, &dipole.x, tmpchar2,
                           &dipole.y, tmpchar3, &dipole.z);

                    continue;
                }
                if (strcmp("Traceless", tmpchar1) == 0
                    && strcmp("Quadrupole", tmpchar2) == 0
                    && strcmp("moment", tmpchar3) == 0) {
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    sscanf(line, "%s%lf%s%lf%s%lf", tmpchar1, &quadrupole.xx, tmpchar2,
                           &quadrupole.yy, tmpchar3, &quadrupole.zz);
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    sscanf(line, "%s%lf%s%lf%s%lf", tmpchar1, &quadrupole.xy, tmpchar2,
                           &quadrupole.xz, tmpchar3, &quadrupole.yz);
                    break;
                }
            }
        }
        fprintf(fpout, "%16.7E%16.7E%16.7E\n", dipole.x, dipole.y, dipole.z);
        fprintf(fpout, "%16.7E%16.7E%16.7E\n", quadrupole.xx, quadrupole.yy,
                quadrupole.zz);
        fprintf(fpout, "%16.7E%16.7E%16.7E\n", quadrupole.xy, quadrupole.xz,
                quadrupole.yz);
    }

}

double translate(char *str)
{
    int i;
    int pos = 0;
    int flag = 0;
    char tmpc1[MAXCHAR];
    char tmpc2[MAXCHAR];
    double f, e;

    for (i = 0; i < strlen(str); i++) {
        if (str[i] == 'D' || str[i] == 'd' || str[i] == 'E' || str[i] == 'e') {
            pos = i;
            flag = 1;
            continue;
        }
        if (flag == 0)
            tmpc1[i] = str[i];
        if (flag == 1)
            tmpc2[i - pos - 1] = str[i];

    }

    f = atof(tmpc1);
    e = atof(tmpc2);
    return f * pow(10.0, e);
}

void rgesp()
{
    int rflag = 0;
    int count = 0;
    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL)
            break;
        sscanf(line, "%s%s%s%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4, tmpchar5,
               tmpchar6, tmpchar7, tmpchar8);
        if (rflag == 3 && strcmp(tmpchar1, "ESP") == 0 && strcmp(tmpchar2, "VALUES") == 0) {
            rflag = 4;
            sscanf(line, "%s%s%s%s%s%s%s%s%d", tmpchar1, tmpchar2, tmpchar3, tmpchar4,
                   tmpchar5, tmpchar6, tmpchar7, tmpchar8, &npoint);
            if (npoint >= maxesp) {
                maxesp = npoint;
                printf("\nInfo: The number of ESP exceeded MAXESP; automatically increasing to %d.",
                       maxesp);
                esp = (ESP *) erealloc(esp, maxesp * sizeof(ESP));
            }
            continue;
        }
        if (rflag == 2 && strcmp(tmpchar1, "TRACELESS") == 0
            && strcmp(tmpchar2, "QUADRUPOLE") == 0) {
            rflag = 3;
            continue;
        }
        if (rflag == 1 && strcmp(tmpchar1, "DIPOLE") == 0) {
            rflag = 2;
            continue;
        }
        if (rflag == 0 && strncmp(line, " ATOMIC COORDINATES", 19) == 0) {
            rflag = 1;
            continue;
        }
        if (rflag == 1) {
            cord[natom][0] = translate(tmpchar2);
            cord[natom][1] = translate(tmpchar3);
            cord[natom][2] = translate(tmpchar4);
            natom++;
            if (natom >= MAX_ATOM_CENTER) {
                eprintf("The number of atomic centers (%d) exceeded MAX_ATOM_CENTER.\n"
                        "Increase MAX_ATOM_CENTER in espgen.c and rebuild.", natom);
            }
        }
        if (rflag == 2) {
            dipole.x = translate(tmpchar2);
            dipole.y = translate(tmpchar4);
            dipole.z = translate(tmpchar6);
            dipole.total = translate(tmpchar8);
        }
        if (rflag == 3) {
            quadrupole.xx = translate(tmpchar2);
            quadrupole.yy = translate(tmpchar4);
            quadrupole.zz = translate(tmpchar6);

            fgets(line, MAXCHAR, fpin);
            sscanf(line, "%s%s%s%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4,
                   tmpchar5, tmpchar6, tmpchar7, tmpchar8);
            quadrupole.xy = translate(tmpchar2);
            quadrupole.xz = translate(tmpchar4);
            quadrupole.yz = translate(tmpchar6);
        }
        if (rflag == 4) {
            if (count < npoint) {
                esp[count].esp = translate(tmpchar1);
                esp[count].x = translate(tmpchar2);
                esp[count].y = translate(tmpchar3);
                esp[count].z = translate(tmpchar4);
                count++;
            } else
                rflag = -1;
        }

    }
    fprintf(fpout, "%5d%5d%5d\n", natom, npoint, 0);
    for (i = 0; i < natom; i++)
        fprintf(fpout, "%32.7E%16.7E%16.7E\n", cord[i][0], cord[i][1], cord[i][2]);

    for (j = 0; j < npoint; j++)
        fprintf(fpout, "%16.7E%16.7E%16.7E%16.7E\n", esp[j].esp, esp[j].x, esp[j].y,
                esp[j].z);
    if (method == 1) {
        fprintf(fpout, "%16.7E%16.7E%16.7E\n", dipole.x, dipole.y, dipole.z);
        fprintf(fpout, "%16.7E%16.7E%16.7E\n", quadrupole.xx, quadrupole.yy,
                quadrupole.zz);
        fprintf(fpout, "%16.7E%16.7E%16.7E\n", quadrupole.xy, quadrupole.xz,
                quadrupole.yz);
    }
}

int main(int argc, char *argv[])
{
    if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
        if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-H") == 0)) {
            printf("[31mUsage: espgen -i  [0m input file name \n"
                   "[31m              -o  [0m output file name \n"
                   /*         "[31m              -dq [0m yes or no (optional, printing out dipole and\n"
                      "                   quadrupole moments or not, default is no\n" */
                );
            exit(0);
        }
        if (argc != 7 && argc != 5) {
            printf("[31mUsage: espgen -i  [0m input file name \n"
                   "[31m              -o  [0m output file name \n"
                   /*         "[31m              -dq [0m yes or no (optional, printing out dipole and\n"
                      "                   quadrupole moments or not, default is no\n" */
                );
            exit(1);
        }
    } else {
        if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-H") == 0)) {
            printf("Usage: espgen -i input file name\n");
            printf("              -o output file name \n");
/*
    printf("Usage: espgen -i  input_file_name  -o output file name \n");
    printf("              -dq yes or no (optional, printing out dipole and\n");
    printf("                  quadrupole moments or not, default is no\n");
*/
            exit(1);
        }
        if (argc != 7 && argc != 5) {
/*
    printf("Usage: espgen -i  input_file_name  -o output_file_name \n");
    printf("              -dq yes or no (optional, printing out dipole and\n");
    printf("                  quadrupole moments or not, default is no\n");
*/
            printf("Usage: espgen -i input file name\n");
            printf("              -o output file name \n");
            exit(1);
        }
    }

/* allocate memory for *esp */
    maxesp = MAXESP;
    esp = (ESP *) ecalloc(maxesp, sizeof(ESP));

    method = 0;
    for (i = 1; i < argc; i += 2) {
        if (strcmp(argv[i], "-i") == 0)
            strcpy(ifilename, argv[i + 1]);
        if (strcmp(argv[i], "-o") == 0)
            strcpy(ofilename, argv[i + 1]);
        if (strcmp(argv[i], "-dq") == 0) {
            if (strcmp("YES", argv[i + 1]) == 0 || strcmp("yes", argv[i + 1]) == 0
                || strcmp("1", argv[i + 1]) == 0)
                method = 1;
            if (strcmp("Y", argv[i + 1]) == 0 || strcmp("y", argv[i + 1]) == 0)
                method = 1;
            if (strcmp("NO", argv[i + 1]) == 0 || strcmp("no", argv[i + 1]) == 0
                || strcmp("0", argv[i + 1]) == 0)
                method = 0;
            if (strcmp("N", argv[i + 1]) == 0 || strcmp("n", argv[i + 1]) == 0)
                method = 0;
        }
    }
    fpin = efopen(ifilename, "r");
    fpout = efopen(ofilename, "w");
    nset = 0;
    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL)
            break;
        strcpy(tmpchar1, "");
        strcpy(tmpchar2, "");
        strcpy(tmpchar3, "");
        strcpy(tmpchar4, "");
        strcpy(tmpchar5, "");
        sscanf(&line[0], "%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4, tmpchar5);

        // First check if it is a GAMESS output file
        if (strncmp(line, " $DATA", 6) == 0) {
            format = 3;
            break;
        }

/*	adding code to detect Gaussion version	*/
/*      adding more code to deal with different Gaussiaon 09 subversions*/
        if (strncmp(line, " ESP FILE -", 11) == 0) {
            format = 2;
            break;
        }
        if (strncmp(line, " Gaussian 09", 12) == 0)
            g09_index = 1;
        if (strcmp("--", tmpchar1) == 0 && strcmp("Stationary", tmpchar2) == 0
            && strcmp("point", tmpchar3) == 0 && strcmp("found.", tmpchar4) == 0)
            opt_index = 1;
        if (strcmp("ESP", tmpchar1) == 0 && strcmp("Fit", tmpchar2) == 0
            && strcmp("Center", tmpchar3) == 0)
            esp_index = 1;
        if (strcmp("Electrostatic", tmpchar1) == 0 && strcmp("Properties", tmpchar2) == 0
            && strcmp("(Atomic", tmpchar3) == 0 && strcmp("Units)", tmpchar4) == 0)
            nset++;
        if (esp_index == 1 && espvalue_index == 0)
            if (strcmp("Fit", tmpchar2) == 0 && tmpchar1[0] >= '0' && tmpchar1[0] <= '9'
                && ((tmpchar3[0] >= '0' && tmpchar3[0] <= '9') || tmpchar3[0] == '-'))
                espvalue_index = 1;
    }

    if (format == 1) {
        if (esp_index == 0) {
            if (g09_index == 0)
                fprintf(stderr,
                        "Error: No ESP fitting centers and fitting values exist, adding 'iop(6/33=2) iop(6/42=6)' to the key word list\n");
            else
                fprintf(stderr,
                        "Error: No ESP fitting centers and fitting values exist, adding 'iop(6/33=2) iop(6/42=6) iop(6/50=1)' to the keyword list\n");
            exit(1);
        }
        if (esp_index == 1 && espvalue_index == 0) {
            fprintf(stderr,
                    "Error: the ESP fitting centers exist, but the fitting values are missing\n");
            if (g09_index == 1)
                fprintf(stderr,
                        "It is recommened to generate esp file for resp fitting from the gesp file generated by adding keyword 'iop(6/50=1) in G09 input'\n");
            exit(1);
        }
    }
    if (format == 1)
        rglog();
    if (format == 2)
        rgesp();
    if (format == 3)
        rgamessdat();

    fclose(fpin);
    fclose(fpout);
    if (fail == 1) {
/*		strcat(tmpchar, ofilename);
		esystem(tmpchar);
*/
        eprintf("Cannot run espgen properly.");

    }
/*
	free(esp);
*/
    printf("\n");
    return (0);
}
