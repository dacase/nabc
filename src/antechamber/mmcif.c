/* Component mmcif files */
/* This file doubles as a header file and as an included (!) source file. */
/* Yes, we are on the road to multiple header files and separate compilation. */

#ifndef MMCIF_H
#define MMCIF_H

#define MAX_CIF_BLOCK_CODE_LENGTH 32

extern char blockId[MAX_CIF_BLOCK_CODE_LENGTH + 1];

#endif                          // MMCIF_H


#ifndef MMCIF_C__YES_THAT_IS_C
#define MMCIF_C__YES_THAT_IS_C

char blockId[MAX_CIF_BLOCK_CODE_LENGTH + 1] = "Undefined_block_code__use-bk";

/*  Read a "component cif" file from the RCSB components library  */

#include <assert.h>
#include "../cifparse/cifparse.h"
#include "eprintf.h"

int rmmcif(char *filename, char *blockId, int *atomnum, ATOM atom[], int *bondnum,
           BOND bond[], CONTROLINFO * cinfo, MOLINFO * minfo, int flag)
{
/*if flag =1, read in atom type, if flag ==0, do not read in atom type */

    int i, numatom, numbond, iblock, iCat;
    int col2, col3, col6, col7, col8;
    FILE *fpin;

    fpin = efopen(filename, "r");
    initial((*cinfo).maxatom, atom, (*minfo).resname);
    numatom = 0;
    numbond = 0;

    ndb_cif_init();
    ndb_cif_read_file(fpin);

    /*  Check for target data block.  */
    iblock = -1;
    for (i = 0; i < cifFiles.numDatablock; i++) {
        if (!strcmp(blockId, cifFiles.datablocks[i].datablockName)) {
            iblock = i;
            break;
        }
    }
    if (iblock == -1) {
        eprintf("Target data block (%s) is not in file.", blockId);
    }

    /* get chem_comp overall info: */
    iCat = get_category_index(iblock, "chem_comp");
    if (iCat == -1) {
        eprintf("Target category (%s) is not in data block (%d).", "atom_site",
                iblock);
    }

    /* get the formal charge into the antechamber minfo struct: */
    col2 = get_column_index(iblock, iCat, "pdbx_formal_charge");
    assert(col2 >= 0);
    (*minfo).usercharge =
        atoi(cifFiles.datablocks[iblock].categories[iCat].rows[0].columns[col2]);

    /* Exit (for now!?) if this is a "linking" -type residue */
    col2 = get_column_index(iblock, iCat, "type");
    assert(col2 >= 0);
    if (strstr
        (cifFiles.datablocks[iblock].categories[iCat].rows[0].columns[col2], "LINKING")
        || strstr(cifFiles.datablocks[iblock].categories[iCat].rows[0].columns[col2],
                  "linking")) {
        eprintf("Residue (%s) has a type of LINKING.", blockId);
    }

    /* get chem_comp_atom atomic info: */
    iCat = get_category_index(iblock, "chem_comp_atom");
    if (iCat == -1) {
        eprintf("Target category (%s) is not in data block (%d).", "atom_site",
                iblock);
    }
    /* fprintf(stderr, "Read %d rows in category %d\n",
       cifFiles.datablocks[iblock].categories[iCat].numRow, iCat); */

    col2 = get_column_index(iblock, iCat, "atom_id");   /* atom name */
    assert(col2 >= 0);
    col3 = get_column_index(iblock, iCat, "comp_id");   /* res. name */
    assert(col3 >= 0);

#define IDEAL                   /* undefine to get the regular "model" coordinates */
#ifdef IDEAL
    col6 = get_column_index(iblock, iCat, "pdbx_model_Cartn_x_ideal");  /* x-coord */
    assert(col6 >= 0);
    col7 = get_column_index(iblock, iCat, "pdbx_model_Cartn_y_ideal");  /* y-coord */
    assert(col7 >= 0);
    col8 = get_column_index(iblock, iCat, "pdbx_model_Cartn_z_ideal");  /* z-coord */
    assert(col8 >= 0);
#else
    col6 = get_column_index(iblock, iCat, "model_Cartn_x");     /* x-coord */
    assert(col6 >= 0);
    col7 = get_column_index(iblock, iCat, "model_Cartn_y");     /* y-coord */
    assert(col7 >= 0);
    col8 = get_column_index(iblock, iCat, "model_Cartn_z");     /* z-coord */
    assert(col8 >= 0);
#endif

    /*  --- need some better error processing here!  */

    for (i = 0; i < cifFiles.datablocks[iblock].categories[iCat].numRow; i++) {

        /* check for missing coordinates; then skip this atom  */
        if (strlen(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col6]) ==
            0
            || strlen(cifFiles.datablocks[iblock].categories[iCat].rows[i].
                      columns[col7]) == 0
            || strlen(cifFiles.datablocks[iblock].categories[iCat].rows[i].
                      columns[col8]) == 0) {
            /*  fprintf( stderr, "skipping %d\n", i );  */
            continue;
        }

        strcpy(atom[numatom].name,
               cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col2]);
        strcpy(atom[numatom].aa,
               cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col3]);
        atom[numatom].resno = 1;
        atom[numatom].x =
            atof(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col6]);
        atom[numatom].y =
            atof(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col7]);
        atom[numatom].z =
            atof(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col8]);

        numatom++;

    }

    fclose(fpin);
    *atomnum = numatom;
    *bondnum = numbond;
/*
	if(overflow_flag == 0) {
		tmpf1 = 0;
		tmpf2 = 0;
		for(i=0;i<numatom;i++) {
			if(atom[i].charge > 0)
				tmpf1+=atom[i].charge;
			else
				tmpf2+=atom[i].charge;
		}
		printf("\nchargep is %9.4lf", tmpf1);
		printf("\nchargen is %9.4lf", tmpf2);
	}
*/
    return 0;
}

#if 0

void wmmcif(char *filename, int atomnum, ATOM atom[], int bondnum, BOND bond[],
            AROM arom[], CONTROLINFO cinfo, MOLINFO minfo)
{
    int flag;
    int bondi, bondj;
    FILE *fpout;
    if ((fpout = fopen(filename, "w")) == NULL) {
        fprintf(stdout, "Warning: Cannot open file (%s) to write.\n", filename);
        return;
    }
    fclose(fpout);
}

#endif

#endif                          // MMCIF_C__YES_THAT_IS_C
