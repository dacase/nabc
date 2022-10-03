#include <stdio.h>
#include "cifparse.h"

void main( int argc, char **argv )
{

ndb_cif_init();
ndb_cif_read_file( stdin );

//  make any needed edits here

ndb_cif_write_file( stdout );
ndb_cif_close();

}
