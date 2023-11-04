
		ringsize[ ring ] = split( ringinfo[ringname[ring]], ringatoms, " " ) - 3;
		assert( ringsize[ring] > 0 );
		assert( ringsize[ring] < 26 );

		ringnum[ ring ] = r.tresnum;
		ringstr[ ring ] = atof( ringatoms[ ringsize[ring] + 1 ]);
		ringcen[ ring ].x = ringcen[ ring ].y = ringcen[ ring ].z = 0.0;
		chi_anis[ ring ] = atof( ringatoms[ ringsize[ring] + 2 ]);
		chi_iso[ ring ] = atof( ringatoms[ ringsize[ring] + 3 ]);

		for( i=1; i<=ringsize[ ring ]; i++ ) {
			found_ring_atom = 0;
			for( a in r ){
				if( a.atomname == ringatoms[ i ] ){
					ringcen[ ring ] += a.pos;
					ringpos[ atring + i ] = a.pos;
					// index atoms with string of the type
					//    "resname:tresnum:atname":
					ringaname[ atring + i ] = r.resname + ":" + 
						sprintf( "%d", r.tresnum ) + ":" + a.atomname;
					found_ring_atom = 1; break;
				}
			}
			if( !found_ring_atom ){
				fprintf( stderr, 
					"Unable to locate ring atom %s in residue %s %d\n",
					ringatoms[ i ], r.resname, r.tresnum );
				exit(1);
			}
		}
		size = ringsize[ ring ];
		ringcen[ ring ] /= size;
