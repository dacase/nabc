#!/usr/bin/env perl 

#  Parse a "ringfrag" file to determine which rings are in the quantum
#     fragments; the parse the SHIFTS output file (.emp) to compute the ring
#     current terms from the rings in the mm regions.  Put everything into
#     and rdb file (on stdout) that could be used to correct the afnmr
#     output results.

#     Major limitions (for now): only works for nucleic acids; only works
#     for proton shifts; requires a "ringfrag" file that would for now be
#     obtains by grep "C1'" 1bna???.pqr; then editing that file to have
#     two columns, one with the fragment number (from the file name) and one
#     with the residue number of the ring that is the qm fragment pqr file.

#     Note also: use the -details flag in the SHIFTS run to get the required
#     details about individual ring current contributions

sub resatom { @GA = split( ':', $a); 
              @GB = split( ':', $b); 
              $GA[0] <=> $GB[0]  or $GA[1] cmp $GB[1] or $a cmp $b; }

open(RINGFRAG, "ringfrag");

while( <RINGFRAG> ){
	@F = split(' ');
	$ringfrag{$F[0] . ":" . $F[1]} = 1;
}
close(RINGFRAG);

print "res	atomname	mmshhm\n4N	8	10N\n";

while( <> ){
	@F = split(' ');

	$proton = $F[4] if( $F[0] eq "Detailed" && $F[1] eq "contributions");
	$shhm{ $proton . $F[2] } += $F[3] if( $F[0] eq "ring:");
}

foreach $key (sort keys %shhm){
	@G = split( ':', $key);
	$keyp = $key;
	$keyp =~ s/:/\t/g;
	if( exists $ringfrag{$G[1] . ":" . $G[3]} ){
		printf( STDERR  "skipping: %s\t%8.3f\n",  $keyp,$shhm{$key});
	} else {
		$mmshhm{$G[1] . ":" . $G[2]} += $shhm{$key};
	}
}

foreach $key (sort resatom keys %mmshhm){
	$keyp = $key;
	$keyp =~ s/:/\t/g;
	printf( "%s\t%8.3f\n",  $keyp,$mmshhm{$key}) 
}
	
	
