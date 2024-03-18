package Prediction_hn;

############################
# unspectacular constructor
############################

sub new() {
	my $class = shift;
	my $type = ref( $class ) || $class;
	my $self = { 
		"filename" => "",
		"length1" => 0,
		"length2" => 0,
		"from" => 0,
		"to" => 0,

		"matrixHNHB" => [[]],
		"matrixHNHBW" => [[]],
		"matrixHNHBT" => [[]],
		"matrixHNSHIFT" => [[]],
		"matrixRing" => [[]],
		"matrixEL" => [[]],
		"matrixMA" => [[]],
		"matrixRC" => [[]],
		"matrixCHN" => [[]],
		"sum_HNHB" => [],
		"sum_HNHBW" => [],
		"sum_HNHBT" => [],
		"sum_HNSHIFT" => [],
		"sum_Ring" => [],
		"sum_EL" => [],
		"sum_MA" => [],
		"sum_RC" => [],
		"sum_CHN" => [],
		"mean_HNHB" => [],
		"mean_HNHBW" => [],
		"mean_HNHBT" => [],
		"mean_HNSHIFT" => [],
		"mean_Ring" => [],
		"mean_EL" => [],
		"mean_MA" => [],
		"mean_RC" => [],
		"mean_CHN" => [],
		"residue_C" => [],
		"residue_HB" => [],
	};
	return bless $self, $type;
}

#################################
# read data from AMBER snapshots
#################################

sub readFile() {
	print "@_\n";
	my $self = shift;
	my $infile_base = shift;
	my $from = shift;
	my $to=shift;
	$self->{filename} = $infile_base;
	#$self->{count} = $to - $from;
	$self->{from} = $from;
	$self->{to} = $to;
	my $count = 0;
	if ( $infile_base ) {     # filename defined
	for ( my $i = $from; $i <= $to; $i+=100 ) { # loop over all snapshots
		$count++;
		$infile = "$infile_base.$i";
		if ( -e $infile ) {   # file exists
			my $l;

			print "reading $infile...\n";  # read $infile, remove protection groups
			open( INFH, "$infile" );
			open( OUTFH, ">temp.pdb" );
			open( OUTFH1, ">temp1.pdb" );
			while ( <INFH> ) {
				s/^TER//;
				s/ASH/ASP/;
				s/HIP/HIS/;
				s/HID/HIS/;
				s/HIE/HIS/;
				s/CYX/CYS/;
				print OUTFH $_;
				if(/\s*(?:\w+)\s+(?:\d+)\s+(?:[a-z\A-Z\0-9\+]+)\s+([a-z\A-Z\+]+)\s+(?:\d+)\s+(?:[0-9\-\.]+\s+){2}(?:[0-9\-\.]+)/){
					if($1 eq "Na+" || $1 eq "WAT"){
						s/^.*$/END/;
					}
				}
				print OUTFH1 $_;
			}
			close( INFH );
			close( OUTFH );
			close( OUTFH1 );
    
        
	########################
	# get HN classic shifts
	########################    

	system "/local-home/case/shifts/src/shifts '::H' temp1 > temp1.out";
	open( INFH, "temp1.out" ) || die "Serious problems...\n";

	my $r = 0;
	while ( <INFH> ) {  # loop over all lines
		if(/([a-z\0-9]+)\s+(\w+)\s+(\d+)\s+(\w+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)\s+([0-9\-\.]+)/){
			if( $4 eq "H"){
				$r++;
				$self->{residue_C}->[$r] = $3;
				$self->{matrixRing}->[$r]->[$count] = $5;
				$self->{matrixEL}->[$r]->[$count] = $6;
				$self->{matrixMA}->[$r]->[$count] = $7;
				$self->{matrixRC}->[$r]->[$count] = $9;
				$self->{matrixCHN}->[$r]->[$count] = $10;	
			}
		}
	} #while
  	$self->{length1} = $r; 
	close( INFH );
	unlink( "temp1.*" );

	####################
	# get HN HB Effects
	####################

	system "/local-home/case/shifts/src/shifts -qdb -noreslib temp";
	open( INFH, "temp.out1" ) || die "Serious problems...\n";
	$r = 0;
	my $lastres = 0;
    
	$r = 0;
	while( <INFH> ) {
		if(/\s*(\w+\s+)(\w+\s+)(\w+)\s+(\d+)\s+(\w+)\s+([0-9\-\.]+\s+)([0-9\-\.]+\s+)([0-9\-\.]+)/){
			if($3 ne "PRO"){
				$r++;
				$self->{residue_HB}->[$r] = $4;
				$self->{matrixHNHB}->[$r]->[$count] = $6;
				$self->{matrixHNHBW}->[$r]->[$count] = $7;
				$self->{matrixHNHBT}->[$r]->[$count] = $8;
			}
		}
	}
  	$self->{length2} = $r; 
	close( INFH );
	# unlink( "temp.*" );

	############

		} #if ( -e $infile )
	} #for ( my $i = $from; $i <= $to; $i+=100 )
	} #if ( $infile_base )
return 1;
}

#######################################
# write data as table
#######################################

sub writeTable() {
	my $self = shift;
	my $outfile = shift;
	my $from = shift;
	my $to = shift;
 
	$HNHB = "HN_HB";
	$HNHBW = "HN_HBW";
	$HNHBT = "HN_HBT";
	$HNSHIFT = "HN_SHIFT";
	$Ring = "HN_RingC";
	$EL = "HN_EL";
	$MA = "HN_MAniso";
	$RC = "HN_RC";
	$CHN = "HN_Class";

	open( HN_OUTFH, ">$HNHB");
        open( HN_OUTFH1, ">$HNHBW");
        open( HN_OUTFH2, ">$HNHBT");
        open( HN_OUTFH3, ">$Ring");
        open( HN_OUTFH4, ">$EL");
        open( HN_OUTFH5, ">$MA");
        open( HN_OUTFH6, ">$RC");
        open( HN_OUTFH7, ">$CHN");
	open( HN_OUTFH10, ">$HNSHIFT");

        my $count = 0;
        for ( my $i = $from; $i <= $to; $i+=100 ) {
                $count++;
                for ( my $j = 1; $j <= $self->{length1}; $j++ ) {

                	for ( my $k = 1; $k <= $self->{length2}; $k++ ) {

				if( $self->{residue_C}->[$j] == $self->{residue_HB}->[$k] ){

                        		print HN_OUTFH3 $self->{matrixRing}->[$j]->[$count] . " ";
                        		print HN_OUTFH4 $self->{matrixEL}->[$j]->[$count] . " ";
                        		print HN_OUTFH5 $self->{matrixMA}->[$j]->[$count] . " ";
                        		print HN_OUTFH6 $self->{matrixRC}->[$j]->[$count] . " ";
                        		print HN_OUTFH7 $self->{matrixCHN}->[$j]->[$count] . " ";

	                        	print HN_OUTFH $self->{matrixHNHB}->[$k]->[$count] . " ";
       		                 	print HN_OUTFH1 $self->{matrixHNHBW}->[$k]->[$count] . " ";
       		                 	print HN_OUTFH2 $self->{matrixHNHBT}->[$k]->[$count] . " ";

					$self->{matrixHNSHIFT}->[$k]->[$count] = $self->{matrixCHN}->[$j]->[$count] + $self->{matrixHNHBT}->[$k]->[$count]; 
					print HN_OUTFH10 $self->{matrixHNSHIFT}->[$k]->[$count] . " ";

				}

			}
	
		}

                print HN_OUTFH "\n";
       	        print HN_OUTFH1 "\n";
                print HN_OUTFH2 "\n";
          	print HN_OUTFH3 "\n";
    	        print HN_OUTFH4 "\n";
	        print HN_OUTFH5 "\n";
	        print HN_OUTFH6 "\n";
	        print HN_OUTFH7 "\n";
	        print HN_OUTFH10 "\n";

        }

        close ( HN_OUTFH );
        close ( HN_OUTFH1);
        close ( HN_OUTFH2);
        close ( HN_OUTFH3);
        close ( HN_OUTFH4);
        close ( HN_OUTFH5);
        close ( HN_OUTFH6);
        close ( HN_OUTFH7);
        close ( HN_OUTFH10);

 return 1;

}

########
# mean
########

sub mean() {
	my $self = shift;
	my $from = shift;
	my $to = shift;
	
        $HNHB_M = "HN_HB_Mean";
        $HNHBW_M = "HN_HBW_Mean";
        $HNHBT_M = "HN_HBT_Mean";
        $HNSHIFT_M = "HN_SHIFT_Mean";
        $Ring_M = "HN_RingC_Mean";
        $EL_M = "HN_EL_Mean";
        $MA_M = "HN_MAniso_Mean";
        $RC_M = "HN_RC_Mean";
        $CHN_M = "HN_Class_Mean";
        $RES1 = "Residue_Num";

        open( HN_OUTFH, ">$HNHB_M");
        open( HN_OUTFH1, ">$HNHBW_M");
        open( HN_OUTFH2, ">$HNHBT_M");
        open( HN_OUTFH3, ">$Ring_M");
        open( HN_OUTFH4, ">$EL_M");
        open( HN_OUTFH5, ">$MA_M");
        open( HN_OUTFH6, ">$RC_M");
        open( HN_OUTFH7, ">$CHN_M");
        open( HN_OUTFH8, ">$RES1");
        open( HN_OUTFH10, ">$HNSHIFT_M");

	
	for(my $i = 1; $i <= ($self->{length1}); $i++){

        	for(my $k = 1; $k <= ($self->{length2}); $k++){
		
			if( $self->{residue_C}->[$i] == $self->{residue_HB}->[$k] ){

			my $count = 0;
			for(my $j = $from; $j <= $to; $j+=100){
				$count++;
				$self->{sum_CHN}->[$i] += $self->{matrixCHN}->[$i]->[$count];
				$self->{sum_Ring}->[$i] += $self->{matrixRing}->[$i]->[$count];
				$self->{sum_EL}->[$i] += $self->{matrixEL}->[$i]->[$count];
				$self->{sum_MA}->[$i] += $self->{matrixMA}->[$i]->[$count];
				$self->{sum_RC}->[$i] += $self->{matrixRC}->[$i]->[$count];

                        	$self->{sum_HNHB}->[$k] += $self->{matrixHNHB}->[$k]->[$count];
                        	$self->{sum_HNHBW}->[$k] += $self->{matrixHNHBW}->[$k]->[$count];
                        	$self->{sum_HNHBT}->[$k] += $self->{matrixHNHBT}->[$k]->[$count];
                        	$self->{sum_HNSHIFT}->[$k] += $self->{matrixHNSHIFT}->[$k]->[$count];
			}

			$self->{mean_CHN}->[$i] = ($self->{sum_CHN}->[$i])/$count;
			$self->{mean_Ring}->[$i] = ($self->{sum_Ring}->[$i])/$count;
			$self->{mean_EL}->[$i] = ($self->{sum_EL}->[$i])/$count;
			$self->{mean_MA}->[$i] = ($self->{sum_MA}->[$i])/$count;
			$self->{mean_RC}->[$i] = ($self->{sum_RC}->[$i])/$count;

                	print HN_OUTFH3 $self->{mean_Ring}->[$i] . " ";
                	print HN_OUTFH4 $self->{mean_EL}->[$i] . " ";
                	print HN_OUTFH5 $self->{mean_MA}->[$i] . " ";
                	print HN_OUTFH6 $self->{mean_RC}->[$i] . " ";
                	print HN_OUTFH7 $self->{mean_CHN}->[$i] . " ";
			print HN_OUTFH8 $self->{residue_C}->[$i] . " ";
		
                	print HN_OUTFH3 "\n";
                	print HN_OUTFH4 "\n";
                	print HN_OUTFH5 "\n";
                	print HN_OUTFH6 "\n";
                	print HN_OUTFH7 "\n";
                	print HN_OUTFH8 "\n";

                	$self->{mean_HNHB}->[$k] = ($self->{sum_HNHB}->[$k])/$count;
                	$self->{mean_HNHBW}->[$k] = ($self->{sum_HNHBW}->[$k])/$count;
                	$self->{mean_HNHBT}->[$k] = ($self->{sum_HNHBT}->[$k])/$count;
                	$self->{mean_HNSHIFT}->[$k] = ($self->{sum_HNSHIFT}->[$k])/$count;

                	print HN_OUTFH $self->{mean_HNHB}->[$k] . " ";
                	print HN_OUTFH1 $self->{mean_HNHBW}->[$k] . " ";
                	print HN_OUTFH2 $self->{mean_HNHBT}->[$k] . " ";
                	print HN_OUTFH10 $self->{mean_HNSHIFT}->[$k] . " ";

                	print HN_OUTFH "\n";
                	print HN_OUTFH1 "\n";
                	print HN_OUTFH2 "\n";
                	print HN_OUTFH9 "\n";
                	print HN_OUTFH10 "\n";
			
			}
		
		}

        }

        close ( HN_OUTFH );
        close ( HN_OUTFH1);
        close ( HN_OUTFH2);
        close ( HN_OUTFH3);
        close ( HN_OUTFH4);
        close ( HN_OUTFH5);
        close ( HN_OUTFH6);
        close ( HN_OUTFH7);
        close ( HN_OUTFH8);
        close ( HN_OUTFH10);
        return 1;

}
			

#######################################
# functions for access to private data
#######################################

sub matrixHA() {
 my $self = shift;
 return $self->{matrixHA};
}

#sub mean() {
# my $self = shift;
# my $res = shift;
# return $self->{mean}->[$res];
#}

sub dev() {
 my $self = shift;
 my $res = shift;
 return $self->{dev}->[$res];
}

sub error_dist() {
 my $self = shift;
 my $res = shift;
 return $self->{distribution}->[$res];
}

sub length() {
 my $self = shift;
 return $self->{length}
}

sub error_matrixHA() {
 my $self = shift;
 my $res = shift;
 my @slice;
 for ( my $i = $self->{from}; $i < $self->{to}; $i++ ) {
  push @slice, $self->{error_matrixHA}->[$res]->[$i];
 }
 return \@slice;
}

1;
