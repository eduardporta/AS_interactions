#!/usr/bin/perl

use strict;
use warnings;

#CHECK ALL THE FILES NEEDED ARE THERE

unless ($ARGV[4]) {
	die "Usage:\nperl interaction_consequences_AS.pl FILE_PFAM_ANNOTATIONS FILE_PFAM_INTERACTIONS FILE_NETWORK FILE_CORR_UCSCISO_ENTREZ FILE_WITH_AS_EVENTS\n";
}

#GET THE PFAM DOMAINS FOR EACH ISOFORM

open IN, $ARGV[0] or die "CAN'T OPEN FILE \"$ARGV[0]\"\n";

my %annotations;

while (<IN>) {

	chomp;
	my @array = split ("\t", $_);

	my $iso = shift @array;
	
	foreach my $pfam (@array) {		
		my ($pfam_name, $start, $end) = split ("-", $pfam);
		$annotations{$iso}{$pfam_name} = 1;
	}
}

close IN or die;

unless (scalar (keys %annotations) > 0) {
	die "NO ANNOTATIONS IN \"$ARGV[0]\"\nIS IT IN THE RIGHT FORMAT?\n";
}

#STORE DOMAIN-DOMAIN INTERACTIONS

my %ddi;

open IN, $ARGV[1] or die "CAN'T OPEN DDI FILE \"$ARGV[1]\"\n";

while (<IN>) {
	
	chomp;
	my ($domain1, $domain2, $source) = split ("\t", $_);
	
	$ddi{$domain1}{$domain2} = 1;
	$ddi{$domain2}{$domain1} = 1;
}

close IN or die;

unless (scalar (keys %ddi) > 0) {
	die "NO DOMAIN-DOMAIN INTERACTIONS IN \"$ARGV[1]\"\nIS IT IN THE RIGHT FORMAT?\n"
}

#STORE PPI
	
my %ppi;

open IN, $ARGV[2] or die "CAN'T OPEN THE INTERACTOME FILE \"$ARGV[2]\"\n";
	
while (<IN>) {	
	
	chomp;
	my ($entrez1, $entrez2) = split ("\t", $_);
		
	$ppi{$entrez1}{$entrez2} = 1;
	$ppi{$entrez2}{$entrez1} = 1;
}

close IN or die;

#STORE ENTREZ-UCSC

my %corr;
my %names;

open IN, $ARGV[3] or die "CAN'T OPEN FILE \"$ARGV[3]\"\n";

while (<IN>) {
	
	chomp;
	my ($iso, $entrez, $symbol) = split ("\t", $_);
	
	$corr{$entrez}{$iso} = 1;
	$names{$entrez} = $symbol;
}

close IN or die;

#ANALYZE SWITCHES

open IN, $ARGV[4] or die "CAN'T OPEN THE FILE WITH THE SWITCHES \"$ARGV[4]\"\n";

while (<IN>) {

	chomp;
	my ($entrez, $symbol, $iso1, $iso2) = split ("\t", $_);

	#FIRST WE CHECK IF THERE ARE INTERACTIONS FOR THIS GENE
		
	unless (exists $ppi{$entrez}) {
		next;
	}
	
	#IF THERE ARE, WE GET THE DOMAINS FOR THE TWO ISOFORMS
			
	my $ref_dom1 = $annotations{$iso1}; 
	my $ref_dom2 = $annotations{$iso2}; 
	
	my @domains1 = keys %$ref_dom1; #GET THE DOMAINS FOR THE FIRST ISOFORM
	my @domains2 = keys %$ref_dom2; #GET THE DOMAINS FOR THE SECOND ISOFORM
			
	#NOW WE GET THE GENES THAT INTERACT WITH OUR QUERY
	
	my $ref_interactions = $ppi{$entrez};
	
	while (my ($entrez2, $ref_source_int) = each %$ref_interactions) {
		
		#WE LOOK FOR THE PROTEIN ISOFORMS OF THE SECOND GENE
		
		unless (exists $corr{$entrez2}) {
			next;
		}
		
		my $ref_iso3 = $corr{$entrez2};

		#WE PRINT THIS LINE JUST TO SHOW THAT WE ANALYZED THIS INTERACTION
		#IF WE HAVE SOME MORE RESULTS ABOUT THE CONSEQUENCES OF THE SPLICING EVENT FOR THIS INTERACTION
		#THEY WILL BE PRINTED IN THIS SAME LINE
		#IF WE CAN'T FIND ANYTHING ELSE, THIS LINE WILL SERVE AS A RECORD THAT, AT LEAST, WE LOOKED AT
		#THIS INTERACTION
		
		print $_, "\t", $entrez2, "\t", $names{$entrez2};
		
		#NOW FOR EACH ISOFORM OF THIS SECOND GENE (ISO3) WE CHECK IF THERE ARE INTERACTIONS LOST
			
		while (my ($iso3, $x) = each %$ref_iso3) {
					
			my $ref_dom3 = $annotations{$iso3};
			my @domains3 = keys %$ref_dom3;
			
			my @ddi1;
			
			#FIRST WE CHECK IF WE CAN FIND AN INTERACTION WITH THE ORIGINAL ISOFORM (ISO1)
			
			foreach my $pfam1 (@domains1) {
				foreach my $pfam3 (@domains3) {
					if (exists $ddi{$pfam1}{$pfam3}) {
						push (@ddi1, $pfam1."/".$pfam3);
					}
				}
			}
			
			#IF THE FIRST ISOFORM (ISO1) COULDN'T INTERACT WITH THIS GENE (ISO3) TO BEGIN WITH
			#IT'S NOT POSSIBLE FOR THE INTERACTION TO BE LOST, SO WE MOVE TO THE NEXT
			
			unless (scalar @ddi1 > 0) {
				next;
			}
			
			#ON THE OTHER HAND, IF WE HAVE A MATCH (I.E. A COMBINATION OF TWO DOMAINS THAT COULD INTERACT)
			#WE PRINT THIS TO KEEP A RECORD OF THE MATCH AND WE ANALYZE FURTHER
			
			print "\tDDI_match";
			
			#NOW WE LOOK AT WHETHER THERE'S A DOMAIN IN THE SECOND ISOFORM (ISO2) OF THE FIRST
			#GENE THAT COULD INTERACT WITH THIS ISOFORM OF THE SECOND GENE (ISO3)
			
			my @ddi2;
			
			foreach my $pfam2 (@domains2) {
				foreach my $pfam3 (@domains3) {
					if (exists $ddi{$pfam2}{$pfam3}) {
						push (@ddi2, $pfam2."/".$pfam3);
					}
				}
			}
			
			#IF THERE'S AT LEAST ONE POSSIBLE DDI COMBINATION, WE SAY THAT THE INTERACTION IS KEPT
			
			if (scalar @ddi2 > 0) {
				my $ddi = join ("_", @ddi2);
				print "\tKept-", $iso3, "-", $ddi;
			}
			
			#IF THERE ISN'T A DDI COMBINATION ANYMORE, THE INTERACTION IS LOST
				
			else {
				my $ddi = join ("_", @ddi1);
				print "\tLost-",$iso3, "-", $ddi;
			}
		}
			
		#NOW WE LOOK FOR GAINS
		#SAME PROCESS, WE ITERATE THROUGH EACH ISOFORM OF THE SECOND GENE (ISO3s)
			
		while (my ($iso3, $x) = each %$ref_iso3) {
					
			my $ref_dom3 = $annotations{$iso3};
			my @domains3 = keys %$ref_dom3;
			
			my @ddi2;
			
			#THIS IS EXACTLY THE SAME AS BEFORE, EXCEPT THAT NOW WE REVERSE THE PROCESS
			#WE START BY LOOKING AT WHETHER THE INTERACTIONS CAN BE MAPPED AT ISO2
			
			foreach my $pfam2 (@domains2) {
				foreach my $pfam3 (@domains3) {
					if (exists $ddi{$pfam2}{$pfam3}) {
						push (@ddi2, $pfam2."/".$pfam3);
					}
				}
			}
			
			#IF THEY CANNOT BE MAPPED TO ISO2, THEN IT'S IMPOSSIBLE THAT THERE'S AN INTERACTION GAIN
			
			unless (scalar @ddi2 > 0) {
				next;
			}
			
			print "\tDDI_match";
			
			#NOW WE CHECK WHETHER THE FIRST ISOFORM COULD INTERACT WITH THIS GENE
			
			my @ddi1;
			
			foreach my $pfam1 (@domains1) {
				foreach my $pfam3 (@domains3) {
					if (exists $ddi{$pfam1}{$pfam3}) {
						push (@ddi1, $pfam1."/".$pfam3);
					}
				}
			}
			
			#IF THE FIRST ISOFORM COULD ALSO INTERACT WITH THIS GENE, THEN THE INTERACTION IS KEPT
			
			if (scalar @ddi1 > 0) {
				my $ddi = join ("_", @ddi1);
				print "\tKept-", $iso3, "-", $ddi;
			}
			
			#IF THE FIRST ISOFORM, ON THE OTHER HAND, COULDN'T INTERACT, THEN THERE'S AN INTERACTION GAIN
				
			else {
				my $ddi = join ("_", @ddi2);
				print "\tGained-",$iso3, "-", $ddi;
			}
		}
	
		print "\n";
	}
}
	
close IN or die;
					
				
					
