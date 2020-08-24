
#my @id = split (".", $ARGV[0]);

#print "$id[0]\n";

my $id = $ARGV[0];
#print $id;
$id =~ s/.selectedvariants.genotype_modified.txt//g;


open (INPUT, "cat $ARGV[0] | grep -v '+' | ") || die $!;

while (<INPUT>) {
	chomp();
	my @array = split /\t/;
	my @variant = split("[:/]",$array[0]);
	my $chr;
	my $pos;
	my $ref;
	my $alt;
	my $length = @variant;
	if ($length eq 4){
		$chr = $variant[0];
		$chr =~ s/chr//;
		$pos = $variant[1];
		$ref = $variant[2];
		$alt = $variant[3];
	} else {
		$chr = $variant[1];
		$pos = $variant[2];
		$ref = $variant[3];
		$alt = $variant[4];

	}
	$ref = "-" if $ref eq ".";
	$alt = "-" if $alt eq "" or $alt eq "*";
	print "$chr\t$pos\t.\t$ref\t$alt\t.\t.\t$id\t";
	shift @array;
	print join ("\t" , @array);
	print "\n";

}

close INPUT;
