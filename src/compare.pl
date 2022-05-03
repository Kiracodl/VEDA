#!/usr/bin/perl

if (@ARGV >= 1) { 
    $tol = $ARGV[0];
} else {
    $tol = 1e-4;
}

while ($_ = <STDIN>)
{

    if ( $_ =~ /\s*(\S+)\s+(\S+)\s*\|\s*(\S+)\s+(\S+)/)
    {
	chomp($_);

	#differing numerical line

	$1 != 0 && ($x = (($1 - $3)/$1));
	$2 != 0 && ($y = (($2 - $4)/$2));
	if (( abs($x) > $tol) || ( abs($y) > $tol)) {
	    print  $_  . " || " . $y . "\n" ;
	}

	if (( $2 == "Infinity") != ($4 == "Infinity")) {
	    print  $_  . " || " . $y . "\n" ;
	}

    }
    elsif ( $_ =~ /curve /)
    {
	#curve line
	print  $_;
    }
    elsif ( $_ =~ /.+\|.+/ )
    {
	#differing non-numerical line
	print  $_;
    }
#    else {  #non-differing line, don't print
#	print $_;
#    }

}


