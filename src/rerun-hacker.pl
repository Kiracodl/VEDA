#!/usr/bin/perl

#purpose: goes through an xml file and converts all of the <default>
#tags to their corresponding <current> tags.  this allows using an old
#driver file as a starting tool file for a new run.  
#optionally also removes all <output> data.

$line1 = "";
$line2 ="";

while ($line2 = <STDIN>)
{
    if ( ($line1 =~ /default>.+<.default/) && ($line2 =~ /current>(.+)<.current/)) {
	print "<default>$1</default>\n";
	print $line2;

	$line2 = "";
    }
    elsif (($line2 =~ /<output>/) && ($#ARGV > -1))
    {
	print;
	print "</input></run>";
	exit; #exits loop
    }
    else
    {
	print $line1;
    }
    
   $line1 = $line2;
}

print $line1;
