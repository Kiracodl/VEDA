#!/usr/bin/perl

# parameter space exploration tool 
# Daniel Kiracofe, Nov 2008
# purpose: execute a program repeatedly with multiple different input values

# usage:  param-explore.pl RUNSFILE  TEMPLATEFILE EXECUTABLEFILE [ OUTPUTPREFIX  STARTLINE STOPLINE]

# example: param-explore.pl runs.txt  tool.xml ~/bin/az_ddaskr

# RUNSFILE:  This should be a text with each line containing the full parameter set for one run.
# should be a comma separated list of parameter=value.  For example, to make 4 runs with different values of
# two parameters you might use this file:

# PARAM_E = 100, PARAM_nu = 10
# PARAM_E = 100, PARAM_nu = 100
# PARAM_E = 10,  PARAM_nu = 10
# PARAM_E = 10,  PARAM_nu = 100

# the use of the prefix PARAM is not necessary, but recommended to make it easy to recognize what is a parameter
# in your template file.  

# TEMPLATEFILE: The template file is the file that you would ordinarily use as an input to your code,
# except that some of the inputs have been replaced by tags.  For example, if your input file normally looks like this


#            <number id="Etip">
#                <about>
#                    <label> Young's modulus of tip (GPa)</label>
#                    <description/>
#                </about>
#                <current>130</current>
#            </number>

# you might edit it to look like:

#            <number id="Etip">
#                <about>
#                    <label> Young's modulus of tip (GPa)</label>
#                    <description/>
#                </about>
#                <current> PARAM_E </current>
#            </number>

# this param-explore.pl will replace PARAM_E in this file with the values that you have given it in the runs file
# note, this is a literal text substitution.  Make sure that the tag names you pick don't occur anywhere else in the
# file.  E.g. If you pick something like 'o' for your parameter instead of 'PARAM_E', you would a file that looks like this

#            <number id="Etip">
#                <ab100ut>
#                    <label> Y100ung's m100dulus 100f tip (GPa)</label>
#                    <descripti100n/>
#                </ab100ut>
#                <current> PARAM_E </current>
#            </number>

# which is probably not what you want.  Hence the suggestion to use a unique prefix like PARAM_

# also, the script will do simple math on any line where a replacement is made.  e.g. if TAG=Z=10, then
#
# <current>TAG_Z + 10</current>
#
# would be replaced with
#
# <current>20</current>

# EXECUTABLEFILE : the actual program you want to run.  due to bug ("feature") in the program, you might need
# to use an absolute path rather than a relative one (.e.g /home/drk/bin/az_ddaskr instead of ../bin/ddaskr)

# OUTPUTPREFIX : we assume that your code produces only one output file, and that it is placed in the current
# working directory.  These files will be renamed according to OUTPUTPREFIX.  Default is 'output', which
# yields files output1.out output2.out output3.out etc.

#STARTLINE STOPLINE : if these are omitted, then the whole file is executed.  Otherwise, only the lines between STARTLINE and STOPLINE are executed.  This is to allow running from the same runs.txt file on multiple different machines.  If STARTLINE is given but STOPLINE is omitted, STOPLINE is set to 1 million.
#
# 

use strict;
use File::Temp qw/ tempfile tempdir /;

my $runsfile = $ARGV[0];
my $templatefile = $ARGV[1];
my $executable = $ARGV[2];
my ($outputprefix, $newoutput, $startline, $stopline);

if (@ARGV <= 1) { print "usage:  param-explore.pl RUNSFILE  TEMPLATEFILE EXECUTABLEFILE [ OUTPUTPREFIX STARTLINE STOPLINE ]\n"; exit }

if (@ARGV >= 3) { $outputprefix = $ARGV[3] }  else {  $outputprefix = "output"; }

if (@ARGV == 5) {
    #specified start line only
    $startline = $ARGV[4];
    $stopline  = 1000000;
} elsif (@ARGV == 6) {
    $startline = $ARGV[4];
    $stopline  = $ARGV[5];
} else
{
    $startline = 0;
    $stopline  = 1000000;
}

print "startline = " . $startline . " stopline = " . $stopline . "\n";

my($tagline, %tags, @match, $pair, $newdriver, $ndx, $tmpdir);


open( TAGS, $runsfile) || die("could not open runs file");

#fixme, make sure this is a unique directory name

$tmpdir = tempdir( 'param_tmp_XXXX', CLEANUP => 1 )  || die("could not create a temporary working directory");
chdir $tmpdir;

while ( $tagline = <TAGS>)
{
    $ndx++;
    if (( $ndx >= $startline) && ($ndx <= $stopline))
    {
	print "running job # " . $ndx . "\n";
    
	#get tags for this run
	%tags = ();
	@match = split(/,/, $tagline);
	
	foreach $pair ( @match ){
	    $pair =~ /\s*(\S+)\s*=\s*(\S+)\s*/;
	    $tags{$1} = ${2};
	}
	
	#create new driver file
	$newdriver = "driver" . $outputprefix . $ndx . ".txt"; 
	
	tag_replacer( "../" . $templatefile, $newdriver, %tags );
	# run the excutable 
	
        print "running $executable $newdriver \n";
#	system( $executable . " " . $newdriver . ' 2>/dev/null 1>/dev/null');
	system( $executable . " " . $newdriver);
	
	$newoutput = "../" . $outputprefix . $ndx . ".out";
	
	system('mv', $newdriver, '..');
	
	# rename output file and move up one dir
	print "running mv * $newoutput \n";
	system("mv * $newoutput");
#    system("gzip $newoutput"); # for faster downloads
    }
    else
    {
	print "skipping job # " . $ndx . "\n";
    }
}

chdir "..";
rmdir 'tmp';

# function to replace tags
sub tag_replacer {
    my ($infile, $outfile, %tags) = @_;
    my $key;
    my $sum;
    my @keys1;

    open(IN, $infile) || die("could not open template file");
    open(OUT, ">" . $outfile) || die("could not create new driver file");

    while ( $_ = <IN> )
    {
	# if have two tags such that one is a substring of the other (e.g. PARAM_E and PARAM_ETA), we need to
	# try to replace the longer one first. Sort the list by length.  note: this is an inefficient sort b/c it
	#recalculates the length at every comparision
	@keys1 =  sort { length $b <=> length $a }  keys(%tags) ;

	#this is kind of inefficient if we have a lot of tags, but a better way is not coming to mind at the moment
	foreach $key ( @keys1 )
	{
	    if ( $_ =~ /$key/)
	    {
		$_ =~ s/$key/$tags{$key}/g;

		#hack to do simple integer math here.
		if ( $_ =~ /([0-9]+) *\+ *([0-9]+)/ )
		{
		    $sum = $1 + $2;
		    $_ =~ s/([0-9]+) *\+ *([0-9]+)/$sum/;
		}
	    }
	}

	print OUT $_;
    }

    close(IN);
    close(OUT);
}
    
