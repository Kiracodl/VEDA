#!/usr/bin/perl
#this script will run VEDA given an input filename, and rename the randomly generated output file
#to be a deterministic filename

$fn = $ARGV[0];

print $fn;

open(PIPE, "nice az_ddaskr $fn |") || die("didnt open veda");

while ($_ = <PIPE>)
{ 
   if ( $_ =~ /.RAPPTURE-RUN..(run.+xml)/)       
   {
       system("mv $1 $fn.out");
   }
   else
   {
       print $_;
   }
}
 
