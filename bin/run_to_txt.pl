#!/usr/bin/perl
use XML::DOM;

my $file = $ARGV[0];
my $parser = XML::DOM::Parser->new();

my $doc = $parser->parsefile($file);

open FILE, ">", "$ARGV[0].txt" or die $! ;

foreach my $curve ($doc->getElementsByTagName('curve'))
{
    $curvename = $curve->getAttribute('id');

    print FILE $curvename .  "_x " . $curvename . "_y\n";

    $elm = $curve->getElementsByTagName('component')->item(0)->getElementsByTagName('xy')->item(0);

    foreach my $child( $elm->getChildNodes() )
    {
	print FILE $child->getData();
    }
}
