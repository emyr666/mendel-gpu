#!/usr/bin/perl -w
use strict;

if (@ARGV<5){
  print "Usage: <mincounts> <legend in> <hap in> <legend out> <hap out>\n";
  exit(1);
}
my $mincounts=shift;
my $inlegend=shift;
my $inhap=shift;
my $outlegend=shift;
my $outhap=shift;
open(OUTLEGEND,">$outlegend")||die "Can't open $outlegend\n";
open(OUTHAP,">$outhap")||die "Can't open $outhap\n";
open(INLEGEND,"$inlegend")||die "Can't open $inlegend\n";
open(INHAP,"$inhap")||die "Can't open $inhap\n";
my $line1=<INLEGEND>;
print OUTLEGEND "$line1";
while(my $line1=<INLEGEND>){
  my $line2=<INHAP>;
  my $counts=$line2 =~ tr/1//;
  if ($counts>=$mincounts){
    #print STDERR "Counts: $counts\n";
    print OUTLEGEND "$line1";
    print OUTHAP "$line2";
  }
}



