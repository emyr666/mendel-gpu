#!/usr/bin/perl -w
use strict;

if (@ARGV<4){
  print "Usage: <pedfile> <SNP def> <famfile> <bimfile>\n";
  exit 1;
}
my $pedfile=shift;
my $deffile=shift;
my $famfile=shift;
my $bimfile=shift;

sub getzero{
  my $text=shift;
  if ($text eq ""){
    return "0";
  }else{
    return $text;
  }
}

open(INPED,$pedfile)||die "Can't open $pedfile\n";
open(OUTPED,">$famfile")||die "Can't open $famfile\n";
while(my $line=<INPED>){
  chomp($line);
  $line=~s/\ +//g;
  my ($ped,$id,$f,$m,$sex)=split(/,/,$line,-1);
  my $sexcode=0;
  if($sex eq "M") {
    $sexcode=1;
  }elsif($sex eq "F"){
    $sexcode=2;
  }
  print OUTPED "$ped\t$id\t".getzero($f)."\t".getzero($m)."\t$sexcode\t-9\n";
}
close(INPED);
close(OUTPED);

open(INDEF,$deffile)||die "Can't open $deffile\n";
open(OUTDEF,">$bimfile")||die "Can't open $bimfile\n";
my $h1=<INDEF>;
my $h2=<INDEF>;
while(my $line=<INDEF>){
  chomp($line);
  $line=~s/\ +//g;
  my ($snp,$chr,$pos)=split(/,/,$line,-1);
  print OUTDEF "$chr\t$snp\t0\t$pos\tA1\tA2\n";
}
close(INDEF);
close(OUTDEF);
