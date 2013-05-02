#!/usr/bin/perl -w
use strict;

if (@ARGV<4){
  print "<basefile> <partsize> <header=1|noheader=0> <is_haploid>\n";
  exit(1);
}

my $basefile=shift;
my $max=shift;
my $is_header=shift;
my $is_haploid=shift;
#my $sexfile=shift;

# read in the sexes
#my @sexarr=();
#open(IN,$sexfile)||die "Can't find $sexfile\n";
#<IN>;
#while(my $line=<IN>){
#  my ($seq,$sex)=split(/\t/,$line);
#  push(@sexarr,$sex);
#}
#close(IN);

my $i=0;
my $c=0;
my $header;
if ($is_header){
  $header=<STDIN>;
  chomp($header);
}
while(my $line=<STDIN>){
  chomp($line);
  if ($i % $max==0){
    if ($i){
      close(OUT3);
    }
    ++$c;
    my $outfile3=$basefile.$c;
    open(OUT3,">$outfile3");
    print STDERR "Writing to $outfile3\n";
    if ($is_header){
      print OUT3 "$header\n";
    }
  }
  #my $persons=scalar(@sexarr);
  my $persons=length($line);
  for(my $i=0;$i<$persons;++$i){
    my $c=substr($line,$i,1);
    #my $sex=$sexarr[$i];
    #if ($chr eq "X" && $sex eq "M" && $c eq "2"){
    if ($is_haploid && $c eq "2"){
      $c="3";
    }
    print OUT3 "$c";
  }
  print OUT3 "\n";
  ++$i;
}
close(OUT3);
