#!/usr/bin/perl -w
use strict;

if (@ARGV<3){
  print "<basefile><partsize><header=1|noheader=0>\n";
  exit(1);
}

my $basefile=shift;
my $max=shift;
my $is_header=shift;
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
  print OUT3 "$line\n";
  ++$i;
}
close(OUT3);
