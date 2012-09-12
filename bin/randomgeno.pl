#!/usr/bin/perl -w
use strict;

for(my $j=0;$j<70;++$j){
  for(my $i=0;$i<1092;++$i){
    my $g=rand()*4;
    printf "%d",$g;
  }
  print "\n";
}
