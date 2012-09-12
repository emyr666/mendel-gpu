#!/usr/bin/perl -w
use strict;

my $epsilon=.001;
my $log_epsilon=log($epsilon);

if (@ARGV<2){
  print STDERR "Usage: <persons><snps>\n";
  exit(1);
}
sub recode{
  my $val=shift;
  return $val;
  #$val=$val<$epsilon?$log_epsilon:log($val);
  return $val;
}

my $n=shift;
my $p=shift;
#print "$n\t$p\n";
while(my $line=<STDIN>){
  chomp($line);
  my @linearr=split(/\t/,$line);
  my $len=scalar(@linearr);
  for(my $i=0;$i<$len;$i+=4){
    my $hap1a=($linearr[$i]+$linearr[$i+1]);
    my $hap1b=($linearr[$i+2]+$linearr[$i+3]);
    my $hap2a=($linearr[$i]+$linearr[$i+2]);
    my $hap2b=($linearr[$i+1]+$linearr[$i+3]);
    #my $sum = $hap1+$hap2;
    #$hap1/=$sum;$hap2/=$sum;
    #$hap1=$hap1<$epsilon?$log_epsilon:log($hap1);
    #$hap2=$hap2<$epsilon?$log_epsilon:log($hap2);
    if ($i){print "\t";}
    print recode($hap1a)."\t".recode($hap1b)."\t".recode($hap2a)."\t".recode($hap2b);
;
  }
  print "\n";
}
