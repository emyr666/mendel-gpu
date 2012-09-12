##convert.awk
##run with: gunzip -c FILE.vcf.gz | gawk -f convertvcf.awk 
(NR==1){ 
	if ($1!="##fileformat=VCFv4.1") { 
			print "Error! VCF input file is not version 4.1" > "/dev/stderr";
			exit 
	}
} 

(substr($1,1,1)!="#" && substr($1,2,1)!="#") {

	if ($3==".") 
		$3=$1 "_" $2 ; 
	printf $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5; 
	
	split($9,format,":");
	i=0;
	for (i in format)
		if (format[i]=="GL")
			gl_index=i;
	if (i==0) {
		print "Error! Could not find 'GL' in FORMAT column at line " $NR > "/dev/stderr";
		exit
	}

	for(i=10;i<=NF;i++) { 
		split($i,indiv,":");
		split(indiv[gl_index],genos,","); 
		printf "\t" genos[1] "\t" genos[2] "\t" genos[3];
	} 
	printf "\n";
}
