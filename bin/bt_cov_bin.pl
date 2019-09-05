#!/usr/bin/env perl
use strict;
use Getopt::Long;
#use Data::Dumper;

#count num bases >thresholds (inputted) from 
# a file/pipe of output from bedtools coverage -d -split
#bins the coverage (count) over each bp in each interval
#input is
#	bt coverage output using options -split and -d.
#	options with defauts (see optional args in usage below):
#		bins - comma delim str of bin thresholds
#		input file or piped to stdin
#		output file or pipe to stdout
#		flag to disable printing header
#		(unimplemented) include averages and stddevs
#output is a table with columns:
#	chromosome
#	interval start
#	interval stop
#	length of interval
#	subsequent columns for each bin threshold
#default values
#	bins >15,>30,>45,>60,>90,>120,>150,>200,>300,>400,>500,>1000

##default vals
my $delim = "|";
my $binstr_default = "15,30,45,60,90,120,150,200,300,400,500,1000";
my $binstr;
my $noheader;
my $addcut;
my $bin;#truebin
my $infile;
my $outfile;


##run vals
my $bt_arg = {"i"=>\$infile,
			  "o"=>\$outfile,
			  "cut"=>\$binstr,
			  "noheader"=>\$noheader,
			  "addcut"=>\$addcut,
			  "bin"=>\$bin};

GetOptions($bt_arg,"i=s","o=s","cut|c=s","noheader|nh","addcut","bin|b");
#print Data::Dumper->Dump(${$bt_arg->{"i"}});

#usage
if(-t *STDIN && !defined(${$bt_arg->{"i"}})){#requires either stdin or file or die
	die "usage: bt_cov_bin.pl -i input file -bin thresholds -o output file (default stdout)\n\
			-i input file (default stdin)\
			-o output file (default stdout)\
			-cut,c comma delimited list of thresholds (default 15,30,45,60,90,120,150,200,300,400,500,1000)\
			-noheader,nh print header (default prints)\
			-addcut add threshold coverage in each cell separated by delimiter\
			-bin,b report read frequency between each threshold\
		#input file optional though requires either input on stdin or specified file\
		#default reports all reads over that threshold\
		#use bin option report all read frequencies between the two thresholds\n"; 
}


##input, stdin or inputted file
my $in_fh;
if(defined(${$bt_arg->{"i"}})){
	open($in_fh,${$bt_arg->{"i"}}) || die "Cannot open input file: ".${$bt_arg->{"i"}}."\n";
}
else{
	$in_fh = \*STDIN;
}
my $out_fh;
if(defined(${$bt_arg->{"o"}})){

	open($out_fh,">${$bt_arg->{'o'}}") || die "Cannot open output file: ".${$bt_arg->{"o"}}."\n";
}
else{
	$out_fh = \*STDOUT;
}


##bins, defaults (above) or input string by comma
my @cuts;
if(defined(${$bt_arg->{"cut"}})){
	@cuts = split(/\,/,${$bt_arg->{"cut"}});
}
else{
	@cuts = split(/\,/,$binstr_default);
}
unshift(@cuts,0) unless $cuts[0] == 0;#add 0 bin

##count each bp
my $prev = undef;#previous interval location

my $num_0 = 0;#count of bp under 1st bin val

#my @lines = <$in_fh>;

#my @col = split(/\t/,$lines[0]);
#my $len = abs($col[2]-$col[1]);
my $len;
my $start;
my $stop;
my $name;

my @freqs = (0) x @cuts;

#print header
unless(${$bt_arg->{"noheader"}}){#optional header
	print $out_fh "#chr\tstart\tstop\tname\tlength\t".join("\t",@cuts)."\n";
}

#my @outs;

#for(my $i=0;$i<@lines;$i++){
#	$lines[$i] =~ s/\n$//;
my $lastname = undef;
while(my $line=<$in_fh>){
	#my($chr,$start,$stop,$rel_pos,$cov) = split(/\t/,$lines[$i]);
	#my @col = split(/\t/,$lines[$i]);
	$line =~ s/\n$//;
	my @col = split(/\t/, $line);
	my $chr = $col[0];
	$start = $col[1];
	$stop = $col[2];
	$name = $col[3];
    my $cov = $col[-1];
	my $coord=$chr.":::".$start.":::".$stop;#unique

	#transition to new interval
	#if($prev ne $coord && $i != 0){#next interval
	$len = abs($col[2]-$col[1]) unless defined($prev);
 
	if(defined($prev) && $prev ne $coord){#next interval
		#output line
		my @loc = split(/:::/,$prev);
		
		if(defined($addcut)){#add threshold to frequency sep by delim
			for(my $x=0;$x<@freqs;$x++){
				$freqs[$x] = $cuts[$x].$delim.$freqs[$x];
			}
		}	
			
		my $outstr = join("\t",@loc)."\t".$lastname."\t".$len."\t".join("\t",@freqs)."\n";
		print $out_fh $outstr;
		#push(@outs,$outstr);
		
		$len = abs($stop-$start);

		@freqs = (0) x @cuts;
	}

	#tally each bin
	for(my $j=0;$j<@cuts;$j++){#skip coord cols
		
		if($cov > $cuts[$j]){
			
			if(defined(${$bt_arg->{"bin"}})){#bin between thresholds
				$freqs[$j]++ if ($j == @cuts-1 || $cov <= $cuts[($j+1)]);
			}			
			else{#bin over thresholds
				$freqs[$j]++;
			}
		}
	}
	
	$prev = $coord;
    $lastname = $name;
}

#output final interval
if(defined($addcut)){#add threshold to frequency sep by delim
	for(my $x=0;$x<@freqs;$x++){
		$freqs[$x] = $cuts[$x].$delim.$freqs[$x];
	}
}	

my @loc = split(/:::/,$prev);

my $outstr = join("\t",@loc)."\t".$name."\t".$len."\t".join("\t",@freqs)."\n";
#push(@outs,$outstr);
print $out_fh $outstr;

##output, stdout or output file
#my $out_fh;
#if(defined(${$bt_arg->{"o"}})){

#	open($out_fh,">${$bt_arg->{'o'}}") || die "Cannot open output file: ".${$bt_arg->{"o"}}."\n";
#}
#else{
#	$out_fh = \*STDOUT;
#}

#print header
#unless(${$bt_arg->{"noheader"}}){#optional header
#	print $out_fh "#chr\tstart\tstop\tlength\t".join("\t",@cuts)."\n";
#}


#print output
#foreach my $out (@outs){
#	print $out_fh $out."\n";
#}

close($in_fh) if defined(${$bt_arg->{"o"}});
close($out_fh) if defined(${$bt_arg->{"o"}});

exit;


#numerical sort
sub number { $a<=>$b }

