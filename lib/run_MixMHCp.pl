#!usr/bin/perl

############
# Written by David Gfeller
#
# For any question, please contact david.gfeller@unil.ch
#
# To cite MixMHCp (version 2.1), please refer to Bassani-Sternberg M and Gfeller D*, J. Immunol. (2016) and Solleder M et al. (2018)
#
# MixMHCp can be used freely by academic groups for non-commercial purposes (see license).
# The product is provided free of charge, and, therefore, on an "as is"
# basis, without warranty of any kind.
#
# FOR-PROFIT USERS
# If you plan to use MixMHCp (version 2.1) in any for-profit
# application, you are required to obtain a separate  license.
# To do so, please contact eauffarth@licr.org or lfoit@licr.org at the Ludwig Institute for  Cancer Research Ltd.
#
# Copyright (2016) David Gfeller
############


use List::MoreUtils qw(uniq);

use Getopt::Long;
use strict;
use File::Copy;


my ($verbose, $output, $input, $dir, $outdir, $temp_keep, $minncomp, $maxncomp, $logo, $comp, $bs, $bias, $name, $logo_type, $lcore, $trash, $alphabet);  

my $MixMHCp_dir = $ARGV[0];

## Take options from commandline
GetOptions ("i=s" => \$input,    # input file name
            "o=s" => \$outdir,    # output dir
	    "b=s" => \$bias,	 # residue bias list file
	    "l=s" => \$logo,   # Decide if the logos should be drawn
            "m1=i" => \$minncomp,# minimal number of PWMs
            "m2=i" => \$maxncomp,# maximal number of PWMs
            "lt=s" => \$logo_type,# type of logo
            "tm" => \$temp_keep, # Don't delete temporary files
            "n=s" => \$name,
	    "lc=i" => \$lcore,
	    "tr=i" => \$trash,
	    "al=s" => \$alphabet,
	    "v"  => \$verbose);  # verbose

if($maxncomp>50){
    $maxncomp=50;
    print "Maximal number of motifs is set to 50\n";
}
if($maxncomp<=0){
    print "Maximal number of motifs must be bigger then 0"."\n";
    exit(10);
}
if($minncomp<=0){
    print "Minimal number of motifs must be bigger then 0"."\n";
    exit(10);
}

if($minncomp > $maxncomp){
    print "Minimal number of motifs must be smaller or equal to maximal number of motifs\n";
    exit(10);
}

if($lcore<2){
    print "Minimal length for the core is 2"."\n";
    exit(10);
}

#Check the alphabet
&check_alphabet($alphabet);



my $naa_min=100;
my $naa_max=0;

my $run_all=1; # Set to 0 if you do not want to run the C++ code.

#Check if the filename contains bad characters
if (!($input eq "")){
    $input = &check_filename($input);
} else {
    print "ERROR!"."\n";
    exit(10);
}


## Check whether neither input nor directory are specified
if (($input eq "") ){
	print "Missing input (\"-i\")\n";
	exit;
}

if($run_all==1){
my @remove_files=qw(project.txt pipeline.log EM_project.txt logos.html length_distribution.html);
my @remove_dir=qw(alignment Multiple_PWMs responsibility logos_html Seq2Logo LoLa ggseqlogo logos weights KLD);

my $f;
foreach $f (@remove_files){
    if(-e "$outdir/$f"){
	system("rm $outdir/$f");
    }
}
foreach $f (@remove_dir){
    if(-d "$outdir/$f/"){
	system("rm -r $outdir/$f/");
    }
}
}

$bs=0;
if ($bias ne "U"){

    &check_bias($bias, $alphabet);
    system("cp ".$bias." $outdir"."/bias.txt");
    $bs=2;
    
} elsif ($bias eq "U"){
    $bs=1;
    if($alphabet ne "ACDEFGHIKLMNPQRSTVWY"){
	print "Impossible to use Uniprot background frequencies (default) with non amino acid alphabet\n";
	exit;
    }
   
}

######## Here we will need to check the bias files, to make sure they are consistent with the alphabet #############
	    
system("mkdir -p ".$outdir);
system("mkdir -p ".$outdir."/data");
system("mkdir -p ".$outdir."/KLD");
system("mkdir -p ".$outdir."/responsibility");
system("mkdir -p ".$outdir."/Multiple_PWMs");
system("mkdir -p ".$outdir."/weights");
system("mkdir -p ".$outdir."/weights/plots");
#system("mkdir -p ".$outdir."/LoLa");
system("mkdir -p ".$outdir."/ggseqlogo");
if($logo==1){
    system("mkdir -p ".$outdir."/logos");
    system("mkdir -p ".$outdir."/logos_html");
}


my $ct;
my $p;

# Read input
print_and_log("Processing input files...\n");
my @pep=();
open IN, $input, or die;
while(my $l=<IN>){
    if(substr($l, 0, 1) ne ">"){
	$l =~ s/\r?\n$//;
        push @pep, $l;
    }
}
close IN;


@pep=uniq(@pep);
@pep=sort(@pep);

#Check input
my $exit_pep=check_input(@pep);
if($exit_pep ne ""){
    exit;
}

open OUT, ">$outdir/data/peptides.fa";
my @len=();
$ct=1;
foreach $p (@pep){
   printf OUT ">$ct %d\n$p\n", length($p);
   push @len, length($p);
#     printf OUT ">$ct\n$p\n";
    $ct++;
}
close OUT;

@len=uniq(@len);
@len=sort(@len);

#die;
if( !(grep { $_ eq $lcore } @len) ){
	print "Error: no peptides of core length in input..."."\n";
	print "Core length: ".$lcore."\n";
	print "Give a different core length [-lc] or double check data."."\n";
	die;
}

if($run_all==1){
    my $command = " -m1 ".$minncomp." -m2 ".$maxncomp." -d ".$outdir." -b $bs -a ".$alphabet." -lc ".$lcore." -tr ".$trash;
    print_and_log("Running MixMHCp..."."\n");
    #print_and_log($MixMHCp_dir."MixMHCp.x $command >> ".$outdir."pipeline.log 2>> ".$outdir."pipeline.log\n");
    my $exit_status = system($MixMHCp_dir."MixMHCp.x $command >> ".$outdir."pipeline.log 2>> ".$outdir."pipeline.log");
    
    if (!($exit_status == 0)){
	print "Error: MixMHCp failed to execute."."\n";
	exit($exit_status);
    }
}


#my $logo="Seq2Logo";  #This is not working well on the GUI because the path to Seq2Logo and gs cannot be found.
#Moreover, in case of gaps, Seq2Logo.py treats them as no sequences... so this works only well if we do not align the sequences

if($logo==1 && $run_all==1){

    my $input_type = "Protein";
        
    if($logo_type eq "ggseqlogo"){
	my $inputType = 'pwm'; # 'resp'; 
	print_and_log("Generating logos with ggseqlogo... "."\n");
	#print("Rscript $MixMHCp_dir/calling-ggseqlogoMOD.r $MixMHCp_dir $outdir $outdir $alphabet $minncomp $maxncomp $inputType\n");
	system("Rscript $MixMHCp_dir/calling-ggseqlogoMOD.r $MixMHCp_dir $outdir $outdir $alphabet $minncomp $maxncomp $inputType >> $outdir/pipeline.log");    #Need to fix $minncomp.
	#die;

    #} elsif ($logo_type eq "LoLa"){
	
	#print_and_log("Generating logos..."."\n");
	#system("mkdir -p ".$outdir."logos");
	#system("java -Xmx1024m -cp ".$MixMHCp_dir."jar/aida-3.3.jar:".$MixMHCp_dir."jar/biojava-1.4a.jar:".$MixMHCp_dir."jar/brainlib-1.4.jar:".$MixMHCp_dir."jar/freehep-export-2.1.1.jar:".$MixMHCp_dir."jar/freehep-graphics2d-2.1.1.jar:".$MixMHCp_dir."jar/freehep-graphicsio-2.1.1.jar:".$MixMHCp_dir."jar/freehep-graphicsio-emf-2.1.1.jar:".$MixMHCp_dir."jar/freehep-graphicsio-java-2.1.1.jar:".$MixMHCp_dir."jar/freehep-graphicsio-pdf-2.1.1.jar:".$MixMHCp_dir."jar/freehep-graphicsio-ps-2.1.1.jar:".$MixMHCp_dir."jar/freehep-graphicsio-svg-2.1.1.jar:".$MixMHCp_dir."jar/freehep-graphicsio-swf-2.1.1.jar:".$MixMHCp_dir."jar/freehep-graphicsio-tests-2.1.1.jar:".$MixMHCp_dir."jar/freehep-io-2.0.2.jar:".$MixMHCp_dir."jar/freehep-swing-2.0.3.jar:".$MixMHCp_dir."jar/freehep-util-2.0.2.jar:".$MixMHCp_dir."jar/freehep-xml-2.1.MixMHCp.jar:".$MixMHCp_dir."jar/itext-2.0.2.jar:".$MixMHCp_dir."jar/jas-plotter-2.2.jar:".$MixMHCp_dir."jar/jdom-1.0.jar:".$MixMHCp_dir."jar/junit-3.8.2.jar:".$MixMHCp_dir."jar/openide-lookup-1.9-patched-1.0.jar:".$MixMHCp_dir." CreateLogo $outdir/ $input_type png $outdir F >> ".$outdir."pipeline.log");

	
    } elsif ($logo_type eq "Seq2Logo"){

	print_and_log("Generating Seq2Logo... "."\n");
	system("mkdir -p ".$outdir."/Seq2Logo/");
	#The parameters for the Logo (background frequency, pseudo-count, Shannon/KL, format... are set in Seq2Logo.pl
	system("perl ".$MixMHCp_dir."Seq2Logo.pl $outdir  $MixMHCp_dir $naa_min $naa_max $trash >> ".$outdir."pipeline.log 2>> ".$outdir."pipeline.log");
	
    }
    
    print_and_log("Generating HTML tables document..."."\n");
    system("perl ".$MixMHCp_dir."tablegenerator_html.pl $outdir $logo_type $name $naa_min $naa_max $minncomp $maxncomp $trash $lcore"); #Need to fix $minncomp.
    
}

#Plot the peptide length distributions
my $exit_status = system("Rscript $MixMHCp_dir/plot_length.R $outdir $naa_min $naa_max $minncomp $maxncomp $trash"); #Need to fix $minncomp.
if (!($exit_status == 0)){
	print "Error: Rscript failed to plot the peptide length distributions."."\n";
	exit($exit_status);
    }


if ($temp_keep == 0){
   
    #system("rm -r $outdir/LoLa");
    #system("rm -r $outdir/ggseqlogo");
    if ($logo_type eq "Seq2Logo"){
    	system("rm -r $outdir/Seq2Logo");
    }
   
}

print "Done"."\n";
exit;

sub check_filename{
	my $file = $_[0];
	my $new_name;
	
	my @bad_char=('\(', '\)', '\|', '\"', '\;');
	my $c;
	$new_name=$file;
	
	foreach $c (@bad_char){
	    $new_name =~ s/$c/_/g;
	}
	#Replace "\"
	$new_name =~ s/\\/_/g;
	
	if($new_name ne $file){
	    print_and_log("Sorry, we had to rename your file \"$file\" into \"$new_name\"\n");
	    copy($file, $new_name);
	    $file=$new_name;
	}
	return($new_name);
}

sub print_and_log{
    print "$_[0]";
    system("echo \"$_[0]\" >> $outdir/pipeline.log");
}

sub check_input{


    my %aa=();
    my @al=split('', $alphabet);
    my $s;
    foreach $s (@al){
	$aa{$s}=1;
    }
    
    my @pep=@_;
    
    my $p;
    my $le=length($pep[0]);
    my @a=();
    my $s;
    my $exit_pep="";
    my $lp;

    my $min=100;
    my $max=0;
    
    foreach $p (@pep){
	$lp=length($p);
	#if(length($p) != $le){
	#    print_and_log("Peptides of different lengths: $le-mers and $lp-mers\n");
	#    $exit_pep=$p;
	#    last;
	#}
	if(length($p) < 5 ){
	    print_and_log("Peptides should be at least of length 5\n");
	    $exit_pep=$p;
	    last;
	}
	if(length($p) > $lcore+10 ){
	    print_and_log("Peptides should not be longer than 10 amino acids + length core ($lcore)\n");
	    $exit_pep=$p;
	    last;
	}
	if($naa_min>length($p)){
	    $naa_min=length($p);
	}
	if($naa_max<length($p)){
	    $naa_max=length($p);
	}

	
	@a=split('', $p);
	foreach $s (@a){
	    if(!exists $aa{$s}){
		print_and_log("Unknown letter in alphabet: $s in $p\n");
		$exit_pep=$s;
		last;
	    }
	}
    }
    return($exit_pep);
}

sub check_alphabet{

    my $s;
    
    my @allowed=qw(A C D E F G H I K L M N P Q R S T V W Y a b c d e f g h i j k l m n o p q r s t u v w x y z);
    my %allowed_alphabet=();
    foreach $s (@allowed){
	$allowed_alphabet{$s}=1;
    }
    
    my @al=split('', $_[0]);
   
    
    foreach $s (@al){
	if(!exists $allowed_alphabet{$s}){
	    print_and_log("Letter not allowed in alphabet: $s\n");
	    die;
	}
    }
}

sub check_bias{

    my $l;
    my $c=0;
    my @b;
    my @ba;
    open IN, "$_[0]", or die "Bias file $_[0] not found\n";
    while($l=<IN>){
	$l =~ s/\r?\n$//;
	@b=split(' ', $l);
	if($b[1] =~ /^[0-9,.E]+$/ ){
	    if($b[1]<=0){
		print_and_log("Invalid bias values: $l\n");
		exit
	    }
	} else {
	    print_and_log("Invalid bias values: $l\n");
	    exit
	}
	$c=$c+$b[1];
	push @ba, $b[0];
    }
    close IN;

    if($c <0.999 || $c > 1.001){
	print_and_log("Bias values will be automatically normalized to one.\n");
    }
    
    my @a=split('', $_[1]);

    @a=sort @a;
    @ba = sort @ba;

    #print "@a\n@ba\n";

    my $m;
    if(scalar(@a) != scalar @ba){
	print_and_log("Invalid bias file (not the same entries as the alphabet).\n");
	exit;
    }
    
    for(my $i=0; $i<scalar(@a); $i++){
	if($a[$i] ne $ba[$i]){
	    print_and_log("Invalid bias file (not the same entries as the alphabet).\n");
	    exit;
	}
    }
    
}
