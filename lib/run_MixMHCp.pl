#!/usr/bin/perl

############
# Written by David Gfeller
#
# For any question, please contact david.gfeller@unil.ch
#
# To cite MixMHCp, please refer to Bassani-Sternberg M and Gfeller D*, J. Immunol. (2016)
#
# MixMHCp can be used freely by academic groups for non-commercial purposes (see license).
# The product is provided free of charge, and, therefore, on an "as is"
# basis, without warranty of any kind.
#
# FOR-PROFIT USERS
# If you plan to use MixMHCp (version 1.0) in any for-profit
# application, you are required to obtain a separate  license.
# To do so, please contact eauffarth@licr.org or lfoit@licr.org at the Ludwig Institute for  Cancer Research Ltd.
#
# Copyright (2016) David Gfeller
############


use List::MoreUtils qw(uniq);

use Getopt::Long;
use strict;
use File::Copy;


my ($verbose, $output, $input, $dir, $outdir, $temp_keep, $help, $maxncomp, $logo, $comp, $bs, $bias, $name, $logo_type);  

my $MixMHCp_dir = $ARGV[0];

## Take options from commandline
GetOptions ("i=s" => \$input,    # input file name
            "o=s" => \$outdir,    # output dir
	    "b=s" => \$bias,	 # residue bias list file
	    "l=s" => \$logo,   # Decide if the logos should be drawn
            "m=i" => \$maxncomp,# maximal number of PWMs
            "lt=s" => \$logo_type,# type of logo
            "h" => \$help,  	 # Help
            "tm" => \$temp_keep, # Don't delete temporary files
            "n=s" => \$name,
	    "v"  => \$verbose);  # verbose

if($maxncomp>20){
    $maxncomp=20;
}

#Check if the filename contains bad characters
if (!($input eq "")){
    $input = &check_filename($input);
} else {
    print "ERROR!"."\n";
    exit(10);
}

## Print help page
if ($help == 1){
	system("cat ".$MixMHCp_dir."help_page.txt");
	exit;
}

## Check whether neither input nor directory are specified
if (($input eq "") ){
	print "Missing input (\"-i\")\n";
	exit;
}


my @remove_files=qw(project.txt pipeline.log EM_project.txt blosum62.txt logos.html);
my @remove_dir=qw(alignment KLD Multiple_PWMs responsibility logos_html Seq2Logo LoLa logos);

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

system("cp $MixMHCp_dir/blosum62.txt $outdir"."/blosum62.txt");

$bs=0;
if (-e $bias){
    system("cp ".$bias." $outdir"."/bias.txt");
    $bs=2;
} elsif ($bias eq "U"){
    $bs=1;
}

	    
system("mkdir -p ".$outdir);
system("mkdir -p ".$outdir."/data");
system("mkdir -p ".$outdir."/KLD");
system("mkdir -p ".$outdir."/responsibility");
system("mkdir -p ".$outdir."/Multiple_PWMs");
system("mkdir -p ".$outdir."/logos");
system("mkdir -p ".$outdir."/LoLa");

my $ct;
my $p;

# Read input
print_and_log("Processing input files...\n");
my @pep=();
open IN, $input, or die;
while(my $l=<IN>){
    if(substr($l, 0, 1) ne ">"){
	print "Invalid input format (not fasta file): $l\n";
    }
    $l=<IN>;
    $l =~ s/\r?\n$//;
    #chomp($l);
    push @pep, $l;
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
$ct=1;
foreach $p (@pep){
    print OUT ">$ct\n$p\n";
    $ct++;
}
close OUT;


my $command = "0 0 ".$maxncomp." -d ".$outdir." -b U";

print_and_log("Running MixMHCp..."."\n");
my $exit_status = system($MixMHCp_dir."MixMHCp.x $command >> ".$outdir."pipeline.log 2>> ".$outdir."pipeline.log");
if (!($exit_status == 0)){
    print "Error: MixMHCp failed to execute."."\n";
    exit($exit_status);
}



#my $logo="Seq2Logo";  #This is not working well on the GUI because the path to Seq2Logo and gs cannot be found.
#Moreover, in case of gaps, Seq2Logo.py treats them as no sequences... so this works only well if we do not align the sequences

if($logo==1){

   
    my $input_type = "Protein";
    if($logo_type eq "LoLa"){
	
	print_and_log("Generating logos..."."\n");
	system("mkdir -p ".$outdir."logos");
	system("java -Xmx1024m -cp ".$MixMHCp_dir."jar/aida-3.3.jar:".$MixMHCp_dir."jar/biojava-1.4a.jar:".$MixMHCp_dir."jar/brainlib-1.4.jar:".$MixMHCp_dir."jar/freehep-export-2.1.1.jar:".$MixMHCp_dir."jar/freehep-graphics2d-2.1.1.jar:".$MixMHCp_dir."jar/freehep-graphicsio-2.1.1.jar:".$MixMHCp_dir."jar/freehep-graphicsio-emf-2.1.1.jar:".$MixMHCp_dir."jar/freehep-graphicsio-java-2.1.1.jar:".$MixMHCp_dir."jar/freehep-graphicsio-pdf-2.1.1.jar:".$MixMHCp_dir."jar/freehep-graphicsio-ps-2.1.1.jar:".$MixMHCp_dir."jar/freehep-graphicsio-svg-2.1.1.jar:".$MixMHCp_dir."jar/freehep-graphicsio-swf-2.1.1.jar:".$MixMHCp_dir."jar/freehep-graphicsio-tests-2.1.1.jar:".$MixMHCp_dir."jar/freehep-io-2.0.2.jar:".$MixMHCp_dir."jar/freehep-swing-2.0.3.jar:".$MixMHCp_dir."jar/freehep-util-2.0.2.jar:".$MixMHCp_dir."jar/freehep-xml-2.1.MixMHCp.jar:".$MixMHCp_dir."jar/itext-2.0.2.jar:".$MixMHCp_dir."jar/jas-plotter-2.2.jar:".$MixMHCp_dir."jar/jdom-1.0.jar:".$MixMHCp_dir."jar/junit-3.8.2.jar:".$MixMHCp_dir."jar/openide-lookup-1.9-patched-1.0.jar:".$MixMHCp_dir." CreateLogo $outdir/ $input_type png $outdir F >> ".$outdir."pipeline.log");

	
    } elsif ($logo_type eq "Seq2Logo"){

	print_and_log("Generating Seq2Logo... "."\n");
	system("mkdir -p ".$outdir."/Seq2Logo/");
	#The parameters for the Logo (background frequency, pseudo-count, Shannon/KL, format... are set in Seq2Logo.pl
	system("perl ".$MixMHCp_dir."Seq2Logo.pl $outdir  $MixMHCp_dir >> ".$outdir."pipeline.log 2>> ".$outdir."pipeline.log");
	
    }
    
    print_and_log("Generating HTML tables document..."."\n");
    system("perl ".$MixMHCp_dir."tablegenerator_html.pl $outdir $logo_type $name");
    
}

system("rm $outdir/blosum62.txt");


if ($temp_keep == 0){
   system("rm $outdir/project.txt");
   system("rm $outdir/EM_project.txt");
   system("rm -r $outdir/LoLa");
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


    my %aa=("A", 1, "C", 2, "D", 3, "E", 4, "F", 5, "G", 6, "H", 7, "I", 8, "K", 9, "L", 10, "M", 11, "N", 12, "P", 13, "Q", 14, "R", 15, "S", 16, "T", 17, "V", 18, "W", 19, "Y", 20);

    my @pep=@_;
    
    my $p;
    my $le=length($pep[0]);
    my @a=();
    my $s;
    my $exit_pep="";
    my $lp;
    
    foreach $p (@pep){
	$lp=length($p);
	if(length($p) != $le){
	    print_and_log("Peptides of different lengths: $le-mers and $lp-mers\n");
	    $exit_pep=$p;
	    last;
	}
	@a=split('', $p);
	foreach $s (@a){
	    if(!exists $aa{$s}){
		print_and_log("Unknown amino acid $s in $p\n");
		$exit_pep=$s;
		last;
	    }
	}
    }
    return($exit_pep);
}
