#!/usr/bin/perl

use List::MoreUtils qw(uniq);

use strict;


my $outdir = $ARGV[0];
my $logo_type=$ARGV[1];
my $name=$ARGV[2];

my $logo_dir;
if($logo_type eq "LoLa"){
    $logo_dir="logos";
} elsif($logo_type eq "Seq2Logo"){
    $logo_dir="Seq2Logo";
}
my $pwms_per_row = 3;
my $head=<<"END";
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html>
<head>
<style type="text/css">
div.float {
	float: left;
	margin-right: 30px;
	margin-left: 30px;
}
div.float p {
	text-align: center;
}
div.spacer {
	clear: both;
}
</style>
</head>
<body>
<a name="top">
END

my $foot=<<"END";
</body>
</html>
END

my $spacer="\n".'<div class="spacer">&nbsp;</div>'."\n";


#Take all files in responsibility

 

open(OUT, ">".$outdir."/logos.html");
print OUT $head;
    
    
print OUT $spacer.'<HR size="3" color="black">'."\n";

#Get the KLD values
my $KLD_pres=0;
my @KLD=();
if ( -e  "$outdir/KLD/KLD.txt") {
    open IN, "$outdir/KLD/KLD.txt", or die;
    while (my $l=<IN>) {
	#chomp($l);
	$l =~ s/\r?\n$//;
	my @a=split(' ', $l);
	$KLD[$a[0]]=$a[1];
    }
    close IN;
    $KLD_pres=1;
}
    

#Get the single logo

my @single_list=();
if ($logo_type eq "LoLa") {
    @single_list = <$outdir/$logo_dir/LoLa_1-*.png>;
} elsif ($logo_type eq "Seq2Logo") {
    @single_list = <$outdir/$logo_dir/Seq2Logo_1-*.png>;
}

my $single;
my $size;
#Make this is really the single PWM (there should be only one, but issues can arise if some files are called PWM_1-*)
for (my $i = 0; $i < scalar(@single_list); $i++) {
    my @a = split("-", $single_list[$i]);
    my @b = split('\.png', $a[(scalar @a)-1]);
    if ($b[0]>=0 && $b[0]<=10000000) {
	my @c=split('\/', $single_list[$i]);
	$single=$c[(scalar @c)-1];
	$size=$b[0];
    }
}

print OUT '<p><font size="7">'.$name.'</font></p>'."\n";
print OUT '<font size="5"><b>Logo of the single cluster</b></font>'."\n";
if ($KLD_pres==1) {
    print OUT '<p><font size="3"><b>Score = '.$KLD[1].'</b></font></p>'."\n";
}
print OUT '<div class="spacer">&nbsp;</div>'."\n";
print OUT '<div class="float"><img src="'.$logo_dir.'/'.$single.'" width="250" height="150"><p>'.$size.' (1.000)</p></div>'."\n\n";

#Get the multiple logos
    
print OUT '<div class="spacer">&nbsp;</div>'."\n";
print OUT '<font size="5"><b>Logos of multiple clusters</b></font>'."\n";
#print OUT '<div class="spacer">&nbsp;</div>'."\n";

    
my @multiple_list = ();
for (my $sm=2; $sm<=20; $sm++) {
    my @multiple=();
    my @multiple_size=();
    $size=0;

    #Read the first line of the responsibility file
    my @name=();
    if ( -e "$outdir/responsibility/resp_$sm.txt") {
	open IN, "$outdir/responsibility/resp_$sm.txt", or die;
	my $l=<IN>;
	my @a=split(' ', $l);
	for (my $i=1; $i <= $sm; $i++) {
	    push @name, $a[$i];
	}
	close IN;
    } else {
	for (my $i=1; $i<=$sm; $i++) {
	    push @name, $i;
	}
    }
	
    for (my $tsm=1; $tsm<=$sm; $tsm++) {

	if ($logo_type eq "LoLa") {
	    @multiple_list = <$outdir/$logo_dir/LoLa_$sm\_$tsm-*.png>;
	} elsif ($logo_type eq "Seq2Logo"){
	    @multiple_list = <$outdir/$logo_dir/Seq2Logo_$sm\_$tsm-*.png>;
	}
	for (my $i = 0; $i < scalar(@multiple_list); $i++) {
	    my @a = split("-", $multiple_list[$i]);
	    my @b = split('\.png', $a[(scalar @a)-1]);
	    if ($b[0]>=0 && $b[0]<=10000000) {
		my @c=split('\/', $multiple_list[$i]);
		push @multiple, $c[(scalar @c)-1];
		push @multiple_size, $b[0];
		$size=$size+$b[0];
	    }
	}
    }
    if (scalar @multiple == $sm) {
	#print OUT '<div class="spacer">&nbsp;</div>'."\n";
	print OUT '<p><font size="4"><b>'.$sm.' clusters</b></font></p>'."\n";
	if ($KLD_pres==1) {
	    print OUT '<p><font size="3"><b>Score = '.$KLD[$sm].'</b><font></p>'."\n";
	}
	for (my $i = 0; $i < scalar(@multiple); $i++) {
	    my $p= sprintf("%.3f", $multiple_size[$i]/$size);
	    if ($multiple_size[$i]>0) {
		printf OUT '<div class="float"><img src="'.$logo_dir.'/'.$multiple[$i].'" width="250" height="150"><p>'.$name[$i].'  -  '.$multiple_size[$i].' ('.$p.')</p></div>'."\n";
	    } else {
		printf OUT '<div class="float"><img src="" width="250" height="150"><p>'.$name[$i].'  -  '.$multiple_size[$i].' ('.$p.')</p></div>'."\n";
	    }
	}
	print OUT '<div class="spacer">&nbsp;</div>'."\n";
    }
}

print OUT '<div class="spacer">&nbsp;</div>'."\n";
print OUT '</body>'."\n";
print OUT '</html>'."\n";
close OUT;



