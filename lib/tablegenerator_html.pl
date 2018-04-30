#!usr/bin/perl
use List::MoreUtils qw(uniq);

use strict;


my $outdir = $ARGV[0];
my $logo_type=$ARGV[1];
my $name=$ARGV[2];
my $naa_min=$ARGV[3];
my $naa_max=$ARGV[4];
my $ncl=$ARGV[5];
my $trash=$ARGV[6];
my $lcore=$ARGV[7];


my $logo_dir="logos";

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
$KLD_pres=0;

my $best_ncl=0;
if ( -e  "$outdir/KLD/best_ncl.txt") {
    open IN, "$outdir/KLD/best_ncl.txt", or die;
    my $l=<IN>;
    #chomp($l);
    $l =~ s/\r?\n$//;
    my @a=split(' ', $l);
    $best_ncl=$a[4];
    close IN;
}

my @full_size=([]);

my @list=();

for (my $s=$naa_min; $s<=$naa_max; $s++) {
    push @list, "$s";
}


my $logo_path;
my $s;
my $p;
foreach $s (@list) {
    
    my $size;
       
    open(OUT, ">".$outdir."/logos_html/logos_L$s.html");
    $logo_path="../";
    
    print OUT $head;
    print OUT $spacer.'<HR size="3" color="black">'."\n";
    print OUT '<p><font size="7">'.$name.' -- '.$s.'-mers</font></p>'."\n";
    print OUT '<p>'."\n";
    for (my $s1=$naa_min; $s1<=$naa_max; $s1++) {
	print OUT '<a href="./logos_L'.$s1.'.html">'.$s1.'-mers</a>'." | \n";
    }
    print OUT '</p>'."\n";
    print OUT '<p><a href="../length_distribution.html">Length distributions</a></p>'."\n";
    print OUT '<div class="spacer">&nbsp;</div>'."\n";

    if($s==$lcore){
	open(OUT2, ">".$outdir."/logos.html");
	
	print OUT2 $head;
	print OUT2 $spacer.'<HR size="3" color="black">'."\n";
	print OUT2 '<p><font size="7">'.$name.' -- core ('.$s.'-mers)</font></p>'."\n";
	print OUT2 '<p>'."\n";
	for (my $s1=$naa_min; $s1<=$naa_max; $s1++) {
	    print OUT2 '<a href="./logos_html/logos_L'.$s1.'.html">'.$s1.'-mers</a>'." | \n";
	}
	print OUT2 '</p>'."\n";
	print OUT2 '<p><a href="./length_distribution.html">Length distributions</a></p>'."\n";
	print OUT2 '<div class="spacer">&nbsp;</div>'."\n";
	
    }
    
    my @multiple_list = ();
    my @multiple_trash = ();
    for (my $sm=1; $sm<=$ncl+$trash; $sm++) {
	my @multiple=();
	my @multiple_size=();
	my @trash=();
	my @trash_size=();
	$size=0;
	
	
	
	for (my $tsm=1; $tsm<=$sm; $tsm++) {

	    if ($logo_type eq "LoLa") {
		@multiple_list = <$outdir/$logo_dir/LoLa_L$s\_$sm\_$tsm-*.png>;
	    } elsif ($logo_type eq "Seq2Logo") {
		@multiple_list = <$outdir/$logo_dir/Seq2Logo_L$s\_$sm\_$tsm-*.png>;
	    }
	    for (my $i = 0; $i < scalar(@multiple_list); $i++) {
		my @a = split("-", $multiple_list[$i]);
		my @b = split('\.png', $a[(scalar @a)-1]);
		if ($b[0]>=0 && $b[0]<=10000000) {
		    my @c=split('\/', $multiple_list[$i]);
		    push @multiple, $c[(scalar @c)-1];
		    push @multiple_size, $b[0];
		    $full_size[$sm][$tsm]=$full_size[$sm][$tsm]+$b[0];
		    $size=$size+$b[0];
		}
	    }
	}
	if ($logo_type eq "LoLa") {
	    @multiple_trash = <$outdir/$logo_dir/LoLa_L$s\_$sm\_Trash-*.png>;
	} elsif ($logo_type eq "Seq2Logo") {
	    @multiple_trash = <$outdir/$logo_dir/Seq2Logo_L$s\_$sm\_Trash-*.png>;
	}
	for (my $i = 0; $i < scalar(@multiple_trash); $i++) {
	    my @a = split("-", $multiple_trash[$i]);
	    my @b = split('\.png', $a[(scalar @a)-1]);
	    if ($b[0]>=0 && $b[0]<=10000000) {
		my @c=split('\/', $multiple_trash[$i]);
		push @trash, $c[(scalar @c)-1];
		push @trash_size, $b[0];
		$size=$size+$b[0];
	    }
	}
	
	if (scalar @multiple == $sm) {
	    if($sm==1){
		print OUT '<p><font size="4"><b>'.$sm.' motif '."\n";
	    }else{
		print OUT '<p><font size="4"><b>'.$sm.' motifs'."\n";
	    }
	    print OUT '</b></font></p>'."\n";
	    
	    if ($KLD_pres==1) {
		print OUT '<p><font size="3"><b>Score = '.$KLD[$sm].'</b><font></p>'."\n";
	    }
	    
	    for (my $i = 0; $i < scalar(@multiple); $i++) {
		if($size>0){
		    $p=sprintf("%.3f", $multiple_size[$i]/$size);
		} else {
		    $p=0;
		}
		my $ti=$i+1;
		if ($multiple_size[$i]>0) {
		    printf OUT '<div class="float"><img src="'.$logo_path.$logo_dir.'/'.$multiple[$i].'" width="250" height="150"><p>'.$ti.'  -  '.$multiple_size[$i].' ('.$p.')</p></div>'."\n";
		} else {
		    printf OUT '<div class="float"><img src="" width="250" height="150"><p>'.$ti.'  -  '.$multiple_size[$i].' ('.$p.')</p></div>'."\n";
		    #printf OUT '<div class="float">EMPTY<p>'.$ti.'  -  '.$multiple_size[$i].' ('.$p.')</p></div>'."\n";
		}
	    }
	    if (scalar @trash>0) {
		if($size>0){
		    $p= sprintf("%.3f", $trash_size[0]/$size);
		} else {
		    $p=0;
		}
		if ($trash_size[0]>0) {
		    printf OUT '<div class="float"><img src="'.$logo_path.$logo_dir.'/'.$trash[0].'" width="250" height="150"><p>Trash  -  '.$trash_size[0].' ('.$p.')</p></div>'."\n";
		} else {
		    printf OUT '<div class="float"><img src="" width="250" height="150"><p>Trash  -  '.$trash_size[0].' ('.$p.')</p></div>'."\n";
		    #printf OUT '<div class="float">EMPTY<p>Trash  -  '.$trash_size[0].' ('.$p.')</p></div>'."\n";
		}
		    
	    }
	    print OUT '<div class="spacer">&nbsp;</div>'."\n";

	    if($s==$lcore){

		if($sm==1){
		    print OUT2 '<p><font size="4"><b>'.$sm.' motif '."\n";
		}else{
		    print OUT2 '<p><font size="4"><b>'.$sm.' motifs'."\n";
		}
		print OUT2 '</b></font></p>'."\n";
		
		if ($KLD_pres==1) {
		    print OUT2 '<p><font size="3"><b>Score = '.$KLD[$sm].'</b><font></p>'."\n";
		}
		
		for (my $i = 0; $i < scalar(@multiple); $i++) {
		    if($size>0){
			$p= sprintf("%.3f", $multiple_size[$i]/$size);
		    }else {
			$p=0;
		    }
		    my $ti=$i+1;
		    if ($multiple_size[$i]>0) {
		    printf OUT2 '<div class="float"><img src="'.$logo_dir.'/'.$multiple[$i].'" width="250" height="150"><p>'.$ti.'  -  '.$multiple_size[$i].' ('.$p.')</p></div>'."\n";
		} else {
		    printf OUT2 '<div class="float"><img src="" width="250" height="150"><p>'.$ti.'  -  '.$multiple_size[$i].' ('.$p.')</p></div>'."\n";
		}
		}
		if (scalar @trash>0) {
		    if($size>0){
			$p= sprintf("%.3f", $trash_size[0]/$size);
		    }else {
			$p=0;
		    }
		    if($trash_size[0]>0){
			printf OUT2 '<div class="float"><img src="'.$logo_dir.'/'.$trash[0].'" width="250" height="150"><p>Trash  -  '.$trash_size[0].' ('.$p.')</p></div>'."\n";
		    } else {
			printf OUT2 '<div class="float"><img src="" width="250" height="150"><p>Trash  -  '.$trash_size[0].' ('.$p.')</p></div>'."\n";
		    }
		}
		print OUT2 '<div class="spacer">&nbsp;</div>'."\n";
		
	    }
	    
	}
    }
    
    print OUT '<div class="spacer">&nbsp;</div>'."\n";
    print OUT '</body>'."\n";
    print OUT '</html>'."\n";
    close OUT;

    if($s==$lcore){
	print OUT2 '<div class="spacer">&nbsp;</div>'."\n";
	print OUT2 '</body>'."\n";
	print OUT2 '</html>'."\n";
	close OUT;
    }
   
    
}


open OUT, ">".$outdir."/length_distribution.html";
print OUT $head;

print OUT $spacer.'<HR size="3" color="black">'."\n";

print OUT '<p><font size="7">'.$name.' -- Peptide length distributions</font></p>'."\n";

print OUT '<p>'."\n";

print OUT '<a href="./logos_html/logos_L'.$lcore.'.html"> logos ('.$lcore.'-mers core)</a>'."\n";

for(my $i=1; $i<=$ncl; $i++){

    if($i==1){
	print OUT '<p><font size="4"><b>'.$i.' motif '."\n";
    }else{
	print OUT '<p><font size="4"><b>'.$i.' motifs'."\n";
    }
    print OUT '</b></font></p>'."\n";
    
    for(my $j=1; $j<=$i; $j++){
	printf OUT '<div class="float"><img src="./weights/plots/lg_'.$i.'_'.$j.'.png" width="200" height="200"></div>'."\n";
    }
    if($trash==1){
	printf OUT '<div class="float"><img src="./weights/plots/lg_'.$i.'_Trash.png" width="200" height="200"></div>'."\n";
    }
    print OUT '<div class="spacer">&nbsp;</div>'."\n";
    
}

print OUT '<div class="spacer">&nbsp;</div>'."\n";
print OUT '</body>'."\n";
print OUT '</html>'."\n";
close OUT;

close OUT;
