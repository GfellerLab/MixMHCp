
#Create the Seq2Logo logos
#This is useful to create the logos with Seq2Logo.

#######################################################
#THE PATH TO SEQ2LOGO NEEDS TO BE IN YOUR PATH
#THE PATH TO 'gs' NEEDS TO BE HARD-CODED IN Seq2Logo.py
#######################################################

#Take all files in responsibility folder
$d=$ARGV[0]."/responsibility/";
$bkg_path=$ARGV[1];
$naa_min=$ARGV[2];
$naa_max=$ARGV[3];
$trash=$ARGV[4];

opendir($dir, $d);
@file_list=();
while (readdir $dir) {
    $f=$_;
    if ($f =~ /.txt/) {
	push @file_list, $f;
    }
}
closedir $dir;

$bkg="flat";

if ($bkg eq "flat") {
    $bg_freq="--bg $bkg_path/Seq2Logo_Flat_Bgdistr.txt";
} else {
    $bg_freq="";
}

$logo_type="Shannon";
if ($logo_type eq "KL") {
    $I=2;
    $y=2;
} elsif ($logo_type eq "Shannon") {
    $I=1;
    $y=0;
}


foreach $FILE (@file_list) {

    @seq=();
    @class=();
    
    @a=split('\.txt', $FILE);
    $f=$a[0];
    
    $p=rindex($f, "_");
    $Np=substr($f, $p+1, length($f)-$p);

    print "$Np\n";
    
    #Create the logo for all peptides based on a clustering of the peptides
    #There can be a problem since gaps are just removed in Seq2Logo (clearly not optimal), but it works well if no alignment is required.
    #However, the pseudocount can be added, so this is useful to compare from sequences directly analyzed with Seq2Logo

    for ($s1=$naa_min; $s1<=$naa_max; $s1++) {
	
	open IN, $d.$FILE, or die;
	$l=<IN>;
	
	@val=();
	@seq=();
	@class=();
	
	while ($l=<IN>) {
	    #chomp($l);
	    $l =~ s/\r?\n$//;
	    
	    @a=split(' ', $l);

	    if (length($a[0]) == $s1 ) {
	    
		$max=0;
		$pmax=0;
		for ($i=1; $i<=$Np+$trash; $i++) {
		    if ($a[$i]>$max) {
			$max=$a[$i];
			$pmax=$i;
		    }
		    $val[$i]=$val[$i]+$a[$i];
		}
		
		push @seq, $a[0];
		push @class, $pmax;
	    }
	    
	}
	close IN;



	if (! -d "$ARGV[0]/Seq2Logo/") {
	    system("mkdir $ARGV[0]/Seq2Logo");
	}
	for ($i=1; $i<=$Np+$trash; $i++) {
	   
	    if($trash==1 && $i==$Np+$trash){
		open OUT, ">$ARGV[0]/Seq2Logo/cluster_L$s1\_$Np\_Trash.fa";
	    } else {
		open OUT, ">$ARGV[0]/Seq2Logo/cluster_L$s1\_$Np\_$i.fa";
	    }
	    
	    
	    $ct=0;
	    for ($j=0; $j<scalar(@seq); $j++) {
		
		if ($class[$j]==$i) {
		    print OUT ">pep\n$seq[$j]\n";
		    $ct++;
		}
	    }
	    if ($ct==0) {
		print OUT ">pep\nEMPTY\n";
	    }
	    if ($ct>$val[$i]+2 || $ct < $val[$i]-2) {
		print "Many Ambiguous peptides...\n";
	    }
	    
	    close OUT;

	   
	    if($trash==1 && $i==$Np+$trash){
		system("Seq2Logo.py -f $ARGV[0]/Seq2Logo/cluster_L$s1\_$Np\_Trash.fa -o $ARGV[0]/Seq2Logo/cluster_L$s1\_$Np\_Trash -y $y:4.32 -b 0 --format PNG -I $I $bg_freq");
		system("cp $ARGV[0]/Seq2Logo/cluster_L$s1\_$Np\_Trash-001.png $ARGV[0]/Seq2Logo/Seq2Logo_L$s1\_$Np\_Trash-$ct.png");
	    } else {
		system("Seq2Logo.py -f $ARGV[0]/Seq2Logo/cluster_L$s1\_$Np\_$i.fa -o $ARGV[0]/Seq2Logo/cluster_L$s1\_$Np\_$i -y $y:4.32 -b 0 --format PNG -I $I $bg_freq");
		system("cp $ARGV[0]/Seq2Logo/cluster_L$s1\_$Np\_$i-001.png $ARGV[0]/Seq2Logo/Seq2Logo_L$s1\_$Np\_$i-$ct.png");
	    }
	    
	}
    }
	
}
system("mv $ARGV[0]/Seq2Logo/Seq2Logo_L* $ARGV[0]/logos/");
