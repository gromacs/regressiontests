#!/usr/bin/perl -w

use strict;

my $parallel = 0;
my $double   = 0;
my $verbose  = 0;
my $etol     = 0.05;
my $ttol     = 0.001;
my $ref      = "";

# virial - this tests the shifted force contribution.
# However, it is a sum of very many large terms, so it is
# often numerically imprecise.
my $virtol_rel   = 0.01;
my $virtol_abs   = 0.1;

# force tolerance is measured by calculating scalar products.
# tolerance 0.001 means the scalar product between the two
# compared forces should be at least 0.999.
my $ftol_rel     = 0.001;
my $ftol_sprod   = 0.001;

# globals used for programs
my $grompp = "";
my $mdrun  = "";

sub setup_vars {
    my $mdprefix = "";
    if ( $parallel > 0 ) {
	$mdprefix = "mpirun -c $parallel ";
    }
    
    if ( $double > 0 ) {
	$grompp = "grompp_d";
	$mdrun  = $mdprefix . "mdrun_d";
	$ref    = "reference_d";
    }
    else {
	$grompp = "grompp";
	$mdrun  = $mdprefix . "mdrun";
	$ref    = "reference_s";
    }
}

sub check_force()
{
    my $tmp = "checkforce.tmp";
    my $cfor = "checkforce.out";
    my $reftrr = "${ref}.trr";
    system("gmxcheck -f $reftrr -f2 traj -tol $ftol_rel > $cfor 2> /dev/null");    
    `grep "f\\[" $cfor > $tmp`;
    my $nerr_force = 0;
    
    open(FIN,"$tmp");
    while(my $line=<FIN>)
    {
	my @f1=split(" ",substr($line,10,38));
	my @f2=split(" ",substr($line,53,38));
	
	my $l1 = sqrt($f1[0]*$f1[0]+$f1[1]*$f1[1]+$f1[2]*$f1[2]);
	my $l2 = sqrt($f2[0]*$f2[0]+$f2[1]*$f2[1]+$f2[2]*$f2[2]);
	my $sprod = ($f1[0]*$f2[0]+$f1[1]*$f2[1]+$f1[2]*$f2[2])/($l1*$l2);
	
	if( $sprod < (1.0-$ftol_sprod))
	{
	    $nerr_force = $nerr_force + 1;
	}
    }     
    close(FIN);
    unlink($tmp,$cfor);
    
    return $nerr_force;
}

sub check_virial()
{
    my $tmp = "checkvir.tmp";
    my $cvir = "checkvir.out";
    my $refedr = "${ref}.edr";
    system("gmxcheck -e $refedr -e2 ener -tol $virtol_rel -lastener Vir-ZZ > $cvir 2> /dev/null");   
    
    `grep "Vir-" $cvir > $tmp`;
    my $nerr_vir = 0;
    
    open(VIN,"$tmp");
    while(my $line=<VIN>)
    {
	my @v1=split(" ",substr($line,26,14));
	my @v2=split(" ",substr($line,52,13));
	
	my $diff = abs($v1[0]-$v2[0]);
	
	my $norm = abs($v1[0])+abs($v2[0]);
	
	if((2*$diff > $virtol_rel *$norm) && ($diff>$virtol_abs))
	{
	    $nerr_vir = $nerr_vir + 1;
	}
    }     
    close(VIN);
    unlink($tmp,$cvir);
    
    return $nerr_vir;
}

sub check_xvg {
    my $refx = shift;
    my $kkk  = shift;
    my $ndx1 = shift;
    my $ndx2 = shift;
    
    my $nerr = 0;
    if ((-f $refx) && (-f $kkk)) {
	open(EEE,"paste $refx $kkk |");
	my $n = 0;
	my $header = 0;
	while (my $line = <EEE>) {
	    if ((index($line,"#") < 0) && (index($line,"\@") < 0)) {
		chomp($line);
		my @tmp = split(' ',$line);
		my $x1 = $tmp[$ndx1];
		my $x2 = $tmp[$ndx2];
		if ((($x1-$x2)/($x1+$x2)) > $etol) {
		    $nerr++;
		    if (!$header) {
			$header = 1;
			printf("N      Reference   This test\n");
		    }
		    printf("%4d  %10g  %10g\n",$n,$tmp[3],$tmp[7]);
		}
		$n++;
	    }
	}
	close EEE;
    }
    return $nerr;
}

sub test_systems {
    setup_vars();
    my $npassed = 0;
    foreach my $dir ( @_ ) {
	if ( -d $dir ) {
	    chdir($dir);
	    if ($verbose > 0) {
		print "Testing $dir . . . ";
	    }
	    
	    my $nerror = 0;
	    my $ndx = "";
	    if ( -f "index.ndx" ) {
		$ndx = "-n index";
	    }
	    my $par = "";
	    if ($parallel > 1) {
		$par = "-np $parallel";
	    }
	    system("$grompp -maxwarn 10 $ndx $par > grompp.out 2>&1");
	    
	    if (! -f "topol.tpr") {
		print ("No topol.tpr file in $dir. grompp failed\n");	    
		$nerror = 1;
	    }
	    if ($nerror == 0) {
		my $reftpr = "${ref}.tpr";
		if (! -f $reftpr) {
		    print ("No $reftpr file in $dir\n");
		    print ("This means you are not really testing $dir\n");
		    system("cp topol.tpr $reftpr");
		}
		system("gmxcheck -s1 $reftpr -s2 topol.tpr -tol $ttol > checktpr.out 2>&1");
		$nerror = `grep step checktpr.out | grep -v gcq | wc -l`;
		if ($nerror > 0) {
		    print "topol.tpr file different from $reftpr. Check files in $dir\n";
		}
	    }
	    if ($nerror == 0) {
		# Do the mdrun at last!
		system("$mdrun > mdrun.out 2>&1");
		
		# First check whether we have any output
		if ((-f "ener.edr" ) && (-f "traj.trr")) {
		    # Now check whether we have any reference files
		    my $refedr = "${ref}.edr";
		    if (! -f  $refedr) {
			print ("No $refedr file in $dir.\n");
			print ("This means you are not really testing $dir\n");
			system("cp ener.edr $refedr");
		    }
		    my $reftrr = "${ref}.trr";
		    if (! -f $reftrr ) {
			print ("No $reftrr file in $dir.\n");
			print ("This means you are not really testing $dir\n");
			system("cp traj.trr $reftrr");
		    }
		    # Now do the real tests
		    system("gmxcheck -e $refedr -e2 ener -tol $etol -lastener Potential > checkpot.out 2> /dev/null");
		    
		    my $nerr_pot   = `grep step checkpot.out | grep -v gcq | wc -l`;
		    my $nerr_force = check_force();
		    my $nerr_vir   = check_virial();
		
		    my $nerr_xvg   = check_xvg("${ref}.xvg","dgdl.xvg",1,3);
		    
		    $nerror = ($nerr_pot || $nerr_force || 
			       $nerr_vir || $nerr_xvg);
		}
		else {
		    $nerror = 1;
		}
	    }
	    if ($nerror > 0) {
		print "FAILED. Check files in $dir\n";
	    }
	    else {
		my @args = glob("#*# *.out topol.tpr confout.gro ener.edr md.log traj.trr");
		unlink(@args);
		
		if ($verbose > 0) {
		    my $nmdp = `diff grompp.mdp mdout.mdp | grep -v host | grep -v date | grep -v user | grep -v generated | wc -l`;
		    if ( $nmdp > 2) {
			printf("PASSED but check mdp file differences\n");
		    }
		    else {
			print "PASSED\n";
			system("rm -f mdout.mdp");
		    }
		}
		$npassed++;
	    }
	    chdir("..");
	}
    }
    return $npassed;
}

sub my_glob {
    my @kkk = `/bin/ls | grep -v CVS`;
    for my $k ( @kkk ) {
	chomp $k;
    }
    return @kkk;
}

sub cleandirs {
    my $mydir = shift;
    chdir($mydir);
    foreach my $dir ( my_glob() ) {
	if ( -d $dir ) {
	    chdir($dir);
	    print "Cleaning $dir\n"; 
	    my @args = glob("#*# *~ *.out core.* *.xvg topol.tpr confout.gro ener.edr md.log traj.trr" );
	    unlink (@args);
	    chdir("..");
	}
    }
    chdir("..");
}

sub refcleandir {
    my $sdir = shift;
    
    if (-d $sdir ) {
	cleandirs($sdir);
	chdir($sdir);
	my @mydirs = my_glob();
	foreach my $dir ( @mydirs ) {
	    if ( -d $dir ) {
		chdir($dir);
		print "Removing reference files in $dir\n"; 
		unlink ("${ref}.edr","${ref}.tpr","${ref}.trr");
		chdir("..");
	    }
	}
	chdir("..");
    }
}

sub test_dirs {
    my $dirs = shift;
    chdir($dirs);
    my @kernels = my_glob();
    my $nn = $#kernels + 1;
    my $npassed = test_systems(@kernels);
    if ($npassed < $nn) {
	printf("%d out of $nn $dirs tests FAILED\n",$nn-$npassed);
    }
    else {
	printf("All $nn $dirs tests PASSED\n");
    }
    chdir("..");
}

sub test_pdb2gmx {
    my $logfn = "pdb2gmx.log";

    setup_vars();    
    chdir("pdb2gmx");
    open (LOG,">$logfn") || die("Opening $logfn for writing");
    my $npdb_dir = 0;
    my @pdb_dirs = ();
    my $ntest    = 0;
    my $nerror   = 0;
    foreach my $pdb ( glob("*.pdb") ) {
	my $pdir = "pdb-$pdb";
	my @kkk  = split('\.',$pdir);
	my $dir  = $kkk[0];
	$pdb_dirs[$npdb_dir++] = $dir;
	mkdir($dir);
	chdir($dir);
	foreach my $ff ( "G43a1", "oplsaa", "gmx" ) {
	    mkdir("ff$ff");
	    chdir("ff$ff");
	    my @www = ();
	    if ( $ff eq "oplsaa"  ) {
		@www = ( "tip3p", "tip4p", "tip5p" );
	    }
	    else {
		@www = ( "spc", "spce", "tip4p" );
	    }
	    foreach my $dd ( "none", "h", "aromatics" ) {
		mkdir("$dd");
		chdir("$dd");
		foreach my $ww ( @www ) {
		    $ntest++;
		    my $line = "";
		    printf(LOG "****************************************************\n");
		    printf(LOG "** PDB = $pdb FF = $ff VSITE = $dd WATER = $ww\n");
		    printf(LOG "****************************************************\n");
		    mkdir("$ww");
		    chdir("$ww");
		    printf(LOG "****************************************************\n");
		    printf(LOG "**  Running pdb2gmx\n");
		    printf(LOG "****************************************************\n");
		    open(PIPE,"pdb2gmx -f ../../../../$pdb -ff $ff -ignh -vsite $dd -water $ww 2>&1 |");
		    while ($line = <PIPE>) { printf(LOG $line); } close PIPE;
		    printf(LOG "****************************************************\n");
		    printf(LOG "**  Running editconf\n");
		    printf(LOG "****************************************************\n");
		    open(PIPE,"editconf -o b4em -box 5 5 5 -c -f conf 2>&1 |");
		    while ($line = <PIPE>) { printf(LOG $line); } close PIPE;
		    printf(LOG "****************************************************\n");
		    printf(LOG "**  Running grompp\n");
		    printf(LOG "****************************************************\n");
		    open(PIPE,"$grompp -f ../../../../em -c b4em 2>&1 |");
		    while ($line = <PIPE>) { printf(LOG $line); } close PIPE;
		    printf(LOG "****************************************************\n");
		    printf(LOG "**  Running mdrun\n");
		    printf(LOG "****************************************************\n");
		    open(PIPE,"$mdrun 2>&1 |");
		    while ($line = <PIPE>) { printf(LOG $line); } close PIPE;
		    chdir("..");
		}
		chdir("..");
	    }
	    chdir("..");
	}
	chdir("..");
    }
    close LOG;
    
    system("grep 'Potential Energy' pdb2gmx.log > ener.log");
    
    my $nsuccess = `wc -l ener.log | awk '{print \$1}'`;
    chomp($nsuccess);
    if ( $nsuccess != $ntest ) {
	print "Error not all $ntest pdb2gmx tests have been done successfully\n";
	print "Only $nsuccess energies in the log file\n";
    }
    else {
	$nerror = check_xvg("reference.log","ener.log",3,7);
	
	if ( $nerror != 0 ) {
	    print "There were $nerror differences in final energy with the reference file\n";
	}
	else {
	    print "All $ntest pdb2gmx tests PASSED\n";
	    system("rm -rf @pdb_dirs");
	    unlink("ener.log","pdb2gmx.log");
	} 
    }
    if ($nerror > 0) {
	print "pdb2gmx tests FAILED\n";
    }
    chdir("..");
}

sub clean_all {
    cleandirs("simple");
    cleandirs("complex");
    cleandirs("kernel");
    chdir("pdb2gmx");
    system("rm -rf pdb-* ener.log pdb2gmx.log");
    chdir("..");
}

sub usage {
    print "Usage: ./gmxtest.pl [ -np N ] [-verbose ] [ -double ] [ simple | complex | kernel | pdb2gmx | all ]\n";
    print "   or: ./gmxtest.pl clean | refclean | dist\n";
    exit "1";
}

my $kk = 0;
for ($kk=0; ($kk <= $#ARGV); $kk++) {
    my $arg = $ARGV[$kk];
    if ($arg eq "simple") {
	test_dirs("simple");
    }
    elsif ($arg eq "complex") {
	test_dirs("complex");
    }
    elsif ($arg eq "kernel" ) {
	test_dirs("kernel");
    }
    elsif ($arg eq "pdb2gmx" ) {
	test_pdb2gmx();
    }
    elsif ($arg eq "all" ) {
	test_dirs("simple");
	test_dirs("complex");
	test_dirs("kernel");
	test_pdb2gmx();
    }
    elsif ($arg eq "clean" ) {
        clean_all();
    }
    elsif ($arg eq "refclean" ) {
	setup_vars();
	refcleandir("simple");
	refcleandir("complex");
    }
    elsif ($arg eq "dist" ) {
	clean_all();
	chdir("..");
	system("tar --exclude CVS -czvf gmxtest.tgz gmxtest");
	chdir("gmxtest");
    }
    elsif ($arg eq "help" ) {
	usage();
    }
    elsif ($arg eq "-verbose") {
	$verbose = 1;
    }
    elsif ($arg eq "-double") {
	$double = 1;
    }
    elsif ($arg eq "-np" ) {
	if ($kk <$#ARGV) {
	    $kk++;
	    $parallel = $ARGV[$kk];
	    if ($parallel < 0) {
		$parallel = 0;
	    }
	    print "Will test on $parallel processors\n";
	}
    }
    else {
	usage();
    }
}

if ($kk == 0) {
    usage();
}
