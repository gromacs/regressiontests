#!/usr/bin/perl -w

$parallel = 0;
$double   = 0;
$verbose  = 0;
$etol     = 0.05;
$ttol     = 0.001;
$ref      = "";

# virial - this tests the shifted force contribution.
# However, it is a sum of very many large terms, so it is
# often numerically imprecise.
$virtol_rel   = 0.01;
$virtol_abs   = 0.1;

# force tolerance is measured by calculating scalar products.
# tolerance 0.001 means the scalar product between the two
# compared forces should be at least 0.999.
$ftol_rel     = 0.001;
$ftol_sprod   = 0.001;

sub setup_vars {
    if ( $parallel > 0 ) {
	$mdprefix = "mpirun -c $parallel ";
    }
    else {
	$mdprefix = "";
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
    $tmp = "checkforce.tmp";
    $cfor = "checkforce.out";
    system("gmxcheck -f $reftrr -f2 traj -tol $ftol_rel > $cfor 2> /dev/null");    
    `grep "f\\[" $cfor > $tmp`;
    $nerr_force = 0;
    
    open(FIN,"$tmp");
    while($line=<FIN>)
    {
	@f1=split(" ",substr($line,10,38));
	@f2=split(" ",substr($line,53,38));
	
	$l1 = sqrt($f1[0]*$f1[0]+$f1[1]*$f1[1]+$f1[2]*$f1[2]);
	$l2 = sqrt($f2[0]*$f2[0]+$f2[1]*$f2[1]+$f2[2]*$f2[2]);
	$sprod = ($f1[0]*$f2[0]+$f1[1]*$f2[1]+$f1[2]*$f2[2])/($l1*$l2);
	
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
    $tmp = "checkvir.tmp";
    $cvir = "checkvir.out";
    
    system("gmxcheck -e $refedr -e2 ener -tol $virtol_rel -lastener Vir-ZZ > $cvir 2> /dev/null");   
    
    `grep "Vir-" $cvir > $tmp`;
    $nerr_vir = 0;
    
    open(VIN,"$tmp");
    while($line=<VIN>)
    {
	@v1=split(" ",substr($line,26,14));
	@v2=split(" ",substr($line,52,13));
	
	$diff = abs($v1[0]-$v2[0]);
	
	$norm = abs($v1[0])+abs($v2[0]);
	
	if((2*$diff > $virtol_rel *$norm) && ($diff>$virtol_abs))
	{
	    $nerr_vir = $nerr_vir + 1;
	}
    }     
    close(VIN);
    unlink($tmp,$cvir);
    
    return $nerr_vir;
}

sub test_systems {
    setup_vars();
    $npassed = 0;
    foreach $dir ( @_ ) {
	if ( -d $dir ) {
	    chdir($dir);
	    if ($verbose > 0) {
		print "Testing $dir . . . ";
	    }
	    
	    $nerror = 0;
	    $ndx = "";
	    if ( -f "index.ndx" ) {
		$ndx = "-n index";
	    }
	    $par = "";
	    if ($parallel > 1) {
		$par = "-np $parallel";
	    }
	    system("$grompp -maxwarn 10 $ndx $par > grompp.out 2>&1");
	    
	    if (! -f "topol.tpr") {
		print ("No topol.tpr file in $dir. grompp failed\n");	    
		$nerror = 1;
	    }
	    if ($nerror == 0) {
		$reftpr = "${ref}.tpr";
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
		    $refedr = "${ref}.edr";
		    if (! -f  $refedr) {
			print ("No $refedr file in $dir.\n");
			print ("This means you are not really testing $dir\n");
			system("cp ener.edr $refedr");
		    }
		    $reftrr = "${ref}.trr";
		    if (! -f $reftrr ) {
			print ("No $reftrr file in $dir.\n");
			print ("This means you are not really testing $dir\n");
			system("cp traj.trr $reftrr");
		    }
		    # Now do the real tests
		    system("gmxcheck -e $refedr -e2 ener -tol $etol -lastener Potential > checkpot.out 2> /dev/null");
		    
		    $nerr_pot   = `grep step checkpot.out | grep -v gcq | wc -l`;
		    $nerr_force = check_force();
		    $nerr_vir   = check_virial();
		    
		    $nerror = ($nerr_pot || $nerr_force || $nerr_vir);
		}
		else {
		    $nerror = 1;
		}
	    }
	    if ($nerror > 0) {
		print "FAILED. Check files in $dir\n";
	    }
	    else {
		@args = glob("#*# *.out topol.tpr confout.gro ener.edr md.log traj.trr");
		unlink(@args);
		
		if ($verbose > 0) {
		    $nmdp = `diff grompp.mdp mdout.mdp | grep -v host | grep -v date | grep -v user | grep -v generated | wc -l`;
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
    @kkk = `/bin/ls | grep -v CVS`;
    for $k ( @kkk ) {
	chomp $k;
    }
    return @kkk;
}

sub cleandirs {
    $mydir = shift;
    chdir($mydir);
    foreach $dir ( my_glob() ) {
	if ( -d $dir ) {
	    chdir($dir);
	    print "Cleaning $dir\n"; 
	    @args = glob("#*# *~ *.out core.* *.xvg topol.tpr confout.gro ener.edr md.log traj.trr" );
	    unlink (@args);
	    chdir("..");
	}
    }
    chdir("..");
}

sub refcleandir {
    $sdir = shift;
    
    if (-d $sdir ) {
	cleandirs($sdir);
	chdir($sdir);
	@mydirs = my_glob();
	foreach $dir ( @mydirs ) {
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
    $dirs = shift;
    chdir($dirs);
    @kernels = my_glob();
    $nn = $#kernels + 1;
    $npassed = test_systems(@kernels);
    if ($npassed < $nn) {
	printf("%d out of $nn $dirs tests FAILED\n",$nn-$npassed);
    }
    else {
	printf("All $nn $dirs tests PASSED\n");
    }
    chdir("..");
}

sub test_pdb2gmx {
    $logfn = "pdb2gmx.log";

    setup_vars();    
    chdir("pdb2gmx");
    open (LOG,">$logfn") || die("Opening $logfn for writing");
    $npdb_dir = 0;
    @pdb_dirs = ();
    foreach $pdb ( glob("*.pdb") ) {
	$pdir = "pdb-$pdb";
	@kkk  = split('\.',$pdir);
	$dir  = $kkk[0];
	$pdb_dirs[$npdb_dir++] = $dir;
	mkdir($dir);
	chdir($dir);
	foreach $ff ( "G43a1", "oplsaa", "gmx" ) {
	    mkdir("ff$ff");
	    chdir("ff$ff");
	    if ( $ff eq "oplsaa"  ) {
		@www = ( "tip3p", "tip4p", "tip5p" );
	    }
	    else {
		@www = ( "spc", "spce", "tip4p" );
	    }
	    foreach $dd ( "none", "h", "aromatics" ) {
		mkdir("$dd");
		chdir("$dd");
		foreach $ww ( @www ) {
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
    
    $nsuccess = `wc -l ener.log | awk '{print \$1}'`;
    if ( $nsuccess != 54 ) {
	print "Error not all 54 tests have been done successfully\n";
	print "Only $nsuccess energies in the log file\n";
    }
    else {
	open(EEE,"paste reference.log ener.log |");
	$n = 0;
	$nerror = 0;
	$header = 0;
	while ($line = <EEE>) {
	    chomp($line);
	    @tmp = split(' ',$line);
	    if ((($tmp[3]-$tmp[7])/($tmp[3]+$tmp[7])) > $etol) {
		$nerror++;
		if (!$header) {
		    $header = 1;
		    printf("Sim.   Reference   This test\n");
		}
		printf("%4d  %10g  %10g\n",$n,$tmp[3],$tmp[7]);
	    }
	    $n++;
	}
	close EEE;
	
	if ( $nerror != 0 ) {
	    print "There were $nerror differences in final energy with the reference file\n";
	}
	else {
	    print "pdb2gmx test PASSED\n";
	    system("rm -rf @pdb_dirs");
	    unlink("ener.log","pdb2gmx.log");
	} 
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

for ($kk=0; ($kk <= $#ARGV); $kk++) {
    $arg = $ARGV[$kk];
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
