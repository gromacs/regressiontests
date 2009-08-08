#!/usr/bin/perl -w

use strict;

my $parallel = 0;
my $double   = 0;
my $verbose  = 1;
my $etol     = 0.05;
my $ttol     = 0.001;
my $suffix   = '';
my $prefix   = '';
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

# global variables to explain some situations to the user
my $addversionnote = 0;

# trickery for program and reference file names
my $mdprefix = '';
my $ref      = '';
my %progs = ( 'grompp'   => 'grompp',
	      'mdrun'    => 'mdrun',
	      'pdb2gmx'  => 'pdb2gmx',
	      'gmxcheck' => 'gmxcheck',
	      'editconf' => 'editconf' );
sub setup_vars()
{
    # We assume that the name of executables match the 
    # pattern ${prefix}mdrun[_mpi][_d]${suffix} where
    # ${prefix} and ${suffix} are as defined above (or
    # over-ridden on the command line), "_d" indicates a
    # double-precision version, and (only in the case of
    # mdrun) "_mpi" indicates a parallel version compiled
    # with MPI.
    if ( $parallel > 0 ) {
	$progs{'mdrun'} .= "_mpi";
	$mdprefix = "mpirun -c $parallel";
    }
    foreach my $prog ( values %progs ) {
	$prog = $prefix . $prog;
	$prog .= "_d" if ( $double > 0 );
	$prog .= $suffix;
    }
    if ( $double > 0 ) {
	$ref    = "reference_d";
    }
    else {
	$ref    = "reference_s";
    }
}

sub check_force()
{
    my $tmp = "checkforce.tmp";
    my $cfor = "checkforce.out";
    my $reftrr = "${ref}.trr";
    system("$progs{'gmxcheck'} -f $reftrr -f2 traj -tol $ftol_rel > $cfor 2> /dev/null");    
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
    if ($nerr_force == 0) {
      unlink($tmp,$cfor);
    }
    return $nerr_force;
}

sub check_virial()
{
    my $tmp = "checkvir.tmp";
    my $cvir = "checkvir.out";
    my $refedr = "${ref}.edr";
    system("$progs{'gmxcheck'} -e $refedr -e2 ener -tol $virtol_rel -lastener Vir-ZZ > $cvir 2> /dev/null");   
    
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
    if ($nerr_vir == 0) {
      unlink($tmp,$cvir);
    }
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
			print("N      Reference   This test\n");
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
    my $npassed = 0;
    foreach my $dir ( @_ ) {
	if ( -d $dir ) {
	    chdir($dir);
	    if ($verbose > 0) {
		print "Testing $dir . . . ";
	    }
	    
	    my $nerror = 0;
	    my $ndx = "";
	    my $tprversionmismatch;
	    if ( -f "index.ndx" ) {
		$ndx = "-n index";
	    }
	    my $par = "";
	    if ($parallel > 1) {
		$par = "-np $parallel";
	    }
	    system("$progs{'grompp'} -maxwarn 10 $ndx > grompp.out 2>&1");
	    
	    my $error_detail = ' ';
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
		system("$progs{'gmxcheck'} -s1 $reftpr -s2 topol.tpr -tol $ttol > checktpr.out 2>&1");
		$nerror = `grep step checktpr.out | grep -v gcq | wc -l`;
		if ($nerror > 0) {
		    print "topol.tpr file different from $reftpr. Check files in $dir\n";
		}
		system("grep 'reading tpx file (reference_d.tpr) version .* with version .* program' checktpr.out >& /dev/null");
		$tprversionmismatch = (0 == ($? >> 8));
		if ($tprversionmismatch > 0) {
		    print "\nThe GROMACS version being tested is older than the reference version.\nPlease see the note at end of this output.\n";
		    $addversionnote = 1;
		}
	    }
	    if ($nerror == 0) {
		# Do the mdrun at last!
		system("$mdprefix $progs{'mdrun'} > mdrun.out 2>&1");
		
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
		    system("$progs{'gmxcheck'} -e $refedr -e2 ener -tol $etol -lastener Potential > checkpot.out 2> /dev/null");
		    
		    my $nerr_pot   = `grep step checkpot.out | grep -v gcq | wc -l`;
		    chomp($nerr_pot);
		    my $nerr_force = check_force();
		    my $nerr_vir   = check_virial();
		
		    my $nerr_xvg   = check_xvg("${ref}.xvg","dgdl.xvg",1,3);
		    
		    $nerror = ($nerr_pot || $nerr_force || 
			       $nerr_vir || $nerr_xvg);

		    my @error_detail;
		    push(@error_detail, "checkpot.out ($nerr_pot errors)") if ($nerr_pot > 0);
		    push(@error_detail, "checkvir.out ($nerr_vir errors)") if ($nerr_vir > 0);
		    push(@error_detail, "checkforce.out ($nerr_force errors)") if ($nerr_force > 0);
		    push(@error_detail, "${ref}.xvg ($nerr_xvg errors)") if ($nerr_xvg > 0);
		    $error_detail = join(', ', @error_detail) . ' ';
		}
		else {
		    $nerror = 1;
		}
	    }
	    if ($nerror > 0) {
		print "FAILED. Check ${error_detail}files in $dir\n";
	    }
	    else {
		my @args = glob("#*# *.out topol.tpr confout.gro ener.edr md.log traj.trr");
		#unlink(@args);
		
		if ($verbose > 0) {
		    my $nmdp = `diff grompp.mdp mdout.mdp | grep -v host | grep -v date | grep -v user | grep -v generated | wc -l`;
		    if ( $nmdp > 2 && `cat grompp.mdp | wc -l` > 50) {
			# if the input .mdp file is trivially short, then 
			# the above diff test will always fail
			print("PASSED but check mdp file differences\n");
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
	    my @args = glob("#*# *~ *.out core.* *.xvg topol.tpr confout.gro ener.edr md.log traj.trr *.tmp mdout.mdp step*.pdb *~ grompp*" );
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
    my @subdirs = my_glob();
    my $nn = $#subdirs + 1;
    my $npassed = test_systems(@subdirs);
    if ($npassed < $nn) {
	printf("%d out of $nn $dirs tests FAILED\n",$nn-$npassed);
    }
    else {
	print("All $nn $dirs tests PASSED\n");
    }
    chdir("..");
}

sub test_pdb2gmx {
    my $logfn = "pdb2gmx.log";

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
	foreach my $ff ( "G43a1", "oplsaa", "G53a6", "encads" ) {
	    mkdir("ff$ff");
	    chdir("ff$ff");
	    my @water = ();
	    my @vsite = ( "none", "h" );
	    if ( $ff eq "oplsaa"  ) {
		@water = ( "tip3p", "tip4p", "tip5p" );
	    }
	    elsif ( $ff eq "encads" ) {
		@vsite = ( "none" );
		@water = ( "spc" );
	    }
	    else {
		@water = ( "spc", "spce" );
	    }
	    foreach my $dd ( @vsite ) {
		mkdir("$dd");
		chdir("$dd");
		foreach my $ww ( @water ) {
		    $ntest++;
		    my $line = "";
		    print(LOG "****************************************************\n");
		    print(LOG "** PDB = $pdb FF = $ff VSITE = $dd WATER = $ww\n");
		    print(LOG "** Working directory = %s\n",`pwd`);
		    print(LOG "****************************************************\n");
		    mkdir("$ww");
		    chdir("$ww");
		    print(LOG "****************************************************\n");
		    print(LOG "**  Running pdb2gmx\n");
		    print(LOG "****************************************************\n");
		    open(PIPE,"$progs{'pdb2gmx'} -f ../../../../$pdb -ff $ff -ignh -vsite $dd -water $ww 2>&1 |");
		    print LOG while <PIPE>;
		    close PIPE;
		    print(LOG "****************************************************\n");
		    print(LOG "**  Running editconf\n");
		    print(LOG "****************************************************\n");
		    open(PIPE,"$progs{'editconf'} -o b4em -box 5 5 5 -c -f conf 2>&1 |");
		    print LOG while <PIPE>;
		    close PIPE;
		    print(LOG "****************************************************\n");
		    print(LOG "**  Running grompp\n");
		    print(LOG "****************************************************\n");
		    open(PIPE,"$progs{'grompp'} -maxwarn 3 -f ../../../../em -c b4em 2>&1 |");
		    print LOG while <PIPE>;
		    close PIPE;
		    print(LOG "****************************************************\n");
		    print(LOG "**  Running mdrun\n");
		    print(LOG "****************************************************\n");
		    open(PIPE,"$mdprefix $progs{'mdrun'} 2>&1 |");
		    print LOG while <PIPE>;
		    close PIPE;
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
	my $reflog = "${ref}.log";
	if (! -f $reflog) {
	    print "No file $reflog. You are not really testing pdb2gmx\n";
	    system("cp ener.log $reflog");
	}
	else {
	    $nerror = check_xvg($reflog,"ener.log",3,7);
	
	    if ( $nerror != 0 ) {
		print "There were $nerror differences in final energy with the reference file\n";
	    }
	}
	print "All $ntest pdb2gmx tests PASSED\n";
	system("rm -rf @pdb_dirs");
	unlink("ener.log","pdb2gmx.log");
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
sub test_gmx {
  foreach my $p ( values %progs ) {
    my $pp = $p;
    my $tgpp = `which $pp`;
    if (index($tgpp,$pp) < 0) {
      print("ERROR: Can not find $pp in your path.\nPlease source GMXRC and try again.\n");
      exit(1);
    }
  }
}

my $kk = 0;
my @work = ("setup_vars()", "test_gmx()");

for ($kk=0; ($kk <= $#ARGV); $kk++) {
    my $arg = $ARGV[$kk];
    if ($arg eq 'simple') {
	push @work, "test_dirs('simple')";
    }
    elsif ($arg eq 'complex') {
	push @work, "test_dirs('complex')";
    }
    elsif ($arg eq 'kernel' ) {
	push @work, "test_dirs('kernel')";
    }
    elsif ($arg eq 'pdb2gmx' ) {
	push @work, "test_pdb2gmx()";
    }
    elsif ($arg eq 'all' ) {
	push @work, "test_dirs('simple')";
	push @work, "test_dirs('complex')";
	push @work, "test_dirs('kernel')";
	push @work, "test_pdb2gmx()";
    }
    elsif ($arg eq 'clean' ) {
        clean_all();
    }
    elsif ($arg eq 'refclean' ) {
	push @work, "refcleandir('simple')";
	push @work, "refcleandir('complex')";
	push @work, "unlink('pdb2gmx/reference_s.log','pdb2gmx/reference_d.log')";
    }
    elsif ($arg eq 'dist' ) {
	push @work, "clean_all()";
	push @work, "chdir('..')";
	push @work, "system('tar --exclude CVS -czvf gmxtest.tgz gmxtest')";
	push @work, "chdir('gmxtest')";
    }
    elsif ($arg eq 'help' ) {
	push @work, "usage()";
    }
    elsif ($arg eq '-verbose') {
	$verbose++;
    }
    elsif ($arg eq '-noverbose') {
	$verbose = 0;
    }
    elsif ($arg eq '-double') {
	$double = 1;
    }
    elsif ($arg eq '-np' ) {
	if ($kk <$#ARGV) {
	    $kk++;
	    $parallel = $ARGV[$kk];
	    if ($parallel < 0) {
		$parallel = 0;
	    }
	    print "Will test on $parallel processors\n";
	}
    }
    elsif ($arg eq '-suffix' ) {
	if ($kk <$#ARGV) {
	    $kk++;
	    $suffix = $ARGV[$kk];
	    print "Will test using executable suffix $suffix\n";
	}
    }
    elsif ($arg eq '-prefix' ) {
	if ($kk <$#ARGV) {
	    $kk++;
	    $prefix = $ARGV[$kk];
	    print "Will test using executable prefix $prefix\n";
	}
    }
    else {
	push @work, "usage()";
    }
}

if ($kk == 0) {
    $#work = -1;
    push @work, "usage()";
}

if ( 1 == $#work ) {
    # there was no work added, so probably this was a gmxtest.pl clean
    # so don't do setup either
    $#work = -1;
}

# setup_vars() is always the first work to do, so now
# parallel and double will work correctly regardless of
# order on the command line
foreach my $w ( @work ) {
#    print "$w\n";
    eval $w;
}

if ($addversionnote > 0) {
    print << "ENDOFNOTE"

Note about GROMACS refernce versions
------------------------------------
Various different GROMACS versions are used to generate the reference
files for these tests. Because of a known bug with Buckingham
interactions in combination with LJ 1-4 interactions in GROMACS 3.3.x, 
the kernel_[0-3]2[0-4] test references are generated with GROMACS 4.0.5. All 
other kernel test references are generated with GROMACS 3.3. Most non-kernel
test references are generated with GROMACS 3.3.2. See the README file for
more detail.

If you are trying to test a version that precedes some of the above, then
this test set will not achieve your aim.

If you are testing a 3.3/3.3.1/3.3.2/3.3.3 version, then it is expected that 
all kernel_[0-3]2[0-4] tests fail, because of the known bug. This is only a 
problem if you wish to use Buckingham interactions with LJ 1-4 interactions, 
which is very rare. Most users can ignore this - or install GROMACS 4 
instead!

If you are seeing this message and neither of the above conditions is true, 
then you have another problem. If you post to the GROMACS mailing lists,
you must include the version number of GROMACS that you were testing, and
which tests failed and what they reported.
ENDOFNOTE
}
