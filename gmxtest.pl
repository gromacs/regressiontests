#!/usr/bin/perl -w

use strict;

my $parallel = 0;
my $double   = 0;
my $bluegene = 0;
my $verbose  = 5;
my $xml      = 0;
my $etol     = 0.05;
my $ttol     = 0.001;
my $suffix   = '';
my $autosuffix   = '1';
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

# global variables to flag whether to explain some situations to the user
my $addversionnote = 0;
my $only_subdir;
my $tightfactor = 1;

# trickery for program and reference file names
my $mdprefix = 'valgrind --suppressions=/home/rschulz/download/gromacs5.0/cmake/legacy_and_external.supp --xml=yes --xml-file=valgrind.xml';
my $mdparams = '';
my $ref      = '';
my $mpirun   = 'mpirun';
my %progs = ( 'grompp'   => 'grompp',
	      'mdrun'    => 'mdrun',
	      'pdb2gmx'  => 'pdb2gmx',
	      'gmxcheck' => 'gmxcheck',
	      'editconf' => 'editconf' );
sub setup_vars()
{
    # We assume that the name of executables match the pattern 
    # ${prefix}mdrun[_mpi][_d]${suffix} where ${prefix} and ${suffix} are 
    # as defined above (oro ver-ridden on the command line), "_d" indicates 
    # a double-precision version, and (only in the case of mdrun) "_mpi" 
    # indicates a parallel version compiled with MPI.
    if ( $parallel > 0 ) {
	if ($autosuffix) {
	    $progs{'mdrun'} .= "_mpi";
	}
	if ( $bluegene > 0 )
	{
	    # edit the next line if you need to customize the call to mpirun
	    $mdprefix = "$mpirun -np $parallel -exp_env GMX_NO_SOLV_OPT -exp_env NOASSEMBLYLOOPS -exp_env GMX_NOOPTIMIZEDKERNELS -exp_env GMX_NB_GENERIC";
	} else {
	    # edit the next line if you need to customize the call to mpirun
	    $mdprefix = "$mpirun -np $parallel -wdir `pwd`";
	}
    }
    foreach my $prog ( values %progs ) {
	$prog = $prefix . $prog;
	if ($autosuffix) {
	    $prog .= "_d" if ( $double > 0 );
	}
	$prog .= $suffix;
    }
    $ref = 'reference_' . ($double > 0 ? 'd' : 's');
    
    # now do -tight stuff
    foreach my $var ( $etol, $ttol, $virtol_rel, $ftol_rel, $ftol_sprod ) {
	$var *= $tightfactor;
    }
}

# Wrapper function to call system(), and then perhaps a callback based on the
# value of the return code. When no callback exists, a generic error is
# displayed
sub do_system
{
    my $command = shift;
    my $normalreturn = shift;
    my $callback = shift;
    $normalreturn = 0 unless(defined $normalreturn);

    if ( $verbose > 2 ) {
	print "$command\n";
    }
    my $returnvalue = system($command) >> 8;
    if ($normalreturn != $returnvalue)
    {
	if (defined $callback)
	{
	    &$callback($returnvalue);
	}
	else 
	{
	    print "\nAbnormal return value for '$command' was $returnvalue\n";
	}
    }
    return $returnvalue;
}

sub check_force()
{
    my $tmp = "checkforce.tmp";
    my $cfor = "checkforce.out";
    my $reftrr = "${ref}.trr";
    my $nerr_force = 0;
    do_system("$progs{'gmxcheck'} -f $reftrr -f2 traj -tol $ftol_rel > $cfor 2> /dev/null", 0,
	      sub { print "\ngmxcheck failed on the .edr file while checking the forces\n"; $nerr_force = 1; });
    `grep "f\\[" $cfor > $tmp`;
    
    open(FIN,"$tmp");
    while(my $line=<FIN>)
    {
	my @ll=split("[()]",$line);
	my @f1=split(" ",$ll[1]);
	my @f2=split(" ",$ll[3]);
	
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
    my $nerr_vir = 0;

    do_system("$progs{'gmxcheck'} -e $refedr -e2 ener -tol $virtol_rel -lastener Vir-ZZ > $cvir 2> /dev/null", 0,
	sub { print "\ngmxcheck failed on the .edr file while checking the virial\n"; $nerr_vir = 1; });
    
    `grep "Vir-" $cvir > $tmp`;
    
    open(VIN,"$tmp");
    while(my $line=<VIN>)
    {
	my @v1=split(" ",substr($line,26,14)); #TODO replace substr with split to make more reliable
	my @v2=split(" ",substr($line,52,13)); #if again reactiving check_virial
	
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
    my $pdb2gmx_test_names = shift;

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
		my $error;
		my $hasPdb2gmx_test_name = defined $pdb2gmx_test_names && defined $$pdb2gmx_test_names[$n];
		print XML "<testcase name=\"$$pdb2gmx_test_names[$n]\">\n" if ($xml && $hasPdb2gmx_test_name);
		my $tol;
		if ($x1+$x2==0) { 
		    $error = abs($x1-$x2);
		    $tol = $ttol;
		} else {
		    $error = abs(($x1-$x2)/($x1+$x2));
		    $tol = $etol;
		}
		if ($error > $tol) {
		    $nerr++;
		    if (!$header) {
			$header = 1;
			print("Here follows a list of the lines in $refx and $kkk which did not\npass the comparison test within tolerance $etol\nIndex  Reference   This test       Error  Description\n");
		    }
		    printf("%4d  %10g  %10g  %10g  %s\n",$n+1,$tmp[3],$tmp[7], $error, 
			   $hasPdb2gmx_test_name ? $$pdb2gmx_test_names[$n] : 'unknown');
		    printf(XML "<error message=\"Reference: %g Result: %g Error: %g\"/>\n",$tmp[3],$tmp[7], $error) 
			if ($xml && $hasPdb2gmx_test_name);
		    
		}
		print XML "</testcase>\n" if ($xml && $hasPdb2gmx_test_name);
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
	    next if ( defined $only_subdir && $dir !~ $only_subdir);
	    chdir($dir);
	    if ($verbose > 1) {
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
	    do_system("$progs{'grompp'} -maxwarn 10 $ndx > grompp.out 2>&1");
	    
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
		    link('topol.tpr', $reftpr);
		} else {
		    do_system("$progs{'gmxcheck'} -s1 $reftpr -s2 topol.tpr -tol $ttol > checktpr.out 2>&1", 0, 
			sub { print "Comparison of input .tpr files failed!\n"; $nerror = 1; });
		    $nerror |= `grep step checktpr.out | grep -v gcq | wc -l`;
		    if ($nerror > 0) {
			print "topol.tpr file different from $reftpr. Check files in $dir\n";
		    }
		    do_system("grep 'reading tpx file (reference_[sd].tpr) version .* with version .* program' checktpr.out > /dev/null 2>&1", 1,
			      sub {
				  print "\nThe GROMACS version being tested may be older than the reference version.\nPlease see the note at end of this output.\n";
				  $addversionnote = 1;
			      } );
		}
	    }
	    my @error_detail;
	    if ($nerror == 0) {
		# Do the mdrun at last!

		# BlueGene mpirun needs to be told the current working
		# directory on the command line, or with env var MPIRUN_CWD,
		# so after the chdir we need to deal with this.
		my $local_mdprefix = $mdprefix;
		if ($bluegene > 0 ) {
		    $local_mdprefix .= " -cwd `pwd`";
		}
		my $local_mdparams = $mdparams;
		if (system("grep ns_type.*simple grompp.mdp > /dev/null")==0) {
		    $local_mdparams .= " -pd"
		}
		$nerror = do_system("$local_mdprefix $progs{'mdrun'} $local_mdparams > mdrun.out 2>&1", 0,
		    sub { push(@error_detail, ("mdrun.out", "md.log")); } );
		
		# First check whether we have any output
		if ((-f "ener.edr" ) && (-f "traj.trr")) {
		    # Now check whether we have any reference files
		    my $refedr = "${ref}.edr";
		    if (! -f  $refedr) {
			print ("No $refedr file in $dir.\n");
			print ("This means you are not really testing $dir\n");
			link('ener.edr', $refedr);
		    } else {
			# Now do the real tests
			do_system("$progs{'gmxcheck'} -e $refedr -e2 ener -tol $etol -lastener Potential > checkpot.out 2> /dev/null", 0,
				  sub {
				      if($nerror != 0) {
					  print "\ngmxcheck failed on the .edr file, probably because mdrun also failed";
				      }
				  });
			my $nerr_pot = `grep step checkpot.out | grep -v gcq | wc -l`;
			chomp($nerr_pot);
			push(@error_detail, "checkpot.out ($nerr_pot errors)") if ($nerr_pot > 0);

			my $nerr_vir   = 0; #TODO: check_virial();
			push(@error_detail, "checkvir.out ($nerr_vir errors)") if ($nerr_vir > 0);

			$nerror |= $nerr_pot | $nerr_vir;
		    }
		    my $reftrr = "${ref}.trr";
		    if (! -f $reftrr ) {
			print ("No $reftrr file in $dir.\n");
			print ("This means you are not really testing $dir\n");
			link('traj.trr', $reftrr);
		    } else {
			# Now do the real tests
			my $nerr_force = check_force();
			push(@error_detail, "checkforce.out ($nerr_force errors)") if ($nerr_force > 0);
			$nerror |= $nerr_force;
		    }
		    # This bit below is only relevant for free energy tests
		    my $refxvg = "${ref}.xvg";
		    my $nerr_xvg = check_xvg($refxvg,'dgdl.xvg',1,3);
		    push(@error_detail, "$refxvg ($nerr_xvg errors)") if ($nerr_xvg > 0);
		    $nerror |= $nerr_xvg;
		}
		else {
		    $nerror = 1;
		}
		$error_detail = join(', ', @error_detail) . ' ';
	    }
	    print XML "<testcase name=\"$dir\">\n" if ($xml);
	    if ($nerror > 0) {
		print "FAILED. Check ${error_detail}files in $dir\n";
		if ($xml) {
		    print XML "<error message=\"Erorrs in ${error_detail}\">\n";
		    print XML "<![CDATA[\n";
		    foreach my $err (@error_detail) {
			my @err = split(/ /, $err);
			my $errfn = $err[0];
			print XML "$errfn:\n";
			open FH, $errfn or die("FAILED to open $errfn");
			while(my $line=<FH>) {
			    print XML $line;
			}
			print XML "\n--------------------------------\n";
			close FH;
		    }
		    print XML "]]>\n";
		    print XML "</error>";
		}
	    }
	    else {
		my @args = glob("#*# *.out topol.tpr confout.gro ener.edr md.log traj.trr");
		#unlink(@args);
		
		if ($verbose > 0) {
		    if (`cat grompp.mdp | wc -l` < 50) {
			# if the input .mdp file is trivially short, then 
			# the diff test below will always fail, however this
			# is normal and expected for the usefully-short
			# kernel test .mdp files, so we don't compare the
			# .mdp files in this case
			print "PASSED\n";
		    }
		    else {
			my $mdp_result = 0;
			foreach my $reference_mdp ( 'grompp.mdp', 'grompp4.mdp', 'grompp41.mdp' ) {
			    if (-f $reference_mdp && `diff $reference_mdp mdout.mdp | grep -v host | grep -v date | grep -v user | grep -v generated | wc -l` > 2 ) {
				$mdp_result++;
			    } 
			    else {
				$mdp_result -= 100;
				last;
			    }
			}
			if ($mdp_result > 0) {
			    print("PASSED but check mdp file differences\n");
			}
			else {
			    print "PASSED\n";
			    unlink("mdout.mdp");
			}
		    }
		}
		$npassed++;
	    }
	    print XML "</testcase>\n" if ($xml);
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
	    my @args = glob("#*# *~ *.out core.* field.xvg dgdl.xvg topol.tpr confout.gro ener.edr md.log traj.trr *.tmp mdout.mdp step*.pdb *~ grompp[A-z]* *.cpt" );
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
    print XML "<testsuite name=\"$dirs\">\n" if ($xml);
    my $npassed = test_systems(@subdirs);
    print XML "</testsuite>\n" if ($xml);
    if ($npassed < $nn) {
	printf("%d out of $nn $dirs tests FAILED\n",$nn-$npassed);
    }
    else {
	print("All $nn $dirs tests PASSED\n");
    }
    chdir("..");
}

#format for cfg files is:
#- command line 
#- list of output files which should be compared (one per line)
#- emtpy line
sub test_tools {
    chdir("tools");
    my $ncfg = 0;
    my $nerror_cfg = 0;

    foreach my $cfg ( glob("*.cfg") ) { #loop over config files
	$ncfg++;
	open(FIN,$cfg);
	my @cfg_name = split(".cfg",$cfg);
	my $cfg_name = $cfg_name[0];
	mkdir($cfg_name);
	chdir($cfg_name);
	my $ncmd = 0;
	my $nerror_cmd = 0;
	if ($verbose > 1) {
	    print "Testing $cfg_name . . . \n";
	}
	while(my $line=<FIN>) {  #loop over commands (seperated by empty line)
	    $ncmd++;
	    chomp($line);
	    my $cmdline = $line;
	    my @ofiles;
	    my $error = 0;
	    while(my $line=<FIN>) {  #add output fiels
		chomp($line);
		if ($line eq "") { last; }
		push(@ofiles, $line);
	    }
	    mkdir($ncmd);
	    chdir($ncmd);

	    do_system("$cmdline > $cfg_name.out 2>&1", 0,
		      sub { print "\n'".$cmdline."' failed"; $error = 1; });
	    
	    if ($error==0) {
		foreach my $of (@ofiles) {
		    if (! -f "$of.ref") {
			print "No file $of.ref. You are not really testing $cfg_name\n";
			link($of, "$of.ref");
		    }
		    else {
			my $nerror = check_xvg("$of.ref",$of,1,3);  #TODO: check all columns
			if ( $nerror != 0 ) {
			    print "There were $nerror differences in $of output file\n";
			    $error += 1;
			}
		    }
		}
	    }
	    if ($error > 0) {
		print "FAILED $cfg_name test $ncmd\n";
		$nerror_cfg++;
		$nerror_cmd++;
	    }
	    chdir("..");
	}
	if ($nerror_cmd>0) {
	    print "$nerror_cmd out of $ncmd $cfg_name tests FAILED\n";
	} elsif ($verbose>1) {
	    print "All $ncmd $cfg_name tests PASSED\n";
	}
	chdir("..");
    }
    if ($nerror_cfg>0) {
	print "$nerror_cfg out of $ncfg tools tests FAILED\n";
    } else {
	print "All $ncfg tools tests PASSED\n";
    }
}

sub test_pdb2gmx {
    my $logfn = "pdb2gmx.log";

    chdir("pdb2gmx");
    open (LOG,">$logfn") || die("FAILED: Opening $logfn for writing");
    my $npdb_dir = 0;
    my @pdb_dirs = ();
    my $ntest    = 0;
    my $nerror   = 0;
    my @pdb2gmx_test_names;
    foreach my $pdb ( glob("*.pdb") ) {
	my $pdir = "pdb-$pdb";
	my @kkk  = split('\.',$pdir);
	my $dir  = $kkk[0];
	$pdb_dirs[$npdb_dir++] = $dir;
	mkdir($dir);
	chdir($dir);
	foreach my $ff ( "gromos43a1", "oplsaa", "gromos53a6" ) {
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
		    push @pdb2gmx_test_names, "$pdb with $ff using vsite=$dd and water=$ww";
		    my $line = "";
		    print(LOG "****************************************************\n");
		    print(LOG "** PDB = $pdb FF = $ff VSITE = $dd WATER = $ww\n");
		    printf(LOG "** Working directory = %s\n",`pwd`);
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
		    open(PIPE,"$mdprefix $progs{'mdrun'} $mdparams 2>&1 |");
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
    
    do_system("grep 'Potential Energy' pdb2gmx.log > ener.log");
    
    my $nsuccess = `wc -l ener.log | awk '{print \$1}'`;
    chomp($nsuccess);
    if ( $nsuccess != $ntest ) {
	print "Error not all $ntest pdb2gmx tests have been done successfully\n";
	print "Only $nsuccess energies in the log file\n";
	$nerror = 1;
    }
    else {
	my $reflog = "${ref}.log";
	if (! -f $reflog) {
	    print "No file $reflog. You are not really testing pdb2gmx\n";
	    link('ener.log', $reflog);
	}
	else {
	    print XML "<testsuite name=\"pdb2gmx\">\n" if ($xml);
	    $nerror = check_xvg($reflog,"ener.log",3,7,\@pdb2gmx_test_names);
	    print XML "</testsuite>\n" if ($xml);
	    if ( $nerror != 0 ) {
		print "There were $nerror/$ntest differences in final energy with the reference file\n";
	    }
	}
    }
    if (0 == $nerror) {
	print "All $ntest pdb2gmx tests PASSED\n";
	`rm -rf @pdb_dirs`;
	unlink("ener.log","pdb2gmx.log");
    }
    else {
	print "pdb2gmx tests FAILED\n";
    }
    chdir("..");
}

sub clean_all {
    cleandirs("simple");
    cleandirs("complex");
    cleandirs("kernel");
    cleandirs("freeenergy");
    chdir("pdb2gmx");
    unlink("ener.log", "pdb2gmx.log");
    `rm -rf pdb-*`;
    chdir("..");
}

sub usage {
    print "Usage: ./gmxtest.pl [ -np N ] [-verbose ] [ -double ] [ -bluegene ]\n [ -prefix xxx ] [ -suffix xxx ] [ -reprod ]\n [ simple | complex | kernel | freeenergy | pdb2gmx | all ]\n";
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
    elsif ($arg eq 'freeenergy' ) {
	push @work, "test_dirs('freeenergy')";
    }
    elsif ($arg eq 'pdb2gmx' ) {
	push @work, "test_pdb2gmx()";
    }
    elsif ($arg eq 'tools' ) {
	push @work, "test_tools()";
    }
    elsif ($arg eq 'all' ) {
	push @work, "test_dirs('simple')";
	push @work, "test_dirs('complex')";
	push @work, "test_dirs('kernel')";
	push @work, "test_dirs('freeenergy')";
	push @work, "test_pdb2gmx()";
	push @work, "test_tools()";
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
	push @work, "system('tar --exclude CVS --exclude .git --exclude .gitattributes --exclude foreach.sh -czvf regressiontests.tgz regressiontests')";
	push @work, "chdir('regressiontests')";
    }
    elsif ($arg eq 'help' ) {
	usage();
    }
    elsif ($arg eq '-verbose') {
	$verbose++;
    }
    elsif ($arg eq '-noverbose') {
	$verbose = 0;
    }
    elsif ($arg eq '-xml') {
	$verbose = 0;
        $xml = 1;
    }
    elsif ($arg eq '-double') {
	$double = 1;
    }
    elsif ($arg eq '-bluegene') {
	$bluegene = 1;
	print "Will test BlueGene\n";
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
    elsif ($arg eq '-nosuffix' ) {
	$autosuffix = 0;
    }
    elsif ($arg eq '-prefix' ) {
	if ($kk <$#ARGV) {
	    $kk++;
	    $prefix = $ARGV[$kk];
	    print "Will test using executable prefix $prefix\n";
	}
    }
    elsif ($arg eq '-mpirun' ) {
	if ($kk <$#ARGV) {
	    $kk++;
	    $mpirun = $ARGV[$kk];
	}
    }
    elsif ($arg eq '-mdparam' ) {
	# The user is responsible for providing sensible values for
	# this flag when they want them
	if ($kk <$#ARGV) {
	    $kk++;
	    $mdparams .= $ARGV[$kk];
	    print "Will test using 'mdrun $mdparams'\n";
	}
    }
    elsif ($arg eq '-only' ) {
	# only test a subdirectory if it matches the following
	# regular expression
	if ($kk <$#ARGV) {
	    $kk++;
	    print "Will only test subdirectories matching '$ARGV[$kk]'\n";
	    $only_subdir = qr/$ARGV[$kk]/;
	}
    }
    elsif ($arg eq '-tight' ) {
	$tightfactor *= 0.1;
	print "Will test with tightness increased\n";
    }
    else {
	usage();
    }
}

if ($kk == 0) {
    $#work = -1;
    usage();
}

if ( 1 == $#work ) {
    # there was no work added, so probably this was a 'gmxtest.pl clean'
    # so don't do setup either
    $#work = -1;
}

if ($xml) {
    my $xmlfn = "gmxtest.xml";
    open (XML,">$xmlfn") || die("FAILED: Opening $xmlfn for writing");
    print XML '<?xml version="1.1" encoding="UTF-8"?>';
    print XML "\n<testsuites>\n";
}
# setup_vars() is always the first work to do, so now
# parallel and double will work correctly regardless of
# order on the command line
map { eval $_ } @work;

print XML "</testsuites>\n" if ($xml);

if ($addversionnote > 0) {
    print << "ENDOFNOTE"

Note about GROMACS reference versions
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
