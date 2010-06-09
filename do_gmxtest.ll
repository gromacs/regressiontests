#!/bin/bash
# @ job_name = thegmxtest
# @ step_name = mdrun
# @ output = $(job_name).$(step_name).$(cluster).out
# @ error = $(job_name).$(step_name).$(cluster).err
# @ wall_clock_limit = 0:30:00
# @ notification = complete
# @ job_type = bluegene
# @ class = bg32_100
# @ bg_size = 32
# @ group = bg01
# @ environment = COPY_ALL;
# @ queue

#export GMX_NOOPTIMIZEDKERNELS=1
#export NOASSEMBLYLOOPS=1
#export GMX_NB_GENERIC=1
#export GMX_NO_SOLV_OPT=1

# Find out the partition name so that we can use mpirun -nofree
# to avoid reallocating an identical partition for each mpirun
# invocation of mdrun for each test, which takes forever.
export MPIRUN_PARTITION=$(llq -b -j $LOADL_STEP_ID  | grep RMP | perl -nle  'print for m/\b(RMP\S*)\b/g')
export MPIRUN_NOFREE=1

#VERBOSE=-verbose -verbose
VERBOSE=-verbose
#DOUBLE=
#DOUBLE=-double
#TIGHT=-tight -tight
BLUEGENE=-bluegene
GMX_SUFFIX=

./gmxtest.pl -np 1 kernel ${VERBOSE} ${DOUBLE} ${BLUEGENE} -prefix ${HOME}/progs/bin/ -suffix ${GMX_SUFFIX} ${TIGHT}
./gmxtest.pl -np 1 simple ${VERBOSE} ${DOUBLE} ${BLUEGENE} -prefix ${HOME}/progs/bin/ -suffix ${GMX_SUFFIX} ${TIGHT}
./gmxtest.pl -np 1 complex ${VERBOSE} ${DOUBLE} ${BLUEGENE} -prefix ${HOME}/progs/bin/ -suffix ${GMX_SUFFIX} ${TIGHT}

# useful invocations for testing PME kernels
#./gmxtest.pl -np 1 kernel ${VERBOSE} ${DOUBLE} ${BLUEGENE} -prefix ${HOME}/progs/bin/ -suffix ${GMX_SUFFIX} ${TIGHT} -only "kernel(010|330|[3][1][012])"
#./gmxtest.pl -np 1 kernel ${VERBOSE} ${DOUBLE} ${BLUEGENE} -prefix ${HOME}/progs/bin/ -suffix ${GMX_SUFFIX} ${TIGHT} -only "kernel(010|330|[3][1][012])" -mdparam '-debug 2'
