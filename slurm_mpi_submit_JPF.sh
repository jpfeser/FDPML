#!/bin/bash -l

# Template:  OpenMPI, Default (OpenIB Infiniband) Variant
# Revision:  $Id: openmpi-ib.qs 577 2015-12-21 14:39:43Z frey $
#
# Usage:
# 1. Modify "NPROC" in the -pe line to reflect the number
#    of processors desired.
# 2. Modify the value of "MY_EXE" to be your MPI program and
#    "MY_EXE_ARGS" to be the array of arguments to be passed to
#    it.
# 3. Uncomment the WANT_CPU_AFFINITY line if you want Open MPI to
#    bind workers to processor cores.  Can increase your program's
#    efficiency.
# 4. Uncomment the SHOW_MPI_DEBUGGING line if you want very verbose
#    output written to the Grid Engine output file by OpenMPI.
# 5. If you use exclusive=1, please be aware that NPROC will be
#    rounded up to a multiple of 20.  In this case, set the
#    WANT_NPROC variable to the actual core count you want.  The
#    script will "load balance" that core count across the N nodes
#    the job receives.
# 6. Jobs default to using 1 GB of system memory per slot.  If you
#    need more than that, set the m_mem_free complex.
# 




#----------------------------------
#                Farber parameters
#SBATCH --partition=standard
# SBATCH --partition=_workgroup_
#
# CPU Allocations
#SBATCH --array=1-1
#SBATCH --ntasks=36
# OR BOTH
# SBATCH --nodes=<nhosts>
# SBATCH --tasks-per-node=<nproc-per-node>
#  this is for multithreading (for mpi only use 1)
#SBATCH --cpus-per-task=1
#
# MEMORY Allocations (suffix K/M/G/T), default 1G per cpu
# SBATCH --mem=8G
#SBATCH --mem-per-cpu=1024M
# scratch space required (optional)
# SBATCH --tmp=24G
#
# TIME Allocation (d-hh:mm:ss format)
#SBATCH --time=0-01:00:00

# SBATCH --job-name=openmpi_job

#SBATCH --mail-user='9735252392@vtext.com'
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_90
#----------------------------------

#      Load any packages you need to run it
vpkg_require openmpi/4.0.2:intel
vpkg_require gnuplot
OPENMPI_FLAGS='-np 36' # DANGER. must match scheduler flag above.

# note: choose num. of tasks = NLlist*Nqlist*Nslist

file_input='.false.'
#default values
Llist_file=Llist.csv
slist_file=slist.csv
qlist_file=qlist.csv
nL=1
ns=1
nq=1
L=1
s=1
qp=1

echo "My Task ID is : " ${SLURM_ARRAY_TASK_ID}

if [ ${file_input} == '.true.' ]; then 
	nL=`head -1 $Llist_file`
	ns=`head -1 $slist_file`
	nq=`head -1 $qlist_file`
	
	echo "Lpoint = " $nL
	echo "spoint = " $ns
	echo "qpoint = " $nq
	
	# use task ID to determine current value of L, s, qp
	L=`expr ${SLURM_ARRAY_TASK_ID} % $nL`
	s=`expr ${SLURM_ARRAY_TASK_ID} / $nL`
	
	s=`expr $s + 1`
	
	if [ $L == 0 ] ; then
		L=$nL
		s=`expr $s - 1`
	fi
	
	qp=`expr $s % $nq`
	s=`expr $s / $nq`
	
	s=`expr $s + 1`
	
	if [ $qp == 0 ] ; then
		qp=$nq
		s=`expr $s - 1`
	fi

fi


#
# The MPI program to execute:
#
BIN_DIR='/home/1627/fdpml_learning/FDPML/bin'
MY_EXE="FDPML.out"

TMP_DIR='/lustre/scratch/jpfeser/FDPML_tmp'
DOMAIN="${TMP_DIR}/Domain.dat"
MASS="${TMP_DIR}/mass.dat"

echo "SLURM parameters:"
echo "  mpirun        = "`which mpirun`
echo "  nhosts        = $NHOSTS"
echo "  nproc         = $NSLOTS"
echo "  executable    = $MY_EXE"
echo "  MPI flags     = $OPENMPI_FLAGS"
echo "-- begin OPENMPI run --"

# use timing for convergence
rm *plot


#	NAMELIST /filenames/ mass_file, domain_file, flfrc1, flfrc2, mass_input &
#			 /system/ simulation_type, PD, LPML, periodic, crystal_coordinates, &
#					  asr, wavetype, q, mode, sigmamax, expense_estimate &
#			 /qlists/ q_from_file, q_file &
#			 /solver/ tol, maxit &
#			 /restartoptions/ restart, iterpause, tmp_dir, &
#			 /postprocessing/ calc_TC, scattered_energy, scattering_Xsec &
#			 /plots/ plot_K, plot_sig, plot_uinc, plot_uscat, plottingmode &
#			 /calibrate/ qlist_file, slist_file, Llist_file, Lpoint, spoint, file_input, &
							qpoint

cat > incard.${SLURM_ARRAY_TASK_ID} << EOF
  &filenames
	mass_file = '${MASS}',
	domain_file = '${DOMAIN}',
	flfrc1 = '/home/1627/fdpml_learning/gendomain/GaAs444.fc', 
	flfrc2 = '/home/1627/fdpml_learning/gendomain/GaAs444_massErAs.fc',
	mass_input = .true.
/
  &system
	simulation_type = 'nanoparticle',
	PD(1) = 11, 11, 11
	periodic = .false.
	crystal_coordinates = .false.
	asr = 'simple'
	wavetype = 'full'
	q(1) = 0.0, 0.0, 0.2
	mode = 3
    LPML = 10
	sigmamax = 10
/
   &qlists 
/
  &solver
	tol = 3.5e-8
	maxit = 1000000
/
  &restartoptions
	tmp_dir='${TMP_DIR}' 
/
  &postprocessing
	calc_TC = .true.
	scattered_energy= .false.
/
  &plots
	plot_K = .false.
	plot_sig = .false.
	plot_uinc = .false.
	plot_uscat = .true.
	plottingmode = 3
/
	&calibrate
	file_input = ${file_input}
	qlist_file = '${qlist_file}'
	Llist_file = '${Llist_file}'
	slist_file = '${slist_file}'
	Lpoint=$L
	spoint=$s
	qpoint=$qp
/ 
EOF


mpirun -quiet ${OPENMPI_FLAGS} $BIN_DIR/$MY_EXE < incard.${SLURM_ARRAY_TASK_ID} > outcard.${SLURM_ARRAY_TASK_ID}

echo "Hello"

# check to see if it completed without errors
rc=$?

exit $rc

echo "-- DONE --"
