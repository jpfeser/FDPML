#!/bin/bash

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
# tells cluster to assign (20) processors
# $ -pe mpi 1
# tells cluster to allocate 1GB memory PER processor
# $ -l m_mem_free=3G
# tells cluster to give you exclusive access to node
#$ -l exclusive=1
# send script to the standby queue
# $ -l standby=1,h_rt=4:00:00
# setup messaging about (b)egin, (e)nd, (a)bort, and (s) of your program
# $ -o mpi_submit_script.sh.out
# $ -m beas
# send messages to this email address
# -M 9735252392@vtext.com
#----------------------------------

#      Load any packages you need to run it
vpkg_require openmpi/intel64
vpkg_require gnuplot/4.6
OPENMPI_FLAGS='-np 20'

L=`expr $SGE_TASK_ID % 11`
s=`expr $SGE_TASK_ID / 11`

s=`expr $s + 1`

if [ $L == 0 ] ; then
		L=11
        s=`expr $s - 1`
fi

qp=`expr $s % 50`
s=`expr $s / 50`

s=`expr $s + 1`

if [ $qp == 0 ] ; then
	qp=50
	s=`expr $s - 1`
fi

echo "Lpoint = " $L
echo "spoint = " $s
echo "qpoint = " $qp

PMLL=20

#
# The MPI program to execute:
#
MY_EXE="FDPML.out"

TMP_DIR='/lustre/scratch/rohitk/FDPML_tmp'
DOMAIN="${TMP_DIR}/Domain.dat"
MASS="${TMP_DIR}/mass.dat"

echo "GridEngine parameters:"
echo "  mpirun        = "`which mpirun`
echo "  nhosts        = $NHOSTS"
echo "  nproc         = $NSLOTS"
echo "  executable    = $MY_EXE"
echo "  MPI flags     = $OPENMPI_FLAGS"
echo "-- begin OPENMPI run --"

cat > test_input.${SGE_TASK_ID} << EOF
  &filenames
	domain_file = '${DOMAIN}',
	mass_file = '${MASS}', 
	flfrc1 = '/home/1628/QuantumEspresso/GaAs/results/GaAs444.fc', 
	flfrc2 = '/home/1628/QuantumEspresso/GaAs/results/GaAs444.fc'
	mass_input = .false.
/
  &system
	simulation_type = 'interface'
	PD(1) = 11, 11, 20
	LPML = ${PMLL}
	periodic = .true.
	crystal_coordinates = .false.
	asr = 'simple'
	wavetype = 'half'
	q(1) = 0.0, 0.0, 0.1
	mode = 3
	sigmamax = 5
/
  &qlists
	q_from_file = .false.
	q_file = 'q_file.csv'
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
	scattered_energy= .true.
	scattering_Xsec = .true.
/
  &plots
	plot_K = .false.
	plot_sig = .false.
	plot_uinc = .false.
	plot_uscat = .true.
	plottingmode = 3
/
  &calibrate
	file_input = .false.
	qlist_file = 'qlist.csv'
	Llist_file = 'Llist.csv'
	slist_file = 'slist.csv'
	Lpoint=$L
	spoint=$s
	qpoint=$qp
/
EOF



mpirun -quiet ${OPENMPI_FLAGS} $MY_EXE < test_input.${SGE_TASK_ID} > output_test.o${SGE_TASK_ID}

echo "Hello"

# check to see if it completed without errors
rc=$?

exit $rc

echo "-- DONE --"
