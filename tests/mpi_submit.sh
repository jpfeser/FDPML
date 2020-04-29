#!/bin/bash

#----------------------------------
#                Farber parameters
# tells cluster to assign (20) processors
#$ -pe mpi 20
# tells cluster to allocate 1GB memory PER processor
# $ -l m_mem_free=3G
# tells cluster to give you exclusive access to node
# $ -l exclusive=1
# send script to the standby queue
#$ -l standby=1,h_rt=4:00:00
#----------------------------------

#      Load any packages you need to run it
vpkg_require openmpi/intel64
vpkg_require gnuplot/4.6

TMP_DIR='/lustre/scratch/rohitk/FDPML/'

BIN_DIR='/home/1628/FDPML'
FDPML_bin='${BIN_DIR}/FDPML.out'
ref_bin='${BIN_DIR}/tests/calculate_error.out'

echo "	"
echo "Starting tests"
echo "	"

for NPROC in 1 5 10 20; do

OPENMPI_FLAGS='np -${NPROC}'
echo "Running first test on ${NPROC} processor(s)"
echo "	"
echo "-----GridEngine parameters-----"
echo "  mpirun        = "`which mpirun`
echo "  nhosts        = $NHOSTS"
echo "  nproc         = $NSLOTS"
echo "  executable    = $FDPML_bin"
echo "  MPI flags     = $OPENMPI_FLAGS"
echo "-------------------------------"

cat > input.txt << EOF
$(cat ${BIN_DIR}/tests/test_1/input.txt)
&restartoptions
	tmp_dir = ${TMP_DIR}
/
EOF

mpirun -quiet ${OPENMPI_FLAGS} $MY_EXE < input.txt > output.txt

echo "	"
printf "cleaning input and output files... "
rm input.txt output.txt

echo "done"

echo "	"

echo "checking output against reference solutions"

cat > input.txt << EOF
&filenames
	flfrc='${BIN_DIR}/tests/test_1/Si_q2.fc'
	tmp_dir='${TMP_DIR}'
	ref_filename='${BIN_DIR}/'
/
&system
	PD(1) = 11, 11, 20
	LPML = 20
	periodic = .true.
	crystal_coordinates = .false.
/
&simulation
	nprocs = ${NPROC}
/
EOF

./${ref_bin} < input.txt 

if [ $? -eq 0 ]; then
    echo "Test #1 passed on ${NPROC} processor(s)"
    echo "	"
else
    echo "Test #1 failed on ${NPROC} processor(s)"
	echo "Information on the failed test"
	echo "	"
	echo "$(cat ${BIN_DIR}/tests/test_1/error_message.txt)"
fi

echo "	"
echo "cleaning input and output files"

rm input.txt

done 
