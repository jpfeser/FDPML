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

# specify and create your directories (Please dont change these)
TMP_DIR='scratch/'

mkdir -p ${TMP_DIR}


printf "cleaning temporary scratch files... "

rm scratch/*

echo "done"

BIN_DIR='/home/1628/FDPML'
FDPML_bin=${BIN_DIR}/bin/FDPML.out
ref_bin=${BIN_DIR}/tests/bin/calculate_error.out

echo "	"
echo "Starting tests"
echo "	"

# Run through all the folders and run each test on 1, 5, 10 and 220 processors

for test in test_*/ ; do

for NPROC in 1 5 10 20; do

OPENMPI_FLAGS="-np $NPROC"

echo "	"

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
$(cat ${BIN_DIR}/tests/${test}/input_FDPML.txt)
EOF

mpirun -quiet ${OPENMPI_FLAGS} $FDPML_bin < input.txt > output.txt

echo "	"
printf "cleaning input and output files... "
rm input.txt output.txt

echo "done"

echo "	"

sleep 10

echo "checking output against reference solutions"

cat > input.txt << EOF
$(cat ${BIN_DIR}/tests/${test}/input_error.txt)
&simulation
	nprocs = ${NPROC}
/
EOF

mpirun -quiet -np 1 ${ref_bin} < input.txt

# check if the run was successful

if [ $? -eq 0 ]; then
    echo "Test #1 passed on ${NPROC} processor(s)"
    echo "	"
else
    echo "Test #1 failed on ${NPROC} processor(s)"
	echo "Information on the failed test"
	echo "	"
	echo "$(cat ${BIN_DIR}/tests/${test}/error_message.txt)"
fi

echo "	"
printf "cleaning input and output files... "

rm input.txt

echo "done"

echo "	"

printf "cleaning temporary scratch files... "

rm scratch/*

echo "done"

sleep 10

done 

done
