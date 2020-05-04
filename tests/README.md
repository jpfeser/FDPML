## Unit testing / Example scripts for FDPML

This module runs FDPML script for some base cases and compares the solutions against reference soltions

### Installation

Make sure the intel fortran is in your path. Run `make all` to compile the script

### Run

To execute just run submit.sh script.

### Creating additional tests (Instructions to developers)

To create additional tests 

1. Create a folder named test_<test no>
1. Create a README.md that describes the test
1. Upload any force constant files, domain specification files ( `mass.dat` and `domain.dat` ) files in this folder
1. Create a text file with FDPML input cards called `input_FDPML.txt` (rename the `tmp_dir` variable to `scratch/`)
1. Create a text file for calculate_error.f90 code named `input_error.txt` (for reference check `test_1/input_error.txt`)
1. Create a text file names `error_message.txt`. This is the error that will be printed out on console log if the test fails
1. Upload a reference solution for the test. Ideally, run the test on FDPML with a single processor and move the solution from your tmp_dir to this folder
 
