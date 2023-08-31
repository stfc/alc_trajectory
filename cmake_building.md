# Building and testing
The following notes describe the steps to build **ALC_TRAJECTORY** using the [**CMake**](https://cmake.org) platform and Fortran compilers. In addition, we show how to run the testing infrastructure.

## Software required
The user must have access to the following software (locally):

* GNU-Fortran (7.2.0) or Intel-Fortran (16.0.1)
* Cmake (3.10.2)
* Make (3.82)
* Git (2.25.1)


Information in parenthesis indicates the minimum version tested during the development of the code. The specification for the minimum versions is not fully rigorous but indicative, as there could be combinations of other minimum versions that still work. Our tests indicate that versions of Intel compiler older than 16.0.1 exhibit problems and should be avoided.

## Building the code
We assume the code has been downloaded to folder ***alc_trajectory*** at the location */home/username/codes*, in the remote machine *wherever* of the account *"username"*. Please refer to file [README.md](./README.md) for instructions to download the code.

### Building manually
For manual compilation with CMake, the user must create a folder where to compile the code. It is good practise to name this folder using the word *"build"* together with any other specification that indicates the type of compilation. For example, if the user aims to build ALC_TRAJECTORY using the GNU-Fortran compiler in Debugging mode, the folder could be named *"build-gnu-debug"*:
```sh
username@wherever:/home/username/codes/alc_trajectory$ mkdir build-gnu-debug
username@wherever:/home/username/codes/alc_trajectory$ cd build-gnu-debug
username@wherever:/home/username/codes/alc_trajectory/build-gnu-debug$ FC=gfortran cmake ../ -DCMAKE_BUILD_TYPE=Debug
```
For the successful execution of the last step, the user must ensure to have access to the minimum version of the required software, as per specification above. In case the user opts to utilise the Inter compiler, *FC=ifort* must be set instead. If successful, the user will identify the following files:

CMakeCache.txt &nbsp; ***CMakeFiles*** &nbsp; Makefile &nbsp; ***bin*** &nbsp; cmake_install.cmake &nbsp; ***modules***

Finally, the user must compile the code as follows
```sh
username@wherever:/home/username/codes/alc_trajectory/build-gnu-debug$ make
```
If compilation is successful, the executable *alc_trajectory.x* will have been generated inside folder ***bin***. Extra implemented options are listed below:
* *-DBUILDER="string"*      (optional) String must be the name of the individual who builds the code.
* *-DWITH_TESTING=ON/OFF*. Must be set to "ON" for testing purposes (see below in section **Testing the code with CMake**). By default this option is set to OFF.

There are two available options for *-DCMAKE_BUILD_TYPE*, namely: *Debug* and *Release*. For code development purposes, we strongly recommend to compile the code with the *Debug* option, independently of the compiler. For the purposes o
f running the  code only, the user is advised to use the option *Release*. The pre-defined flags options for compilation are define in file *cmake/flags.cmake* and depend on the compiler, as we detailed in the following.

#### GNU-Fortran compiler
If the compiler is GNU-Fortran, the pre-defined compilation flags for the *Release* options are:
```sh
"-Ofast -ftree-vectorize -funroll-loops -ffast-math"
```
In contrast, if the selected option is *Debug*, the predefined compilation flags are:
```sh
"-g -Wextra -Wuse-without-only -frecord-gcc-switches -O0 -std=f2008 -pedantic -fbacktrace -fcheck=all -finit-integer=2147483647
-finit-real=snan -finit-logical=true -finit-character=42 -ffpe-trap=invalid,zero,overflow -fdump-core -fstack-protector-all -Wall -pipe"
```
for *gFortran* versions older than 6.5. If the compiler version is 7.2 or newer, the flag *-finit-derived* is also added to the list above to allow default initialization of derived-type variables.

#### Intel-Fortran compiler

For *Intel-Fortran* compiler (ifort), the pre-defined flags for option *Debug* are:
```sh
"-g -O0 -stand f08 -traceback -C -fp-stack-check -ftrapuv -init=snan -init=arrays"
```
whereas for the *Release* option, we have defined only the flag *-Ofast*.

### Building the code automatically
Inside folder **tools**, the user will find the following shell files (sh-files) that automatically build and compile the code:

* gnu-build-debug.sh: *Debug* option with *gFortran*. Generated folder is ***build-gnu-debug***.
* gnu-build.sh: &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Release* option with *gFortran*. Generated folder is ***build-gnu***.
* intel-build-debug.sh: *Debug* option with *Intel-Fortran*. Generated folder is ***build-intel-debug***.
* intel-build.sh:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; *Release* option with *Intel-Fortran*. Generated folder is ***build-intel***.

As an example, we will consider the sh-file *"gnu-build-debug.sh"* with the following content:
```sh
#!/usr/bin/env bash

folder="build-gnu-debug"
rm -rf $folder && mkdir $folder && cd $folder
FC=gfortran cmake ../  -DCMAKE_BUILD_TYPE=Debug  -DWITH_TESTING=Off
make
```

To execute this script, the user must proceed as follows:
```sh
username@wherever:/home/username/codes/alc_trajectory$ sh tools/gnu-build-debug.sh
```
which builds and compiles the code within folder ***build-gnu-debug***.

## Testing the code
### Building manually including the testing infrastructure
To have access to the testing infrastructure for ALC_TRAJECTORY, the user must add *-DWITH_TESTING=ON* when building with cmake. Based on the example above to building-only the code, the user must execute the following commands
```sh
username@wherever:/home/username/codes/alc_trajectory$ mkdir test-gnu-debug
username@wherever:/home/username/codes/alc_trajectory$ cd test-gnu-debug
username@wherever:/home/username/codes/alc_trajectory/test-gnu-debug$ FC=gfortran cmake ../ -DCMAKE_BUILD_TYPE=Debug -DWITH_TESTING=ON
```
where folder ***test-gnu-debug*** is generated to build, compile and test the code (any other name can be chosen).  If successful, the user will identify the following filing structure:

CMakeCache.txt &nbsp; ***CMakeFiles*** &nbsp; CTestTestfile.cmake &nbsp; DartConfiguration.tcl &nbsp; Makefile &nbsp; ***Testing*** &nbsp; ***bin*** &nbsp; cmake_install.cmake &nbsp; ***modules***

In addition to the files and folder for the building-only case above, the user will find:
* CTestTestfile.cmake: includes the relevant testing commands required for testing.
* DartConfiguration.tcl: describes the system on which tests are performed.
* ***Testing***: contains sub-folders ***new*** (where tests will be executed), ***reference*** (with the reference data for each test) and ***Temporary*** (for testing record purposes).

Results for each test run in ***new*** will be compared with reference data in ***reference***. Finally, the user must compile the code
```sh
username@wherever:/home/username/codes/alc_trajectory/test-gnu-debug$ make
```
If compilation is successful, executable *alc_trajectory.x* will have been generated in folder ***bin***. The list of available tests can be displayed by executing
```sh
username@wherever:/home/username/codes/alc_trajectory/test-gnu-debug$ ctest -N
```
### Running tests manually
Once the code has been compiled, all tests can be run using ctest (the testing tool of CMake) simply by executing
```sh
username@wherever:/home/username/codes/alc_trajectory/test-gnu-debug$ ctest
```
Each test will be executed using the executable *alc_trajectory.x* in each of the sub-folders ***Testing/new/testX***, where the generated files will be compared with the references files of sub-folders ***Testing/reference/testX***. A file diagnose.log is generated inside each ***Testing/new/testX*** sub-folder. This file provides a summary of the execution by reporting which of the generated files has/has not passed the test. To run a specific test, for example test10, the user must execute
```sh
username@wherever:/home/username/codes/alc_trajectory/test-gnu-debug$ ctest -R test10
```
Assuming we deliberately make test10 to fail, the generated output to screen will show as follows
```sh
Test project /home/username/codes/alc_trajectory/test-gnu-debug
    Start 10: test10
1/1 Test #10: test10 ...........................***Failed    0.43 sec

0% tests passed, 1 tests failed out of 1

Total Test time (real) =   0.47 sec

The following tests FAILED:
         10 - test10 (Failed)
Errors while running CTest
```
The user can opt to add *"--output-on-failure"* to print the outcome of *diagnose.log*. In case the test is completed but there is an inconsistency with the reference data on obtains
```sh
username@wherever:/home/username/codes/alc_trajectory/test-gnu-debug$ ctest --output-on-failure -R test10
Test project home/username/codes/alc_trajectory/test-gnu-debug
    Start 10: test10
1/1 Test #10: test10 ...........................***Failed    0.62 sec
*** Check for generated files against reference data ***
SUCCESS !!! OUTPUT passed the test
SUCCESS !!! TRACK_CHEMISTRY passed the test
FAILURE !!! RDF has NOT passed the test
SUCCESS !!! OCF passed the test
SUCCESS !!! OCF_AVG passed the test
SUCCESS !!! TCF passed the test
SUCCESS !!! TCF_AVG passed the test
SUCCESS !!! RES_TIMES passed the test
*** Test FAILED ****************************


0% tests passed, 1 tests failed out of 1

Total Test time (real) =   0.64 sec

The following tests FAILED:
         10 - test10 (Failed)
Errors while running CTest
```
which shows there was a problem with the data generated for the RDF file. If, in contrast, the test fails to run, one gets the following:
```sh
username@wherever:/home/username/codes/alc_trajectory/test-gnu-debug$ ctest --output-on-failure -R test10
Test project home/username/codes/alc_trajectory/test-gnu-debug
    Start 10: test10
1/1 Test #10: test10 ...........................***Failed    0.17 sec
*******************************************************
** ERROR !!! ERROR !!! ERROR !!! ERROR !!! ERROR !!! **
*******************************************************
**  Please see details in the OUTPUT. If the         **
**  error message makes no sense, check that the     **
**  input files are free of non-ASCII characters     **
*******************************************************
ERROR STOP
*** EXECUTION FAILED ***

0% tests passed, 1 tests failed out of 1

Total Test time (real) =   0.23 sec

The following tests FAILED:
         10 - test10 (Failed)
Errors while running CTest
```
where the error banner is also displayed. Option *--output-on-failure* can also be use when running all the tests with ctest.

### Building the code and running tests automatically
Inside folder ***tools***, the user will find the following sh-files that automatically i) build and compile the code and ii) run the whole set of tests:

* gnu-test-debug.sh: Building and compilation with the *Debug* option using *gFortran*. Tests run with ctest. Generated folder is ***test-gnu-debug***.
* gnu-test.sh: &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Building and compilation with the *Release* option using *gFortran*. Tests run with ctest. Generated folder is ***test-gnu***.
* intel-test-debug.sh: Building and compilation with the *Debug* option using *Intel-Fortran*. Tests run with ctest. Generated folder is ***test-intel-debug***.
* intel-test.sh:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Building and compilation with the *Release* option using *Intel-Fortran*. Tests run with ctest. Generated folder is ***test-intel***.

For all these sh-files, the option *"--output-on-failure"* is included for ctest, as explained before. As an example, we will consider the sh-file *gnu-test-debug.sh*, which has the following content:
```sh
#!/usr/bin/env bash
folder="test-gnu-debug"
rm -rf $folder && mkdir $folder && cd $folder
FC=gfortran cmake ../  -DCMAKE_BUILD_TYPE=Debug  -DWITH_TESTING=ON
make
ctest --output-on-failure
```
To execute this script, the user must proceed as follows:
```sh
username@wherever:/home/username/codes/alc_trajectory$ sh tools/gnu-test-debug.sh
```
which builds, compiles and runs the set of tests automatically.
