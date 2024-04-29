## About the code
**ALC_TRAJECTORY** is an open-source software that allows computing: 

* Orientational Correlation Functions (OCFs),
* Radial Distribution Functions (RDFs), and 
* Mean Square Displacements (MSDs) 

from Molecular Dynamics (MD) trajectories, both for reactive and non-reactive atomic systems ***all at once***. By "reactive" we refer to systems in gas and condensed phase where the constituent atomic species change their chemistry composition along the trajectory, forming and breaking bonds, as it occurs in anion  [[1]](https://pubs.acs.org/doi/10.1021/acs.jpcc.8b10298) and proton exchange membranes [[2]](https://pubs.acs.org/doi/10.1021/acs.jpclett.1c04071?ref=PDF), for example. For reactive systems **ALC_TRAJECTORY** also allow computing:

* the location of the changing chemical species,
* residence times, and
* Transfer Correlation Functions (TCFs). 

As for the OCF, RDF and MSD, these last capabilities can be computed  ***all at once***, thus offering a novel tool to analyse reactive systems where the changing chemical species play a decisive rol in the computation of these quantities. Moreover, the simple and flexible input structure is convenient to indentify chemical species in different type of systems and environments. Although the applicability of ***ALC_TRAJECTORY*** is general, its development and optimization have been focused to analyse systems within the size range of Density Functional Theory (DFT) simulations.  

The development of this code started in February 2023 at the Ada Lovelace Centre (ALC) of the Science and Technology Facilities Council (STFC). **ALC_TRAJECTORY** is a serial code written in modern Fortran according to the 2008 standards. Its structure for development (and maintenance) follows the Continuous Integration (CI) practice and it is integrated within the GitLab DevOps of the STFC.
In the root folder, the user will find several Markdown files, which are intended to provide help with the compilation and execution as well as guidance with the multiple available functionalities.  

## Disclaimer
The ALC does not fully guarantee the code is free of errors and assumes no legal responsibility for any incorrect outcome or loss of data.

## Contributors
### Original author
Ivan Scivetti (SCD, STFC)
### Scientific and project support
Gilberto Teobaldi (SCD, STFC)

## Structure of files and folders
ALC_TRAJECTORY contains the following set of files and folders (in italic-bold):

* [***CI-tests***](./CI-tests): contains the tests files (in .tar format) needed for CI purposes. The user should execute the available scripts of the [***tools***](./tools) folder to run the test automatically and verify the code has been installed properly (see the [cmake_building.md](./cmake_building.md) file for instructions).
* [***cmake***](./cmake): contains the specification for the compilation flags depending on the Fortran compiler and version, including options for debugging.
* [***examples***](./examples): example cases to help the user to become familiarised with the code. The SETTINGS files are described in detail.  
* [***scripts***](./scripts): contains scripts for data processing.
* [***source***](./source): contains the source code. Files have the *.F90* extension
* [***tools***](./tools): shell files for building, compiling and testing the code automatically.
* [.gitignore](./.gitignore): instructs Git which file to ignore for development and integration.
* [gitlab-ci.yml](gitlab-ci.yml): settings for automatic CI building and testing.
* [CI_instructions.md](./CI_instructions.md): instructions to new developers for CI.
* [CMakeList.txt](./CMakeList.txt): sets the framework for code building and testing with CMake. This file must ONLY be modified to add test cases.
* [LICENSE](./LICENSE): BSD 3-Clause License for ALC_TRAJECTORY. 
* README.md: this file.
* [cmake_building.md](./cmake_building.md): steps to build, compile and run tests using the CMake platform.
* [coding_protocol.md](./coding_protocol.md): details the instructions of the adopted protocol for code development with Fortran.
* [use_code.md](./use_code.md): provides instructions for use together with a detailed description of the implemented capabilties. 

## Dependencies
The user must have access to the following software (locally):

* GNU-Fortran (7.2.0) or Intel-Fortran (16.0.1)
* Cmake (3.10.2)
* Make (3.82)
* git (2.25.1)

Information in parenthesis indicates the minimum version tested during the development of the code. The specification for the minimum versions is not fully rigorous but indicative, as there could be combinations of other minimum versions that still work. Our tests indicate that versions of Intel compiler older than 16.0.1 exhibit problems and should be avoided.

## Getting started

### Obtaining the code
The user with account *"username"* can clone the code locally (in machine *"wherever"*) by executing the following command with the SSH protocol
```sh
username@wherever:/home/username/codes$ git clone git@gitlab.stfc.ac.uk:alc_trajectory/alc_trajectory.git
```
which generate the ***alc_trajectory*** folder as the root directory. Alternatively, the code can be downloaded from any of the available assets.

### Building and testing the code with CMake
Details can be found in file [cmake_building.md](./cmake_building.md)

### Making use of the software
Once the code has been installed and tested, the user should create a folder where to run the code from. In such folder, the MD trajectory must be copied to the TRAJECTORY file. The user also needs to provide the SETTINGS file with instructions for the type of analysis to execute. Instructions of the implemented capabilities can be found in the [use_code.md](./use_code.md) file. The SETTINGS files in folder [***examples***](./examples) offer explanatory templates, which are intended to help new users in the setting of input directives for execution. In each of the directories, the user will also find the corresponding TRAJECTORY files.

