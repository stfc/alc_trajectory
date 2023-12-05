#!/bin/bash

### ####################################################################################################
## Bash file to run a test. Function CompareFiles compares the generated 
## output files with the reference output files cpied in the reference directory
## #####################################################################################################

  CompareFile()
{
  ## Function to compare files. In case of OUTPUT, banner and appendix are removed from the comparison
  local genfile="$1"
  local ref="$2"

  if [  -f ${fileref} ]; then

    if [ $genfile = "OUTPUT" ] ;then
      nl=$(wc -l ${genfile} | awk '{ print $1 }')
      nappex=10
      nbanner=14
      nlow="$((nl - nappex))"
      ntop="$((nlow - nbanner))"
      (head -n $nlow $genfile | tail -n $ntop) &> out1.txt
      (head -n $nlow $fileref | tail -n $ntop) &> out2.txt
      diff out1.txt out2.txt  &> diff.log
      rm out1.txt out2.txt
    else
      diff $genfile $fileref &> diff.log
    fi

    if [ ! -s diff.log ]; then
      echo  "SUCCESS !!! ${genfile} passed the test"
    else
      echo  "FAILURE !!! ${genfile} has NOT passed the test"
    fi
    rm diff.log
  fi
}

# # Define the list of output files
   declare -a FileArray=("OUTPUT" "TAGGED_TRAJECTORY" \\
                         "TRACK_CHEMISTRY"  "UNCHANGED_CHEMISTRY" \\
                         "RDF" "MSD" "MSD_AVG" "OCF" "OCF_AVG"  \\
                         "TCF" "TCF_AVG" "RES_TIMES" "COORD_DISTRIBUTION" \\
                         "INTRAMOL_DISTANCES" "INTRAMOL_ANGLES" "DISTANCE_SHORTEST_PAIR" \\
                         "INTERMOL_DISTANCES_NN1" "INTERMOL_DISTANCES_NN2" "INTERMOL_ANGLES_NN")

  sep="/"	
  error="success"

  # Delete all output files before they are computed. This is because the test includes the
  # reference outputs files inside the testX.tar, which is also copied (and unpacked) to the "reference" directory.
  for i in "${FileArray[@]}"; do
    filename=${i}
    if [ -f ${filename} ]; then
      rm -rf $filename
    fi
  done

## Execute the job
$1

## Check if job has run ot not
if [ $? -ne 0 ]; then
  status="FAILED"
  echo "*** EXECUTION FAILED ***"  | tee -a diagnose.log 
  exit 125
else
  status="PASSED"
fi

## If the jobs for has executed, it is time to check the output with reference data
if [ "$status" = "PASSED" ]; then

  echo "*** Check for generated files against reference data ***"  | tee -a diagnose.log

  for i in "${FileArray[@]}"; do
    genfile=${i}
    fileref="$2$sep${i}"

    if [ -f ${genfile} ]; then
      CompareFile "$genfile" "$fileref" | tee -a diagnose.log
      check=$(awk '{w=$1} END{print w}' diagnose.log)
      if [ "$check" = "FAILURE" ] ; then
        error="fail"
      fi
    else
      if [ -f ${fileref} ]; then
        echo  "ERROR !! ${genfile} has NOT been generated"
        error="fail"
      fi
    fi
  done

  if [ $error = "fail" ] ; then
    echo "*** Test FAILED ****************************"  | tee -a diagnose.log
    exit 125
  else
    echo "*** Test was SUCCESSFUL ********************"  | tee -a diagnose.log
  fi  

fi
