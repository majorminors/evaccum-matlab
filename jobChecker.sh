#!/bin/bash
directory=$1

runChecks () {

for job in $directory/Job*/; do
    echo "

==============================
this is:                      
                              
$job                          
==============================

"
   echo "
   log:
   "
   # search for words that might indicated a problem
   if tail $job/Task1.log | grep -i -q  "exit"; then
       echo "
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  FINISHED  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"
   fi
   if tail $job/Task1.log | grep -i -q "error"; then
       echo "
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!  ERROR  !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
"
   fi
   # print the tail of the log
   tail $job/Task1.log
   echo "
   diary:
   "
   # print the tail of the diary
   tail $job/Task1.diary.txt
   echo "

------------------------------
finito                             
------------------------------

"
done

}


if [[ "$2" == "view" ]]; then
    runChecks
elif [[ "$2" == "check" ]]; then
    # only output if there is a warning or error
    if runChecks | grep -i -q "error"; then echo "there is an errored job"; fi
    if runChecks | grep -i -q "finished"; then echo "there is a finished job"; fi
    echo "there is nothing else to report"
fi

