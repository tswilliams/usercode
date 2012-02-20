#!/bin/bash
qsub  $*

while [ "$?" != "0" ]; do 
   echo " ** qsub had non-zero exit code. Trying to resubmit ... ** " >&2
   qsub $*
done
