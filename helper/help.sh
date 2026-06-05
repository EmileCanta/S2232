#!/bin/bash

HELPER=/home/emile/postdoc/analysis/s2232/helper/rcmp/RCMPHelper.cxx 
HELPERNAME=RCMP
HELPERDIR=/home/emile/postdoc/analysis/s2232/helper/rcmp
OPTIONS="--max-workers 20"
ANALYSISDIR=/home/emile/postdoc/data/s2232/sorted/analysis 
HISTDIR=/home/emile/postdoc/analysis/s2232/histograms/rcmp
LOGDIR=/home/emile/postdoc/analysis/s2232/helper/logs

firstRun=1 
lastRun=1

# Force rebuild if helper changed
if [ $HELPER -nt ${HELPERDIR}/${HELPERNAME}Helper.so ]; then
    echo "Rebuilding helper..."
    rm -f ${HELPERDIR}/${HELPERNAME}Helper.so ${HELPERDIR}/${HELPERNAME}Helper.o
fi

if [ $# -eq 0 ] ; then
   echo "usage: $0 <first run> <last run>"
   exit
fi

if [ $# -ge 1 ] ; then
   firstRun=$1
fi

if [ $# -ge 2 ] ; then
   lastRun=$2
else
   lastRun=$firstRun
fi

killproof() {
	procs=`pgrep proofserv.exe`

	if [[ $? -ne 0 ]] ; then
		exit
	fi
	# get all proof workers PID and PPID
	mapfile -t workers < <( ps --no-headers -fp $procs | tr -s ' ' | cut -d ' ' -f2,12 )

	# loop over all workers
	for t in ${!workers[@]} ; do
		# check if PPID exists
		if [ ! -e /proc/$(echo ${workers[$t]} | cut -d ' ' -f2) ] ; then
			# kill workers where parent does not exist
			kill -9 $(echo ${workers[$t]} | cut -d ' ' -f1)
		fi
	done
}

# loop over all runs from first to last
for run in `seq $firstRun $lastRun` ; do
   # check if the analysis-file of the first subrun exists
   if [ ! -e ${ANALYSISDIR}/analysis${run}_000.root ] ; then  #UNCOMMENT THIS UP TO THE fi LINE!!!
		echo "no ${ANALYSISDIR}/analysis${run}_000.root"
      continue
   fi
   # loop over all subruns that haven't been changed in the last three minutes
   # (only important when using this during an experiment)
   #for analysisFile in `find ${ANALYSISDIR} -maxdepth 1 -amin +3 -name "analysis${run}_???.root"` ; do
   analysisFiles=$(ls ${ANALYSISDIR}/analysis${run}_???.root)
   logFile="$LOGDIR/helper${run}.txt"
   # print the sort command to the log-file
   echo "grsiframe $OPTIONS $HELPER $analysisFiles" | tee $logFile
   # start sorting and append all output to the log-file
   # (2>&1 redirects stderr to stdout so we get error messages in the log file as well)
   grsiframe $OPTIONS $HELPER $analysisFiles 2>&1 | tee -a $logFile
   mv $HELPERNAME${run}_*.root $HISTDIR
   mv *.log $LOGDIR 
done
