#!/bin/bash

HELPER=/home/sydney/S2232_analysis/RCMPHelpers/RCMPHelper.cxx #Change for each histogram helper
HELPERNAME=RCMP #Same as the prefix in your your .hh file!! 
OPTIONS="--max-workers 31" #Change the 31 to whatever number of cores you want your computer
ANALYSISDIR=/media/sydney/GRINYER_28/S2232_analysis/AnalysisTrees #Direct to where your analysis root files are! #/home/sydney/Desktop/28226E    
HISTDIR=/home/sydney/S2232_analysis/Histograms/RCMP
LOGDIR=/home/sydney/S2232_analysis/Logs/HistogramLogs

firstRun=1 
lastRun=1

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
done
