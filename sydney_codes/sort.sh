DATADIR=/media/sydney/GRINYER_29/MidasFiles #/media/sydney/RCMP_SG/S2232_midas_files
SORTOPTIONS="--recommended"
CALFILES=/home/sydney/S2232_analysis/cal_files/S2232Cal.cal  #~ is shortcut for /home/sydney
ANALYSISDIR=/media/sydney/GRINYER_29/S2232_analysis/AnalysisTrees
FRAGMENTDIR=/media/sydney/GRINYER_29/S2232_analysis/FragmentTrees
LOGDIR=/media/sydney/GRINYER_29/S2232_analysis/Logs

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

# loop over all runs from first to last
for run in `seq $firstRun $lastRun` ; do
	# check if the midas-file of the first subrun exists
	if [ ! -e ${DATADIR}/run${run}_000.mid ] ; then
		echo "${DATADIR}/run${run}_000.mid does not exist"
		continue
	fi
	# loop over all subruns that haven't been changed in the last three minutes
	# (only important when using this during an experiment)
	for midasFile in `find ${DATADIR} -maxdepth 1 -name "run${run}_???.mid"` ; do
		if [ ! -e $midasFile ] ; then
			continue
		fi
		subrun=$(basename $midasFile | cut -d '_' -f2 | cut -d '.' -f1)
		analysisFile="$ANALYSISDIR/analysis${run}_${subrun}.root"
		# if the analysis file exists, we don't re-run the analysis
		if [ -e $analysisFile ] ; then
			if [ $analysisFile -nt $midasFile ] ; then
				continue
			fi
		fi
		# touching the file means it exists now and if we run this script in multiple terminal
		# the other scripts won't try and start sorting this file
		touch $analysisFile
		logFile="$LOGDIR/log${run}_${subrun}.txt"
		# print the sort command to the log-file
		echo "grsisort $SORTOPTIONS $CALFILES $midasFile" | tee $logFile
		# start sorting and append all output to the log-file
		# (2>&1 redirects stderr to stdout so we get error messages in the log file as well)
		grsisort $SORTOPTIONS $CALFILES $midasFile 2>&1 | tee -a $logFile
		mv fragment${run}_${subrun}.root $FRAGMENTDIR
		mv analysis${run}_${subrun}.root $ANALYSISDIR
	done
done
