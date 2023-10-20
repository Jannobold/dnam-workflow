#!/usr/bin/env bash

logfile=$1
input1=$2
workinggds=$3
rawunfiltered=$4
rawfiltered=$5
outlyx=$6
bscon=$7
outliernames=$8
input2=$9
rawfilteredall=${10}
normunfiltered=${11}
qualout=${12}

normoutlier=true
while [ "$normoutlier" = true ]
do
	qcoutlier=true
	while [ "$qcoutlier" = true ]
	do
		if [[ ! -f $outlyx ]];
		then
			Rscript scripts/1_bigmelon_workflow1_loadqc.r $logfile $input1 $workinggds $rawunfiltered $rawfilteredall $outlyx $bscon $outliernames
			if [[ -s $outliernames ]]; 
			then
				while read p
				do
					plength=${#p}
					p=${p:1:$plength-2}
					mv ${input1}/${p}* ${input2} 
				done < "$outliernames"
				cp $rawfilteredall $rawfiltered
			else
				qcoutlier=false
				cp $rawfilteredall $rawfiltered
			fi
		else
			rm $workinggds
			rm $rawunfiltered
			rm $rawfiltered
			Rscript scripts/1_bigmelon_workflow1_loadqc.r $logfile $input1 $workinggds $rawunfiltered $rawfiltered $outlyx $bscon $outliernames
			if [[ -s $outliernames ]]; 
			then
				while read p
				do
					plength=${#p}
					p=${p:1:$plength-2}
					mv ${input1}/${p}* ${input2} 
				done < "$outliernames"
			else
				qcoutlier=false
			fi
		fi
	done
	Rscript scripts/2_bigmelon_workflow2_normqc.r $logfile $workinggds $normunfiltered $rawfiltered $qualout
	if [[ -s $qualout ]];
	then
		while read p
		do
			plength=${#p}
			p=${p:1:$plength-2}
			mv ${input1}/${p}* ${input2}
		done < "$qualout"
		rm $normunfiltered
	else
		normoutlier=false
	fi
done