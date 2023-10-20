#!/usr/bin/env bash

celltype=$1
workinggds=$2
output=$3
idatdir=$4
logfile=$5
samplesheet=$6

if [[ "${celltype}" == "brain" ]];
then
	Rscript scripts/4_minfi_estimatecelltypecomposition.r $logfile $idatdir $output $samplesheet
elif [[ "${celltype}" == "buccal" ]];
then
	Rscript scripts/4_epidish_estimatecelltypecomposition.r $logfile $workinggds $output $idatdir $samplesheet $celltype

elif [[ "${celltype}" == "saliva" ]];
then
	Rscript scripts/4_epidish_estimatecelltypecomposition.r $logfile $workinggds $output $idatdir $samplesheet $celltype

elif [[ "${celltype}" == "blood" ]];
then
	Rscript scripts/4_epidish_estimatecelltypecomposition.r $logfile $workinggds $output $idatdir $samplesheet $celltype
else
	echo "Cell Type estimation for this cell type has not been implemented yet."
fi