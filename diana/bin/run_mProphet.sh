#!/bin/bash
#
#	usage:
#	run_mProphet.sh FILE_BASE
#
#

if [ $# == 0 ] 
then
	echo "usage: 
 run_mProphet.sh FILE_BASE"
	exit
fi

FILE_BASE=$1

./prep4mProphet.py $FILE_BASE.csv $FILE_BASE.decoy.csv $FILE_BASE.irt.csv


R --slave --args \
	bin_dir=./mProphet/ \
	data_file=$FILE_BASE.csv.to.mProph \
	workflow=LABEL_FREE \
	num_xval=5 \
	run_log=FALSE \
	write_classifier=1 \
	write_all_pg=1 \
	help=0 \
	project=$FILE_BASE \
	< ./mProphet/mProphet.R

csv2esv -f ${FILE_BASE}_raw_stat.xls
esv clean ${FILE_BASE}_raw_stat.xls.esv filtered.esv 'qvalue < 0.01'
