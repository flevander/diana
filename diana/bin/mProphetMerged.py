#!/usr/bin/env python

import os
import sys
import time
from datetime import datetime
from dianaLib import *

dryrun = False

if 'DIANA_DIR' not in os.environ:
        print "  'DIANA_DIR' environmental variable not set. Exiting..."
        exit(1)

dianaDir = os.environ['DIANA_DIR']

usage = """ mProphetMerged.py SWATH_LIST_FILE PROJECT_NAME"""

if len(sys.argv) < 3:
	print usage
	exit(1)
else:
	swathListFile   = sys.argv[-2]
	project 		= sys.argv[-1]


swathFiles = []
if swathListFile.lower().endswith("wiff") or swathListFile.lower().endswith("raw"):
	swathFiles.append(swathListFile)
else:
	swathFiles = readList(swathListFile)

mergedMprophIn = "merged.mProph.in.csv"
out = open(mergedMprophIn, "w")

writtenHeader = False

for swathFile in swathFiles:
	fileBase = base(swathFile)
	sf = open("%s/%s_all_targets.to.mProph.csv" % (fileBase, fileBase), "r")
	
	header = sf.readline()
	if not writtenHeader:
		out.write("%s\tfile\n" % (header[:-1]))
		writtenHeader = True
	
	for line in sf.readlines():
		out.write("%s|%s\t%s\n" % (fileBase, line[:-1], fileBase))
	
	sf.close()

out.close()




args = 	"--slave --args bin_dir=%s/bin/mProphet/ data_file=%s " % (dianaDir, mergedMprophIn)
args += "workflow=LABEL_FREE num_xval=5 run_log=FALSE write_classifier=1 "
args += "write_all_pg=0 help=0 project=%s < %s/bin/mProphet/mProphet.R" % (project, dianaDir)

if dryrun:
	print "    DEBUG      "
	print "   =======     "
	print " wanted to run:"
	print "R", args
else:
	os.system("R " + args)
