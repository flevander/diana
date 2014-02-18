
import sys
import os 

dianaDir = None
for a in sys.argv:
	if a.lower().startswith("--diana-dir="):
		dianaDir = a[12:]

if dianaDir is None:
	if 'DIANA_DIR' not in os.environ:
		print "  'DIANA_DIR' environmental variable not set, and not given on command line. Exiting..."
		exit(1)
	else:
		dianaDir = os.environ['DIANA_DIR']
