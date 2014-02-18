

import os
import sys
import time
import subprocess
import shlex
import threading 
from dianaDir import *

#if 'DIANA_DIR' not in os.environ:
#	print "'DIANA_DIR' environmental variable not set, exiting"
#	exit(1)

#dianaDir 	= os.environ['DIANA_DIR']



# CONF READING
def readConf(path, dryrun = True):
	conf = {}
	f = open(path, "r")
	for line in f:
		if '=' in line:
			kv = line.split('=', 2)
			conf[kv[0].strip()] = kv[1].strip()
	f.close()
	
	if dryrun:
		for key in conf:
			print key, "->", conf[key]
	
	return conf


def require(conf, name):
	if name not in conf:
		print "Required parameter '%s' not found in conf file '%s'. Exiting" % (name, confPath)
		exit(2)
	else:
		return conf[name]

		
def default(conf, name, x):
	if name not in conf:
		return x
	else:
		return conf[name]




# LIST FILE
def readList(path):
	f = open(path, "r")
	items = [ x.strip() for x in f.readlines() ]
	f.close()
	return items



# FILE PATH MANIPULATION
def withoutExt(f):
	return ".".join(f.split(".")[:-1])

def base(f):
	return withoutExt(os.path.basename(f))

def to(f, ext):
	return "%s.%s" % (withoutExt(f), ext)



# COMMAND EXECUTION
INFO 	= ' INFO'
DEBUG 	= ' DEBUG'
WARNING = ' WARNI'
ERROR 	= ' ERROR'


def simpleExe(cmd, dryrun = True):
	if dryrun:
		print DEBUG, cmd
		return (0, "", "")
	try:
		print INFO, "Attempting to run   %s" % cmd
		shlexed = shlex.split(cmd)
		#retcode = os.system(" ".join(shlexed))
		p = subprocess.Popen(shlexed, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout, stderr = p.communicate()

		if p.returncode == 0:
			print "Successfully ran %s" % cmd
		else:
			print WARNING, "%s failed with exit code %d" % (cmd, p.returncode)
		return (p.returncode, stdout, stderr)
		
	except Exception as e:
		print ERROR, "Error running %s: %s" % (cmd, e)
		return (e.errno, "ERROR HAPPENED", "ERROR HAPPENED")


def execute(cmd, args, dryrun = True):
	def printOutput(program, outputType, output):
		if output != "":
			print DEBUG, "%s %s:" % (program, outputType)
			print output
	
	if dryrun:
		print DEBUG, cmd, args
		return 0
	else:
		try:
			print INFO, "Attempting to run   %s %s" % (cmd, args)
			shlexed = [cmd] + shlex.split(args)
			#retcode = os.system(" ".join(shlexed))
			p = subprocess.Popen(shlexed, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			stdout, stderr = p.communicate()

			printOutput(cmd, "stdout", stdout.replace("\r", ""))
			printOutput(cmd, "stderr", stderr.replace("\r", ""))

			if p.returncode == 0:
				print "Successfully ran %s" % cmd
			else:
				print WARNING, "%s failed with exit code %d" % (cmd, p.returncode)
			return p.returncode
			
		except Exception as e:
			print ERROR, "Error running %s: %s" % (cmd, e)
			return e.errno

def dianaBin(cmd):
	return "%s/bin/%s" % (dianaDir, cmd)
	
def dianaRun(cmd, args, dryrun = True):
	return execute(dianaBin(cmd), args, dryrun)

