import argparse
import sys
import os
import shutil

def check_arg(args = None):
	'''Parses command line arguments.'''
	parser = argparse.ArgumentParser(description = 'Pipeline for scRNA-Seq analysis.')
	parser.add_argument('-d', '--directory',
		help = 'Directory where data and scripts are stored',
		required = 'True'
		)	
	parser.add_argument('-l', '--logFileName',
		help = 'Desired Log File Name',
		required = 'True'
		)
	return parser.parse_args(args)

def enterDirectory(directory):
  os.chdir(directory)
  
def createLog(logFileName):
  '''Creates an empty log file where run info will be stored.'''
  logFile = open(logFileName, 'w')
  return logFile

def callSeurat(seuratScript):
  '''Function to call the Seurat script, which is an Rscript.'''
  callingScript = f'Rscript {seuratScript}'
  logFile.write('Running Seurat... ')
  os.system(callingScript)
  logFile.write('COMPLETE')
  logFile.write('\n-------------------------------------------\n')
  
  
# retrieving command line arguments and assigning to variables
args = check_arg(sys.argv[1:])
directory = args.directory
logFileName = args.logFileName

#initializeOutputDirectory(outputDir)
enterDirectory(directory)
logFile = createLog(logFileName)
seurat = callSeurat('seurat.R')

