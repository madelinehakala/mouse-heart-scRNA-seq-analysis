import argparse
import sys
import os
import shutil

def check_arg(args = None):
	'''Parses command line arguments.'''
	parser = argparse.ArgumentParser(description = 'Pipeline for scRNA-Seq analysis.')
	parser.add_argument('-o', '--outputDirectory',
		help = 'Directory where all outputed files and folders will be stored',
		required = 'True'
		)	
	parser.add_argument('-l', '--logFileName',
		help = 'Desired Log File Name',
		required = 'True'
		)
	return parser.parse_args(args)

def initializeOutputDirectory(directory):
  '''Initializes a directory where all outputted files, folders, etc will be stored.'''
  if os.path.exists(directory): # if the directory with that name already exists, remove it
    shutil.rmtree(directory)
  os.makedirs(directory) # then make the directory and move into it
  os.chdir(directory)
  
def createLog(logFileName):
  '''Creates an empty log file where run info will be stored.'''
  logFile = open(logFileName, 'w')
  return logFile

def callSeurat(seuratScript):
  '''Function to call the Seurat script, which is an Rscript.'''
  callingScript = f'Rscript {seuratScript}'
  os.system(callingScript)
  logFile.write('Running Seurat... COMPLETE')
  logFile.write('\n-------------------------------------------\n')
  
  
# retrieving command line arguments and assigning to variables
args = check_arg(sys.argv[1:])
outputDir = args.outputDirectory
logFileName = args.logFileName

#initializeOutputDirectory(outputDir)
logFile = createLog(logFileName)
seurat = callSeurat('seurat.R')
