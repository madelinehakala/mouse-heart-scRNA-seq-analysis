import argparse
import sys
import os
import shutil
import csv
import pandas as pd

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
  logFile.write('-------------------------------------------\n')
  return logFile

def callSeurat(seuratScript):
  '''Function to call the Seurat script, which is an Rscript.'''
  callingScript = f'Rscript {seuratScript}'
  logFile.write('Running Seurat... ')
  os.system(callingScript)
  logFile.write('COMPLETE\n')
  logFile.write('-------------------------------------------\n')

def selectMuscleCells(seuratOutfile, filteredOutfileName, filterValueA, filterValueB):
  expressionData = pd.read_csv(seuratOutfile)
  logFile.write(f'{seuratOutfile} prior to filtering...\n\n')
  logFile.write(f'{expressionData.head()}\n\n')
  muscleContractionFile = open('muscle_contraction_genes.txt', 'r')
  muscleContractionData = muscleContractionFile.read()
  muscleContractionGeneList = muscleContractionData.split(',')
  expressedMuscleContractionGeneList = []
  for gene in muscleContractionGeneList:
    if gene in expressionData.columns:
      expressedMuscleContractionGeneList.append(gene)
  expressionDataMuscleContractionGenes = expressionData[expressedMuscleContractionGeneList]
  expressionDataMuscleContractionGenes = expressionDataMuscleContractionGenes.loc[expressionDataMuscleContractionGenes['Myh6'] >= filterValueA]
  expressionDataMuscleContractionGenes = expressionDataMuscleContractionGenes.loc[expressionDataMuscleContractionGenes['Actc1'] >= filterValueB]
  expressionDataMuscleContractionGenes.to_csv(filteredOutfileName, index = False)
  logFile.write(f'{seuratOutfile} after filtering...\n\n')
  logFile.write(f'{expressionDataMuscleContractionGenes.head()}\n')
  logFile.write('-------------------------------------------\n')
  
  
# retrieving command line arguments and assigning to variables
args = check_arg(sys.argv[1:])
directory = args.directory
logFileName = args.logFileName

#initializeOutputDirectory(outputDir)
enterDirectory(directory)
logFile = createLog(logFileName)
seurat = callSeurat('seurat.R')
zone2MuscleCells = selectMuscleCells('zone2_expression_data.csv', 'zone2_expression_data_genes_of_interest.csv', 0.75, 3)
zone3MuscleCells = selectMuscleCells('zone3_expression_data.csv', 'zone3_expression_data_genes_of_interest.csv', 0.75, 3)

