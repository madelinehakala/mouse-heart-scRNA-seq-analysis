import argparse
import sys
import os
import shutil
import csv
import pandas as pd
import seaborn as sns
import subprocess

def check_arg(args = None):
	'''Parses command line arguments.'''
	parser = argparse.ArgumentParser(description = 'Pipeline for scRNA-Seq analysis.')
	parser.add_argument('-d', '--directory',
		help = 'Directory where data and scripts are stored',
		required = 'True'
		)	
	parser.add_argument('-g', '--genesOfInterestFile',
		help = 'File name for a text file that contains the genes of interest you would like to include in the final dataset',
		required = 'True'
		)	
	parser.add_argument('-s', '--sampleLabel',
		help = 'Provide a label for your samples (ex: heartZone2)',
		required = 'True'
		)	
	parser.add_argument('-l', '--logFileName',
		help = 'Desired Log File Name',
		required = 'True'
		)
	parser.add_argument('-g1', '--markerGene1',
		help = 'Marker gene 1 (required)',
		required = 'False'
		)
	parser.add_argument('-g2', '--markerGene2',
		help = 'Marker gene 2 (not a required argument)',
		required = 'False'
		)
	parser.add_argument('-f', '--filteredOutfileName',
		help = 'Desired Name for the Filtered Outfile (with selected cell population of interest)',
		required = 'False'
		)	
	parser.add_argument('-p', '--plotFileName',
		help = 'Desired Name for the Expression Plot File',
		required = 'False'
		)	
	parser.add_argument('-n', '--normalizedOutfileName',
		help = 'Desired Name for the final outfile with normalized data',
		required = 'False'
		)	
	return parser.parse_args(args)

  
def createLog(logFileName):
  '''Creates an empty log file where run info will be stored.'''
  logFile = open(logFileName, 'w')
  logFile.write('-------------------------------------------\n')
  return logFile

def callSeurat(seuratScript):
  '''Function to call the Seurat script, which is an Rscript.'''
  logFile.write('Running Seurat... ')
  subprocess.call(['Rscript', seuratScript, directory, sampleLabel, markerGene1, markerGene2])
  logFile.write('COMPLETE\n')
  logFile.write('-------------------------------------------\n')
  seuratOutfile = f'{sampleLabel}_expression_data.csv'
  return seuratOutfile

def selectCells(seuratOutfile, filteredOutfileName, genesOfInterest, markerGene1, markerGene2):
  '''Filters the outputed file from Seurat to include only the cell type of interest. Additionally, it selects only the columns with genes of interest. Outputs this data into a new csv.'''
  if markerGene2 == None:
    print('View the Rplots.pdf file for expression visualizations of your marker genes of interest to assess necessary minimum expression cutoffs.')
    filterValueA = float(input(f'Set Minimum Expression Level of {markerGene1}: '))
    expressionData = pd.read_csv(seuratOutfile)
    logFile.write(f'{seuratOutfile} prior to filtering...\n\n')
    logFile.write(f'{expressionData.head()}\n\n')
    logFile.write('-------------------------------------------\n')
    genesOfInterestFile = open(genesOfInterest, 'r')
    genesOfInterestData = genesOfInterestFile.read()
    genesOfInterestList = genesOfInterestData.split(',')
    expressedGenesOfInterestList = []
    for gene in genesOfInterestList:
      if gene in expressionData.columns:
        expressedGenesOfInterestList.append(gene)
    expressionDataGenesOfInterest = expressionData[expressedGenesOfInterestList]
    expressionDataGenesOfInterest = expressionDataGenesOfInterest.loc[expressionDataGenesOfInterest[markerGene1] >= filterValueA]
    expressionDataGenesOfInterest.to_csv(filteredOutfileName, index = False)
    logFile.write(f'{filteredOutfileName} after filtering...\n\n')
    logFile.write(f'{expressionDataGenesOfInterest.head()}\n')
    logFile.write('-------------------------------------------\n')
  else:
    print('View the Rplots.pdf file for expression visualizations of your marker genes of interest to assess necessary minimum expression cutoffs.')
    filterValueA = float(input(f'Set Minimum Expression Level of {markerGene1}: '))
    filterValueB = float(input(f'Set Minimum Expression Level of {markerGene2}: '))
    expressionData = pd.read_csv(seuratOutfile)
    logFile.write(f'{seuratOutfile} prior to filtering...\n\n')
    logFile.write(f'{expressionData.head()}\n\n')
    genesOfInterestFile = open(genesOfInterest, 'r')
    genesOfInterestData = genesOfInterestFile.read()
    genesOfInterestList = genesOfInterestData.split(',')
    expressedGenesOfInterestList = []
    for gene in genesOfInterestList:
      if gene in expressionData.columns:
        expressedGenesOfInterestList.append(gene)
    expressionDataGenesOfInterest = expressionData[expressedGenesOfInterestList]
    expressionDataGenesOfInterest = expressionDataGenesOfInterest.loc[expressionDataGenesOfInterest[markerGene1] >= filterValueA]
    expressionDataGenesOfInterest = expressionDataGenesOfInterest.loc[expressionDataGenesOfInterest[markerGene1] >= filterValueB]
    expressionDataGenesOfInterest.to_csv(filteredOutfileName, index = False)
    logFile.write(f'{filteredOutfileName} after filtering...\n\n')
    logFile.write(f'{expressionDataGenesOfInterest.head()}\n')
    logFile.write('-------------------------------------------\n')
  
def visualizeExpressionData(filteredOutfile, plotFileName):
  '''Visualizes expression data.'''
  filteredExpressionData = pd.read_csv(filteredOutfile)
  sns.set_theme(rc = {'figure.figsize': (20, 6)})
  boxplot = sns.boxplot(filteredExpressionData)
  boxplot.set_xticklabels(boxplot.get_xticklabels(), rotation = 45)
  boxplot.set(ylabel = 'Expression Level')
  figure = boxplot.get_figure()
  figure.savefig(plotFileName)
  
def normalizeExpressionData(filteredOutfile, normalizedOutfileName):
  '''Normalizes the expression data by dividing each individual value by the population average for that gene.'''
  filteredExpressionData = pd.read_csv(filteredOutfile)
  for column in filteredExpressionData.columns:
    try:
      values = filteredExpressionData[column]
      mean = filteredExpressionData[column].mean()
      if mean != 0:
        filteredExpressionData[column] = values / mean
      else:
        filteredExpressionData[column] = values
    except TypeError:
      continue
  filteredExpressionData.to_csv(normalizedOutfileName, index = False)
  logFile.write(f'{normalizedOutfileName} after normalizing expression levels...\n\n')
  logFile.write(f'{filteredExpressionData.head()}\n')
  logFile.write('-------------------------------------------\n')
  
# retrieving command line arguments and assigning to variables
args = check_arg(sys.argv[1:])
directory = args.directory
genesOfInterest = args.genesOfInterestFile
sampleLabel = args.sampleLabel
logFileName = args.logFileName
markerGene1 = args.markerGene1
try:
  markerGene2 = args.markerGene2
except:
  markerGene2 = None
filteredOutfileName = args.filteredOutfileName
plotFileName = args.plotFileName
normalizedOutfileName = args.normalizedOutfileName

# calling functions
logFile = createLog(logFileName)
seurat = callSeurat('seurat.R')
cellsOfInterest = selectCells(seurat, filteredOutfileName, genesOfInterest, markerGene1, markerGene2)
preNormalizationVisualization = visualizeExpressionData(filteredOutfileName, plotFileName)
cellsOfInterestNormalized = normalizeExpressionData(filteredOutfileName, normalizedOutfileName)
postNormalizationVisualization = visualizeExpressionData(normalizedOutfileName, f'normalized_{plotFileName}')
