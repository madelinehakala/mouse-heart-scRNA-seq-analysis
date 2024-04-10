import argparse
import sys
import os
import shutil
import csv
import pandas as pd
import seaborn as sns

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
  '''Filters the outputed file from Seurat to include only muscle cells. Additionally, it selects only the columns with muscle contraction genes of interest. Outputs this data into a new csv.'''
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
  logFile.write(f'{filteredOutfileName} after filtering...\n\n')
  logFile.write(f'{expressionDataMuscleContractionGenes.head()}\n')
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
  
def normalizeExpressionData(filteredOutfile, outfileName):
  '''Calculates the average expression level of each muscle marker gene and outputs to a new file.'''
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
  filteredExpressionData.to_csv(outfileName, index = False)
  logFile.write(f'{outfileName} after normalizing expression levels...\n\n')
  logFile.write(f'{filteredExpressionData.head()}\n')
  logFile.write('-------------------------------------------\n')
  
# retrieving command line arguments and assigning to variables
args = check_arg(sys.argv[1:])
directory = args.directory
logFileName = args.logFileName

enterDirectory(directory)
logFile = createLog(logFileName)
seurat = callSeurat('seurat.R')
zone2MuscleCells = selectMuscleCells('zone2_expression_data.csv', 'zone2_muscle_cell_expression_data.csv', 0.75, 3) # 0.75 is minimum expression level for Myh6
zone3MuscleCells = selectMuscleCells('zone3_expression_data.csv', 'zone3_muscle_cell_expression_data.csv', 0.75, 3) # 3 is minimum expression level for Actc1
zone2Visualization = visualizeExpressionData('zone2_muscle_cell_expression_data.csv', 'zone2_muscle_cell_expression_visualization.png')
zone3Visualization = visualizeExpressionData('zone3_muscle_cell_expression_data.csv', 'zone3_muscle_cell_expression_visualization.png')
zone2MuscleCellsNormalized = normalizeExpressionData('zone2_muscle_cell_expression_data.csv', 'zone2_muscle_cell_normalized_expression_data.csv')
zone3MuscleCellsNormalized = normalizeExpressionData('zone3_muscle_cell_expression_data.csv', 'zone3_muscle_cell_normalized_expression_data.csv')
zone2NormalizedVisualization = visualizeExpressionData('zone2_muscle_cell_normalized_expression_data.csv', 'zone2_muscle_cell_normalized_expression_visualization.png')
zone3NormalizedVisualization = visualizeExpressionData('zone3_muscle_cell_normalized_expression_data.csv', 'zone3_muscle_cell_normalized_expression_visualization.png')
