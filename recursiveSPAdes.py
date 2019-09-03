#!/usr/bin/env python
#Version: 1.0
#Author: Alex Schomaker - alexschomaker@ufrj.br
#LAMPADA - IBQM - UFRJ

'''
Copyright (c) 2014 Alex Schomaker Bastos - LAMPADA/UFRJ

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

from subprocess import Popen
import shlex, os, shutil, FirstBuildChecker

#this function changes an orientation formatted to MIRA into SPAdes format
def spadesOrientation(orientationGiven):
	'''
	orientation (MIRA only): what is the orientation of paired
	reads? ---> <--- / ---> ---> / <--- <--- / <--- ---> (usable with MIRA
	mode. If ommited MIRA will automatically guess this information).
	'''
	dictOfOrientations = { '---> <---':'"fr"', '---> --->':'"ff"', 
						   '<--- <---':'"ff"', '<--- --->':'"rf"'}
	return dictOfOrientations[orientationGiven]

def recursiveSpades(processName = 'teste', inputFile = 'teste.input', kmers = [23,31,43,53,71], processorsToUse = 4, 
                  spadesFolder = 'installed', sizeToLook = 16000, refSeqFile = None, cutoffValue=0.125, 
                  blasteVal = 0.0001, blastHitSizePercentage = 0.60, miraTechnology = 'solexa',
                  buildCloroplast = False, skipTrnaScan = False, circularSize = 50, circularOffSet = 220, ignoreFirstBuildChecks = False,
                  cutoffEquality = 0.60, organismType = 2, blastFolder = 'installed', noExtension = False, coveCutOff = 7, buildBacteria = False,
				  buildArchea = False):
	'''
	Run SPAdes N times looking for target DNA
	'''
	pathToSpades = spadesFolder
	
	print 'Starting recursiveSPAdes phase...'

	bestBuild = None
	
	for currentKmer in kmers: #loops through different k-mers trying to assemble a target DNA
		#create folders for the different kmers if they do not already exist
		print('========== STARTING K-MER %s ==========' % currentKmer)
		print 'Creating folder for kmer = ' + str(currentKmer) + '. \nLog files will be saved there.'
		pathToWork = 'kmer_' + str(currentKmer) + '/'
		if not os.path.exists(pathToWork): os.makedirs(pathToWork)
			
		#copy input file to secondary folder to keep it organized
		destFile = pathToWork + inputFile
		shutil.copyfile(inputFile, destFile)

		#####################################
		########### Run SPAdes!!! ###########
		#####################################

		#create spades yaml file
		spadesInputFile = pathToWork + processName + '-denovo.yaml'

		with open(inputFile, 'r') as soapInputFile: #soapInputFile, just for easy copy pasting from recursiveSOAP, should alterations be needed
			with open(spadesInputFile, 'w') as manifestFile:
				
				manifestFile.write('[\n')

				listOfInputs = {}
				listOfTechnologies = {}
				listOfOrientations = {}
				listOfInputsOrder = {}
				n = 0
				#reading the soap format input file
				for soapLine in soapInputFile:
					if '[LIB]' in soapLine: #start of a new library
						n += 1
						listOfInputs[n] = []
						listOfInputsOrder[n] = []
						listOfTechnologies[n] = miraTechnology
						listOfOrientations[n] = None
					if 'q' == soapLine[0] or 'f' == soapLine[0]: #append the data files
						listOfInputs[n].append(soapLine.replace('\n','').split('=')[-1]) #grab only the file part of this line
						listOfInputsOrder[n].append(soapLine.replace('\n','').split('=')[0])
					elif 'technology' in soapLine:
						listOfTechnologies[n] = soapLine.replace('\n','').split('=')[-1]
					elif 'orientation' in soapLine:
						listOfOrientations[n] = soapLine.replace('\n','').split('=')[-1]

				#time to grab the readgroups from the dict
				for readGroup in listOfInputs:
					manifestFile.write(' {\n') #new dataset
					
					if len(listOfInputs[readGroup]) > 1:
						if listOfOrientations[readGroup] != None:
							manifestFile.write('  orientation: %s,\n' % spadesOrientation(listOfOrientations[readGroup]))
						else:
							manifestFile.write('  orientation: "fr",\n')
						manifestFile.write('  type: "paired-end",\n')
						#write the reads paths
						numberOfReadFile = 0
						for readFile in listOfInputs[readGroup]:
							numberOfReadFile += 1
							if listOfInputsOrder[readGroup][numberOfReadFile - 1] == 'q1':
								manifestFile.write('  left reads: [%s],\n' % readFile)
							else:
								manifestFile.write('  right reads: [%s],\n' % readFile)
					else:
						if listOfTechnologies[readGroup] == 'pacbio':
							manifestFile.write('  type: "pacbio",\n')
						else:
							manifestFile.write('  type: "single",\n')
					
							for readFile in listOfInputs[readGroup]:
								manifestFile.write('  single reads: [%s],\n' % readFile)
					
					manifestFile.write(' },\n') #end read group
				manifestFile.write(']') #end yaml file
					
			
		#create Mira logfile:
		with open(pathToWork + 'spades_' + str(currentKmer) + 'mer.log','w') as spadesLogFile:
			#Run SPAdes and wait for it to finish.
			print 'Running SPAdes...'
			try:
				if pathToSpades == 'installed':
					command = 'spades.py -k %s --only-assembler --careful --dataset %s-denovo.yaml -o %s -t %s' % (currentKmer, processName,
					          processName, processorsToUse)
				else:
					command = '%sbin/spades.py -k %s --only-assembler --careful --dataset %s-denovo.yaml -o %s -t %s' % (pathToSpades,
					          currentKmer, processName, processName, processorsToUse)
				args = shlex.split(command)
				spadesRun = Popen(args, cwd=pathToWork, stdout=spadesLogFile)
				spadesRun.wait()
				#check SPAdes output to see if reference sequence was built
				numberOfReadGroups = n
				soapWasChecked = FirstBuildChecker.checkSoapOutput(processName, pathToWork, sizeToLook, refSeqFile, cutoffValue, blasteVal,
																   blastHitSizePercentage, 'Spades', numberOfReadGroups, buildCloroplast, skipTrnaScan, circularSize,
																   circularOffSet, cutoffEquality, organismType, blastFolder, noExtension,
																   ignoreFirstBuildChecks, coveCutOff, buildBacteria, buildArchea)
			except KeyboardInterrupt:
				return False
			except:
				print ''
				print "An error occured while running SPAdes for DeNovo assembly. Check it's logs for more information."
				print 'Printing last 10 lines of log file and aborting...\n'
				with open(pathToWork + 'spades_' + str(currentKmer) + 'mer.log','r') as logFile:
					content = logFile.readlines()
					for n in xrange(-1,-11,-1):
						print content[n]
				print ''
				raise
				if bestBuild == None:
					return False
				else:
					print 'WARNING: An error occured, but we are procceding with the best build so far anyway...'
					return (bestBuild, str(bestKmer))
			
		#check if target dna was successfully built
		if ignoreFirstBuildChecks == True:
			print 'Ignoring tRNAscan, circular and genomic checks...'
			print ''
			if bestBuild == None:
				bestBuild = soapWasChecked
				bestKmer = currentKmer
			elif len(soapWasChecked.refSeq) > len(bestBuild.refSeq): #lastly, check if it's bigger
				bestBuild = soapWasChecked
				bestKmer = currentKmer
		if soapWasChecked != True and soapWasChecked != False and ignoreFirstBuildChecks == False:
			if bestBuild == None:
				bestBuild = soapWasChecked
				bestKmer = currentKmer
			else:
				'''
				Down here doing comparison between best and current build.
				Firstly check if we have all features built and the best assembly so far doesn't.
				Afterwards check the number of complete genes found. Then tRNAs (stored in len of Assembly).
				If still no better build found, check for which build has the most features built and less splits.
				Lastly check circularization and size of sequence.
				'''
				thisBuildPresentFeatures = soapWasChecked.checkFeatures[0]
				thisBuildImportantFeatures = soapWasChecked.checkFeatures[1]
				thisBuildSplits = len(soapWasChecked.checkFeatures[2])
				thisBuildCompleteGenes = len(soapWasChecked.checkFeatures[3])
				thisBuildValidContigs = soapWasChecked.validContigs
				bestBuildPresentFeatures = bestBuild.checkFeatures[0]
				bestBuildImportantFeatures = bestBuild.checkFeatures[1]
				bestBuildSplits = len(bestBuild.checkFeatures[2])
				bestBuildCompleteGenes = len(bestBuild.checkFeatures[3])
				bestBuildValidContigs = bestBuild.validContigs
				if thisBuildPresentFeatures >= thisBuildImportantFeatures and bestBuildPresentFeatures < bestBuildImportantFeatures and \
				   len(soapWasChecked) >= len(bestBuild) and thisBuildSplits <= bestBuildSplits:
					bestBuild = soapWasChecked
					bestKmer = currentKmer
				elif thisBuildCompleteGenes > bestBuildCompleteGenes and len(soapWasChecked) >= len(bestBuild): #check if it has more completed features
					bestBuild = soapWasChecked
					bestKmer = currentKmer
				elif len(soapWasChecked) > len(bestBuild) and thisBuildCompleteGenes >= bestBuildCompleteGenes: #check if it has more tRNAs
					bestBuild = soapWasChecked
					bestKmer = currentKmer
				elif thisBuildPresentFeatures > bestBuildPresentFeatures and thisBuildSplits <= bestBuildSplits and \
				     len(soapWasChecked) >= len(bestBuild): #compare splits
					bestBuild = soapWasChecked
					bestKmer = currentKmer
				elif thisBuildValidContigs < bestBuildValidContigs and thisBuildCompleteGenes >= bestBuildCompleteGenes and \
				     len(soapWasChecked) >= len(bestBuild) and thisBuildPresentFeatures >= bestBuildPresentFeatures and \
				     thisBuildSplits <= bestBuildSplits:
					bestBuild = soapWasChecked
					bestKmer = currentKmer
				elif soapWasChecked.isCircular() == True and bestBuild.isCircular() == False and thisBuildValidContigs <= bestBuildValidContigs and \
				     thisBuildCompleteGenes >= bestBuildCompleteGenes and len(soapWasChecked) >= len(bestBuild) and thisBuildPresentFeatures >= bestBuildPresentFeatures and \
				     thisBuildSplits <= bestBuildSplits:
					bestBuild = soapWasChecked
					bestKmer = currentKmer
				elif len(soapWasChecked.refSeq.seq) - soapWasChecked.refSeq.seq.lower().count('n') > \
				     len(bestBuild.refSeq.seq) - bestBuild.refSeq.seq.lower().count('n') and thisBuildValidContigs <= bestBuildValidContigs \
				     and thisBuildCompleteGenes >= bestBuildCompleteGenes and len(soapWasChecked) >= len(bestBuild) and \
				     thisBuildPresentFeatures >= bestBuildPresentFeatures and thisBuildSplits <= bestBuildSplits: #lastly, check if it's bigger
					bestBuild = soapWasChecked
					bestKmer = currentKmer
			if currentKmer == kmers[-1]:
				print 'Best kmer chosen = %s' % bestKmer
				return (bestBuild, str(bestKmer))
		elif soapWasChecked == True:
		#procceed to next step...
			print 'Best kmer chosen = %s' % currentKmer
			return (True, str(currentKmer))
		elif currentKmer == kmers[-1] and bestBuild == None: #if we already tried all k-mers, just give up
			return False
		elif currentKmer == kmers[-1]:
			print 'Best kmer chosen = %s' % bestKmer
			return (True, str(bestKmer))
