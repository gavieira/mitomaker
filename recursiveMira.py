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

def recursiveMira(processName = 'teste', inputFile = 'teste.input', kmers = [16,31,43,53,71], processorsToUse = 4, 
                  miraFolder = 'installed', sizeToLook = 16000, refSeqFile = None, cutoffValue=0.125, 
                  blasteVal = 0.0001, blastHitSizePercentage = 0.60, miraGenome = False, miraTechnology = 'solexa',
                  buildCloroplast = False, skipTrnaScan = False, circularSize = 50, circularOffSet = 220, ignoreFirstBuildChecks = False,
                  cutoffEquality = 0.60, organismType = 2, blastFolder = 'installed', noExtension = False, coveCutOff = 8, buildBacteria = False,
				  buildArchea = False):
	'''
	Run Mira N times looking for target DNA
	'''
	pathToMira = miraFolder
	
	print 'Starting recursiveMira phase...'

	if not 'default' in kmers:
		kmers.insert(0,'default') #if default was not inputed by user, do mira assembly with default values first, because they are dependant on technology
	
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
		############ Run Mira!!! ############
		#####################################

		#create mira manifest file
		miraInputFile = pathToWork + processName + '-denovo.manifest'

		with open(inputFile, 'r') as soapInputFile: #soapInputFile, just for easy copy pasting from recursiveSOAP, should alterations be needed
			with open(miraInputFile, 'w') as manifestFile:
				
				manifestFile.write('project = ' + processName + '-denovo\n')
				if miraGenome == False: #run in est mode (transcriptome)
					manifestFile.write('job = est, denovo, accurate\n')
				else: #run in genome mode
					manifestFile.write('job = genome, denovo, accurate\n')

				if currentKmer == 'default': #run first with default value
					manifestFile.write('parameters = -NW:cmrnl=no:cac=no:cnfs=warn -GE:not=' + str(processorsToUse) + '\n\n')
					print ''
					print "Running Mira first with it's default values, which are dependant on technology."
					print ''
				else:
					manifestFile.write('parameters = -NW:cmrnl=no:cac=no:cnfs=warn -GE:not=' + str(processorsToUse) + ' -SK:bph=' + str(currentKmer) + '\n\n')

				listOfInputs = {}
				listOfTechnologies = {}
				listOfOrientations = {}
				listOfInserts = {}
				listOfQualities = {}
				n = 0
				#reading the soap format input file
				for soapLine in soapInputFile:
					if '[LIB]' in soapLine: #start of a new library
						n += 1
						listOfInputs[n] = []
						listOfTechnologies[n] = miraTechnology
						listOfOrientations[n] = 'autopairing'
						listOfInserts[n] = None
						listOfQualities[n] = None
					if 'q' == soapLine[0] or 'f' == soapLine[0]: #append the data files
						listOfInputs[n].append(soapLine.replace('\n','').split('=')[-1]) #grab only the file part of this line
					elif 'technology' in soapLine:
						listOfTechnologies[n] = soapLine.replace('\n','').split('=')[-1]
					elif 'orientation' in soapLine:
						listOfOrientations[n] = soapLine.replace('\n','').split('=')[-1]
					elif 'avg_ins' in soapLine:
						avgInsert = int(soapLine.replace('\n','').split('=')[-1])
						insertLowEnd = avgInsert / 2
						insertHighEnd = avgInsert * 2
						listOfInserts[n] = str(insertLowEnd) + ' ' + str(insertHighEnd) + ' autorefine'
					elif 'default_qual' in soapLine:
						listOfQualities[n] = soapLine.replace('\n','').split('=')[-1]

				#time to grab the readgroups from the dict
				for readGroup in listOfInputs:
					manifestFile.write('readgroup\n')
					manifestFile.write('data =')
					for readFile in listOfInputs[readGroup]:
						manifestFile.write(' ' + readFile)
					manifestFile.write('\n')
					if len(listOfInputs[readGroup]) > 1:
						if listOfOrientations[readGroup] == 'autopairing' or listOfInserts[readGroup] == None:
							manifestFile.write('autopairing\n')
						else:
							manifestFile.write('template_size = ' + listOfInserts[readGroup] + '\n')
							manifestFile.write('segment_placement = ' + listOfOrientations[readGroup] + '\n')
					manifestFile.write('technology = ' + listOfTechnologies[readGroup] + '\n')
					if listOfQualities[readGroup] != None:
						manifestFile.write('default_qual = ' + listOfQualities[readGroup] + '\n')
					manifestFile.write('strain = ' + processName + '_' + str(readGroup) + '\n\n')
					
			
		#create Mira logfile:
		with open(pathToWork + 'mira_' + str(currentKmer) + 'mer.log','w') as miraLogFile:
			#soapTransLogFile = open(pathToWork + 'soap_' + str(currentKmer) + 'mer.log','w')
			#Run Mira (in EST or genome) and wait for it to finish.
			print 'Running Mira...'
			try:
				if pathToMira == 'installed':
					command = 'mira %s-denovo.manifest' % processName
				else:
					command = '%sbin/mira %s-denovo.manifest' %(pathToMira, processName)
				args = shlex.split(command)
				miraRun = Popen(args, cwd=pathToWork, stdout=miraLogFile)
				miraRun.wait()
				#check MIRA output to see if reference sequence was built
				numberOfReadGroups = n
				soapWasChecked = FirstBuildChecker.checkSoapOutput(processName, pathToWork, sizeToLook, refSeqFile, cutoffValue, blasteVal,
																   blastHitSizePercentage, False, numberOfReadGroups, buildCloroplast, skipTrnaScan, circularSize,
																   circularOffSet, cutoffEquality, organismType, blastFolder, noExtension,
																   ignoreFirstBuildChecks, coveCutOff, buildBacteria, buildArchea)
			except KeyboardInterrupt:
				return False
			except:
				print ''
				print "An error occured while running MIRA4 for DeNovo assembly. Check it's logs for more information."
				print 'Printing last 10 lines of log file and aborting...\n'
				with open(pathToWork + 'mira_' + str(currentKmer) + 'mer.log','r') as logFile:
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
