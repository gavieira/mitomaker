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

def recursiveSOAPdenovo(processName = 'teste', shortestContig = 100, inputFile = 'teste.input', kmers = [23,31,43,53,71], processorsToUse = 4, 
                  soapDeNovoFolder = 'installed', sizeToLook = 16000, refSeqFile = None, cutoffValue=0.125,
                  blasteVal = 0.0001, blastHitSizePercentage = 0.60, buildCloroplast = False, skipTrnaScan = False, circularSize = 50,
                  circularOffSet = 220, ignoreFirstBuildChecks = False, cutoffEquality = 0.60, organismType = 2, blastFolder = 'installed',
                  noExtension = False, coveCutOff = 8, buildBacteria = False, buildArchea = False):
	'''
	Run SOAPdenovo2 N times looking for a sequence that ressembles a reference one, or is close to a size we are looking for
	'''
	pathToSOAP = soapDeNovoFolder
	bestBuild = None
	
	print 'Starting recursive SOAPdenovo phase...'
	
	for currentKmer in kmers: #loops through different k-mers trying to 
		#create folders for the different kmers if they do not already exist
		print('========STARTING K-MER %s ==========' % currentKmer)
		print 'Creating folder for kmer = ' + str(currentKmer) + '. \nLog files will be saved there.'
		pathToWork = 'kmer_' + str(currentKmer) + '/'
		if not os.path.exists(pathToWork): os.makedirs(pathToWork)
			
		#copy input file to secondary folder to keep it organized
		destFile = pathToWork + inputFile
		shutil.copyfile(inputFile, destFile)

		#####################################
		###### Run SOAPdenovo-Trans!!! ######
		#####################################

		if int(currentKmer) <= 63: #check if I need to run Trans-63mer or 127mer
			soapVerToRun = '63'
		else:
			soapVerToRun = '127'
		print 'Running SOAPdenovo with ' + soapVerToRun + 'mer version.'
			
		#create SOAPdenovo-Trans logfile:
		try:
			with open(pathToWork + 'soap_' + str(currentKmer) + 'mer.log','w') as soapDeNovoLogFile:
				'''
				Run SOAPdenovo and wait for it to finish.
				'''
				command = '%sSOAPdenovo-%smer all -s %s -K %s -o %s -p %s -L %s' %(pathToSOAP, soapVerToRun, inputFile, 
							currentKmer, processName, processorsToUse, shortestContig)
				args = shlex.split(command)
				soapDeNovo = Popen(args, cwd=pathToWork, stdout=soapDeNovoLogFile, stderr=soapDeNovoLogFile)
				soapDeNovo.wait()
				#check SOAP output to see if reference sequence was built
				soapWasChecked = FirstBuildChecker.checkSoapOutput(processName, pathToWork, sizeToLook, refSeqFile, cutoffValue, blasteVal,
																   blastHitSizePercentage, True, 1, buildCloroplast, skipTrnaScan, circularSize,
																   circularOffSet, cutoffEquality, organismType, blastFolder, noExtension,
																   ignoreFirstBuildChecks, coveCutOff, buildBacteria, buildArchea)
		except KeyboardInterrupt:
			return False
		except:
			print ''
			print "An error occured while running SOAPdenovo for DeNovo assembly. Check it's logs for more information."
			print 'Printing last 10 lines of log file and aborting...\n'
			with open(pathToWork + 'soap_' + str(currentKmer) + 'mer.log','r') as logFile:
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
