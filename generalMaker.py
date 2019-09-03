#!/usr/bin/python
#Version: 1.14
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

import recursiveSOAP, recursiveSOAPdenovo, recursiveMira, miraMapping, mitoBimWrapper, \
	circularizationCheck, tRNAscanChecker, geneChecker, genbankOutput, recursiveSPAdes
import argparse, os, shlex, shutil, sys
from tRNAscanChecker import tRNAconvert, prettyRNAName
from geneChecker import createImageOfAnnotation
from subprocess import Popen
from Bio import SeqIO, SeqFeature, SeqUtils
from Bio.Alphabet import generic_dna, generic_protein

'''
Main script for handling organellar DNA building pipeline.
'''

if __name__ == "__main__":
	'''
	Parse the arguments and start the process!
	'''
	parser = argparse.ArgumentParser(description='Try to build a target DNA with input reads.')
	parser.add_argument('-j', '--jobname', help = 'Job name to be used throughout the project', required=True, dest='processName')
	parser.add_argument('-i', '--input', help='Input file of reads, in SOAPdenovo format', required=True, dest='inputFile')
	parser.add_argument('-l', '--length', help='Shortest contig length to be used in scaffolding. Default = 100', type=int,
						default=100, dest='shortestContig')
	parser.add_argument('-k', '--kmers', help='The number of kmers you want the program to try (default = 23,31)\n\
						   When running mitoMaker default values are 23,31,43,53',\
						default='23,31', dest='kmers')
	parser.add_argument('-p', '--processors', help='Number of threads SOAPdenovo-Trans and Mira will use at most.', type=int,
						default=4, dest='processorsToUse')
	parser.add_argument('-r', '--refseq', help='What reference sequence should we look for, in fasta or genbank? Example: mitochondrial DNA of related species',
						default=None, dest='refSeqFile')
	parser.add_argument('-op', '--optimum', help='What optimum length of sequence should we look for? Ex: 16kb for mtDNA.\nIf refSeq is present, look for a sequence with its size - a cutoff value',
						default=-1, type=int, dest='sizeToLook')
	parser.add_argument('-c', '--mincutoff', help='Minimum size to consider. Default = -1\nUse -1 to automatically set this parameter.', type=int,
						default=-1, dest='minCutOffValue')
	parser.add_argument('-mc', '--maxcutoff', help='Maximum size to consider. Default = -1\nUse -1 to automatically set this parameter.', type=int,
						default=-1, dest='maxCutOffValue')
	parser.add_argument('-e', '--blaste', help='e-Value for blastn program. Default = 5.0', type=float,
						default=5.0, dest='blasteVal')
	parser.add_argument('--miratech', help='Read technology for MIRA assembly, ex: solexa\nNeeded only when using Mira 3.4',
						default='solexa', dest = 'miraTechnology')
	parser.add_argument('-mti', '--mitobimiter', help='Number of iterations to run Mitobim.pl. Default = 10',
						default=10, type=int, dest='mitobimIterations')
	parser.add_argument('--blastpercent', help='Percentage of covered span in blast best hit to be considered good. Default = 0.60',
						default=0.60, type=float, dest='blastHitSizePercentage')
	parser.add_argument('--mapmira3', help='Do mapping with mira3.4. Default = False',
						default=True, dest='useNewMira', action='store_false')
	parser.add_argument('--keepfolders', help='Keep temporary folders once assembly is done. Default = False',
						default=False, dest='keepTmpFolders', action='store_true')
	parser.add_argument('--circularsize', help='Size to consider when checking for circularization. Default = 45',
						default=45, type=int, dest='circularSize')
	parser.add_argument('--circularoffset', help='Offset from start and finish to consider when looking for circularization. Default = 200',
						default=200, type=int, dest='circularOffSet')
	parser.add_argument('--autopairend', help='Have mira4.0, during mapping phase, auto-discover the pairing information? Informational only, paired information is not needed for mapping assembly.',
						default=False, dest='pairedEnd', action='store_true')
	parser.add_argument('--recursivemira', help='Use Mira in recursive mode, by default it will run in genome mode.\nUse --miraest to run it in est mode instead of genome.',
						default=False, dest='recursiveMira', action='store_true')
	parser.add_argument('--soaptrans', help='Use SOAP in Transcriptome mode, by default it will run in DeNovo mode.\n\
						When using mitoMaker.py it will run in Transcriptome mode instead.',\
						default=False, dest='soapTrans', action='store_true')
	parser.add_argument('--forcedenovo', help='Force SOAP usage in DeNovo mode.',
			    default=False, dest='forceDeNovo', action='store_true')
	parser.add_argument('--spades', help='Use SPAdes to assemble genome.',
			    default=False, dest='useSpades', action='store_true')
	parser.add_argument('--miraest', help='Use Mira in recursive mode, running in est mode.',
						default=True, dest='miraGenome', action='store_false')
	parser.add_argument('--chloroplast', help='Are you building a cloroplast? Important for tRNA ang gene checker.',
						default=False, dest='buildCloroplast', action='store_true')
	parser.add_argument('--bacteria', help='Are you building a non organellar bacterial DNA? Important for tRNA checker.',
						default=False, dest='buildBacteria', action='store_true')
	parser.add_argument('--archea', help='Are you building a non organellar archea DNA? Important for tRNA checker.',
						default=False, dest='buildArchea', action='store_true')
	parser.add_argument('--nocopykmers', help='Should mira mapping work with a different  k-mer as the result from SOAP or Mira?',
						default=True, dest='copyKmers', action='store_false')
	parser.add_argument('-notrna', '--skiptrna', help='Should mitomaker ignore tRNAscan-SE checks? Default = False',
						default=False, dest='skipTrnaScan', action='store_true')
	parser.add_argument('--relaxed', help='Ignore tRNAscanCheck, genomic checks and Circularization Check during initial De Novo assembly? Default = False',
						default=False, dest='ignoreFirstBuildChecks', action='store_true')
	parser.add_argument('-ceq', '--cutoffequality', help='Percentage of identity for gene checking to consider that the gene is present. Default = 0.55', type=float,
						default=0.55, dest='cutoffEquality')
	parser.add_argument('-cove', '--covecutoff', help='Cove cutoff for tRNAscan-SE. Default = 7', type=int,
						default=7, dest='coveCutOff')
	parser.add_argument('-o', '--organism', help="What should the genome checking and annotation consider as genetic code type. NCBI's table (integer):\n\
					1. The Standard Code\
					2. The Vertebrate Mitochondrial Code\
					3. The Yeast Mitochondrial Code\
					4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code\
					5. The Invertebrate Mitochondrial Code\
					6. The Ciliate, Dasycladacean and Hexamita Nuclear Code\
					9. The Echinoderm and Flatworm Mitochondrial Code\
					10. The Euplotid Nuclear Code\
					11. The Bacterial, Archaeal and Plant Plastid Code\
					12. The Alternative Yeast Nuclear Code\
					13. The Ascidian Mitochondrial Code\
					14. The Alternative Flatworm Mitochondrial Code\
					16. Chlorophycean Mitochondrial Code\
					21. Trematode Mitochondrial Code\
					22. Scenedesmus obliquus Mitochondrial Code\
					23. Thraustochytrium Mitochondrial Code\
					24. Pterobranchia Mitochondrial Code\
					25. Candidate Division SR1 and Gracilibacteria Code",
						type=int, default=1, dest='organismType')
	parser.add_argument('--skipmitobim', help="Don't run MITObim after mapping assembly? Default = False",
						default=False, dest='skipMitobim', action='store_true')
	parser.add_argument('--skipdenovo', help="Skip DeNovo phase? Default = False \nYou need to input a fasta or genbank file \
						containing a sequence to start the process", \
						default=False, dest='skipFirstStep')
	parser.add_argument('-d', '--dloopsize', help='Expected size of control region for automated annotation. Ignored if --chloroplast flag is on\nor if = -1. Default = 900', type=int,
						default=900, dest='dLoopSize')
	parser.add_argument('--edge', help='Size of edge at the result sequence to be blasted against when looking for missing parts. Default = 50', type=int,
						default=50, dest='edgesToLook')
	parser.add_argument('--noextension', help="Don't try to extend De-Novo assembly? Default = False",
						default=False, dest='noExtension', action='store_true')
	#parser.add_argument('--version', help="Version=1.14", default=False, dest='versionCheck', action='store_true')
	args = parser.parse_args()

	blasteVal = args.blasteVal
	usingOwnGenBankReference = False
	
	print 'Command line: %s' % ' '.join(sys.argv)
	print 'Now running generalMaker.py ...'

	'''
	Read config file and import information.
	'''
	module_dir = os.path.dirname(__file__)
	module_dir = os.path.abspath(module_dir)
	cfg_full_path = os.path.join(module_dir, 'generalMaker.config')

	with open(cfg_full_path,'r') as configFile:
		for line in configFile:
			if '#' != line[0] and line != '\n':
				configPart = line.lower().replace('\n','').replace(' ','').split('=')[0]
				if configPart == 'mitobimfolder':
					mitobimFolder = line.replace('\n','').replace(' ','').split('=')[-1]
				elif configPart == 'mitobim1.7folder':
					newMitobimFolder = line.replace('\n','').replace(' ','').split('=')[-1]
				elif configPart == 'mira4folder':
					pathToNewMira = line.replace('\n','').replace(' ','').split('=')[-1]
				elif configPart == 'mira3folder':
					pathToOldMira = line.replace('\n','').replace(' ','').split('=')[-1]
				elif configPart == 'soaptransfolder':
					soapTransFolder = line.replace('\n','').replace(' ','').split('=')[-1]
				elif configPart == 'soapdenovofolder':
					soapDeNovoFolder = line.replace('\n','').replace(' ','').split('=')[-1]
				elif configPart == 'spadesfolder':
					spadesFolder = line.replace('\n','').replace(' ','').split('=')[-1]
				elif configPart == 'blastfolder':
					blastFolder = line.replace('\n','').replace(' ','').split('=')[-1]
					
	#if config file has 'default' in the folder field, use the default program folders given with the script
	if mitobimFolder.lower() == 'default':
		mitobimFolder = os.path.join(module_dir, 'mitobim1.6/')
		print 'WARNING: MITObim1.6 is still set as default folder. Change it in the config file if you encounter problems running this script.\n'
		
	if newMitobimFolder.lower() == 'default':
		newMitobimFolder = os.path.join(module_dir, 'mitobim1.7/')
		print 'WARNING: MITObim1.7 is still set as default folder. Change it in the config file if you encounter problems running this script.\n'
		
	if pathToNewMira.lower() == 'default':
		pathToNewMira = os.path.join(module_dir, 'mira4/')
		print 'WARNING: MIRA4.0 is still set as default folder. Change it in the config file if you encounter problems running this script.\n'
		
	if pathToOldMira.lower() == 'default':
		pathToOldMira = os.path.join(module_dir, 'mira3/')
		print 'WARNING: MIRA3.4 is still set as default folder. Change it in the config file if you encounter problems running this script.\n'
		
	if soapTransFolder.lower() == 'default':
		soapTransFolder = os.path.join(module_dir, 'soapdenovo-trans/')
		print 'WARNING: SOAPdenovo-Trans is still set as default folder. Change it in the config file if you encounter problems running this script.\n'

	if soapDeNovoFolder.lower() == 'default':
		soapDeNovoFolder = os.path.join(module_dir, 'soapdenovo/')
		print 'WARNING: SOAPdenovo is still set as default folder. Change it in the config file if you encounter problems running this script.\n'
		
	if spadesFolder.lower() == 'default':
		spadesFolder = os.path.join(module_dir, 'spades/')
		print 'WARNING: SPAdes is still set as default folder. Change it in the config file if you encounter problems running this script.\n'

	if blastFolder.lower() == 'default':
		blastFolder = os.path.join(module_dir, 'blast/')
		print 'WARNING: Blast is still set as default folder. Change it in the config file if you encounter problems running this script.\n'

	print 'Program folders:'
	print 'MITObim = %s' % mitobimFolder
	print 'MITObim1.7 = %s' % newMitobimFolder
	print 'MIRA4 folder = %s' % pathToNewMira
	print 'MIRA3 folder = %s' % pathToOldMira
	print 'SOAPtrans folder = %s' % soapTransFolder
	print 'SOAPdenovo folder = %s' % soapDeNovoFolder
	print 'SPAdes folder = %s' % spadesFolder
	print 'Blast folder = %s' % blastFolder
	print ''
	print ''

	if args.refSeqFile == None:
		print 'WARNING: You are not using a reference targeted assembly.\nYou should try and find a closely related reference in either fasta or genbank.'
		print 'Mitomaker works best when a reference is given. If you do not have one, check the references/ folder for a possible reference.'
		print('If you still have to use a non-reference assembly, remember to check the current flags:\n \
			--optimum ; -c ; -mc\nFor more information on each flag, call generalMaker.py with --help\n')
		print 'Current -o flag is set to: %s\nThis will be used for the auto-annotation and genomic check part, if a genbank reference is not given.\n\n' % args.organismType

	#just start the variables for future checking
	firstStep = None #recursiveSOAP or recursiveMIRA
	secondStep = None #miraMapping
	thirdStep = None #mitobim
	fourthStep = None #circularization check
	fifthStep = None #tRNAscan

	#make cutoffvalue tuple
	cutoffValue = (args.minCutOffValue, args.maxCutOffValue)

	#if flag --chloroplast is on and max cutoff value wasn't altered, increase it to cloroplast values
	#if values are -1, automatically find this information out.
	if args.buildCloroplast == True and args.refSeqFile == None:
		if cutoffValue[0] == -1:
			print 'WARNING: You did not specify a minimum target value. Auto determining one...\nIf you want to set it yourself, change the -c flag'
			cutoffValue = (250, cutoffValue[1])
		if cutoffValue[1] == -1:
			print 'WARNING: You did not specify a maximum target value. Auto determining one...\nIf you want to set it yourself, change the -mc flag'
			cutoffValue = (cutoffValue[0], 190000)
		if args.sizeToLook == -1:
			print 'WARNING: You did not specify an optimum target value. Auto determining one...\nIf you want to set it yourself, change the --optimum flag'
			args.sizeToLook = 160000
	elif args.buildCloroplast == False and args.refSeqFile == None:
		if cutoffValue[0] == -1:
			print 'WARNING: You did not specify a minimum target value. Auto determining one...\nIf you want to set it yourself, change the -c flag'
			cutoffValue = (250, cutoffValue[1])
		if cutoffValue[1] == -1:
			print 'WARNING: You did not specify a maximum target value. Auto determining one...\nIf you want to set it yourself, change the -mc flag'
			cutoffValue = (cutoffValue[0], 19000)
		if args.sizeToLook == -1:
			print 'WARNING: You did not specify an optimum target value. Auto determining one...\nIf you want to set it yourself, change the --optimum flag'
			args.sizeToLook = cutoffValue[1] * 0.85
	elif args.refSeqFile != None:
		if args.refSeqFile[-6:] == '.fasta':
			refSeq = SeqIO.read(args.refSeqFile, "fasta", generic_dna)
		else:
			refSeq = SeqIO.read(args.refSeqFile, "genbank", generic_dna)
		refSize = len(refSeq.seq)
		if cutoffValue[0] == -1:
			cutoffValue = (max(1,int(refSize * 0.01)), cutoffValue[1])
		if cutoffValue[1] == -1:
			cutoffValue = (cutoffValue[0], int(refSize * 1.175))
		if args.sizeToLook == -1:
			print 'WARNING: You did not specify an optimum target value. Auto determining one based on reference...\nIf you want to set it yourself, change the --optimum flag'
			args.sizeToLook = refSize

	print 'Minimum size to consider: %s' %cutoffValue[0]
	print 'Maximum size to consider: %s' %cutoffValue[1]
	print 'Optimum size to consider: %s' %args.sizeToLook
	print ''

	if args.buildCloroplast or args.buildBacteria or args.buildArchea:
		args.organismType = 11

	maxReadLen = 100
	with open(args.inputFile, 'r') as soapInputFile:
		for soapLine in soapInputFile:
			if 'max_rd_len' in soapLine:
				maxReadLen = int(soapLine.replace('\n','').split('=')[-1])
				break

	print 'Read length to be used in MITObim: %s\n\n' %maxReadLen

	#Let's call recursiveSOAP or recursiveMIRA (runs SOAP or MIRA multiple times trying to find a referenced DNA!)
	if args.skipFirstStep != False:
		print '--skipdenovo is pointing to a file, going to skip DeNovo step...'
		kmerToUse = args.kmers.lower().split(',')[0]
		print 'Using %s as target k-mer' % kmerToUse
		print ''
		firstStep = (None, kmerToUse)
	else:
		if args.useSpades == True:
			firstStep = recursiveSPAdes.recursiveSpades(processName = args.processName, 
						    inputFile = args.inputFile, kmers = args.kmers.lower().split(','), 
		                    processorsToUse = args.processorsToUse, spadesFolder=spadesFolder, sizeToLook = args.sizeToLook, 
		                    refSeqFile = args.refSeqFile, cutoffValue = cutoffValue, blasteVal = args.blasteVal,
		                    blastHitSizePercentage = args.blastHitSizePercentage,miraTechnology=args.miraTechnology,
		                    buildCloroplast = args.buildCloroplast, skipTrnaScan = args.skipTrnaScan, circularSize = args.circularSize,
		                    circularOffSet = args.circularOffSet, ignoreFirstBuildChecks = args.ignoreFirstBuildChecks,
		                    cutoffEquality = args.cutoffEquality, organismType = args.organismType,blastFolder = blastFolder, 
		                    noExtension = args.noExtension, coveCutOff = args.coveCutOff, buildBacteria = args.buildBacteria, 
		                    buildArchea = args.buildArchea)
		elif (args.soapTrans == False and args.recursiveMira == False) or args.forceDeNovo == True:
			firstStep = recursiveSOAPdenovo.recursiveSOAPdenovo(processName = args.processName, 
				    shortestContig = args.shortestContig, inputFile = args.inputFile, kmers = args.kmers.lower().split(','), 
		                    processorsToUse = args.processorsToUse, soapDeNovoFolder=soapDeNovoFolder, sizeToLook = args.sizeToLook, 
		                    refSeqFile = args.refSeqFile, cutoffValue = cutoffValue, blasteVal = args.blasteVal,
		                    blastHitSizePercentage = args.blastHitSizePercentage, buildCloroplast = args.buildCloroplast,
		                    skipTrnaScan = args.skipTrnaScan, circularSize = args.circularSize, circularOffSet = args.circularOffSet,
		                    ignoreFirstBuildChecks = args.ignoreFirstBuildChecks, cutoffEquality = args.cutoffEquality,
		                    organismType = args.organismType,blastFolder = blastFolder, noExtension = args.noExtension,
				    coveCutOff = args.coveCutOff, buildBacteria = args.buildBacteria, buildArchea = args.buildArchea)
		elif args.soapTrans == True and args.recursiveMira == False:
			firstStep = recursiveSOAP.recursiveSOAP(processName = args.processName, shortestContig = args.shortestContig,
			            inputFile = args.inputFile, kmers = args.kmers.lower().split(','), 
		                    processorsToUse = args.processorsToUse, soapTransFolder=soapTransFolder, sizeToLook = args.sizeToLook, 
		                    refSeqFile = args.refSeqFile, cutoffValue = cutoffValue, blasteVal = args.blasteVal,
		                    blastHitSizePercentage = args.blastHitSizePercentage, buildCloroplast = args.buildCloroplast,
		                    skipTrnaScan = args.skipTrnaScan, circularSize = args.circularSize, circularOffSet = args.circularOffSet,
		                    ignoreFirstBuildChecks = args.ignoreFirstBuildChecks, cutoffEquality = args.cutoffEquality,
		                    organismType = args.organismType,blastFolder = blastFolder, noExtension = args.noExtension,
				    coveCutOff = args.coveCutOff, buildBacteria = args.buildBacteria, buildArchea = args.buildArchea)
		elif args.recursiveMira == True:
			firstStep = recursiveMira.recursiveMira(processName = args.processName, inputFile = args.inputFile,
				    kmers = args.kmers.lower().split(','), processorsToUse = args.processorsToUse, miraFolder = pathToNewMira,
                                    sizeToLook = args.sizeToLook, refSeqFile = args.refSeqFile, cutoffValue = cutoffValue, 
                                    blasteVal = args.blasteVal, blastHitSizePercentage = args.blastHitSizePercentage, 
                                    miraGenome = args.miraGenome, miraTechnology = args.miraTechnology,
		                    buildCloroplast = args.buildCloroplast, skipTrnaScan = args.skipTrnaScan, circularSize = args.circularSize,
		                    circularOffSet = args.circularOffSet, ignoreFirstBuildChecks = args.ignoreFirstBuildChecks,
		                    cutoffEquality = args.cutoffEquality, organismType = args.organismType, blastFolder = blastFolder,
		                    noExtension = args.noExtension, coveCutOff = args.coveCutOff, buildBacteria = args.buildBacteria,
				    buildArchea = args.buildArchea)
	
	#time for second step, mira mapping
	if firstStep != False or args.skipFirstStep != False: #if soap/mira ran until the end
		print 'Starting second step (mira mapping)...'
		#make best_query.fasta hold the best sequence
		if firstStep[0] != True or args.skipFirstStep != False:
			with open('best_query.fasta', 'w') as bestQueryFile:
				#print firstStep[0].refSeq
				if args.skipFirstStep != False:
					if args.skipFirstStep.endswith('.fasta') or args.skipFirstStep.endswith('.fa'):
						seqForBestQuery = SeqIO.read(args.skipFirstStep, "fasta", generic_dna)
					else:
						seqForBestQuery = SeqIO.read(args.skipFirstStep, "genbank", generic_dna)
				else:
					seqForBestQuery = firstStep[0].refSeq
				SeqIO.write(seqForBestQuery, bestQueryFile, 'fasta')
		secondStep = miraMapping.miraMapping(processName = args.processName, processorsToUse = args.processorsToUse,
					inputFile = args.inputFile, miraTechnology = args.miraTechnology.lower(),
					useNewMira = args.useNewMira, pathToNewMira = pathToNewMira, pathToOldMira = pathToOldMira,
					pairedEnd = args.pairedEnd, copyKmers = args.copyKmers, lastKmer = firstStep[1])
		if secondStep == False:
			print 'miraMapping failed. Aborting. Check logs for info.'
		else:
			if args.skipMitobim == False:
			#procceed with MITObim if second step was successful...
				print 'Starting third step (mitobim)...'
				thirdStep = mitoBimWrapper.mitoBimWrapper(mitobimIterations = args.mitobimIterations, processName = args.processName,
							miraTechnology = args.miraTechnology.lower(), mitobimFolder = mitobimFolder, readLen = maxReadLen,
							newMira = args.useNewMira, newMitobimFolder = newMitobimFolder, pathToNewMira = pathToNewMira,
							pathToOldMira = pathToOldMira, kmerUsed = firstStep[1])
			if thirdStep == True:
				print ''
				print 'MITObim finished running.'
				print ''

			'''
			Do circularization check on MITObim or MIRA4 results...
			'''
			print 'Procceding to circularization check...'
			print ''

			#figuring out which result file to use for the sequence (as resultFile)...
			pathOfResult = None
			if thirdStep == None: #MITObim wasn't ran
				print "MITObim wasn't ran..."
				pathOfResult = 'mira_mapping/' + args.processName + '_assembly/' + args.processName + '_d_results/' + args.processName + '_out_AllStrains.unpadded.fasta'
				pathOfMafResult = 'mira_mapping/' + args.processName + '_assembly/' + args.processName + '_d_results/' + args.processName + '_out.maf'
				pathOfCafResult = 'mira_mapping/' + args.processName + '_assembly/' + args.processName + '_d_results/' + args.processName + '_out.caf'
			elif thirdStep == True: #MITObim was succesfully ran
				#let's check from top to bottom for iteration folders, when last iteration is found, grab result from that
				if args.useNewMira == False:
					for iteration in xrange(args.mitobimIterations, 0, -1):
						iterationFolder = 'iteration' + str(iteration) + '/'
						if os.path.exists(iterationFolder):
							pathOfResult = 'iteration' + str(iteration) + '/' + args.processName + '-ReferenceStrain_assembly/' + args.processName + '-ReferenceStrain_d_results/' + args.processName + '-ReferenceStrain_out.unpadded.fasta'
							pathOfMafResult = 'iteration' + str(iteration) + '/' + args.processName + '-ReferenceStrain_assembly/' + args.processName + '-ReferenceStrain_d_results/' + args.processName + '-ReferenceStrain_out.maf'
							pathOfCafResult = 'iteration' + str(iteration) + '/' + args.processName + '-ReferenceStrain_assembly/' + args.processName + '-ReferenceStrain_d_results/' + args.processName + '-ReferenceStrain_out.caf'
							print 'Using iteration ' + str(iteration) + ' for circularization checking.'
							break
				elif args.useNewMira == True:
					for iteration in xrange(args.mitobimIterations, 0, -1):
						iterationFolder = 'iteration' + str(iteration) + '/'
						if os.path.exists(iterationFolder):
							pathOfResult = 'iteration' + str(iteration) + '/' + args.processName + '_1-backbone_assembly/' + args.processName + '_1-backbone_d_results/' + args.processName + '_1-backbone_out_AllStrains.unpadded.fasta'
							pathOfMafResult = 'iteration' + str(iteration) + '/' + args.processName + '_1-backbone_assembly/' + args.processName + '_1-backbone_d_results/' + args.processName + '_1-backbone_out.maf'
							pathOfCafResult = 'iteration' + str(iteration) + '/' + args.processName + '_1-backbone_assembly/' + args.processName + '_1-backbone_d_results/' + args.processName + '_1-backbone_out.caf'
							print 'Using iteration ' + str(iteration) + ' for circularization checking.'
							break
							
			if pathOfResult is None: #if mitobim had a problem, use mira mapping as result
				print '#'*28
				print 'WARNING: There was a problem running MITObim. Going to use MIRA mapping assembly as result.'
				print '#'*28
				pathOfResult = 'mira_mapping/' + args.processName + '_assembly/' + args.processName + '_d_results/' + args.processName + '_out_AllStrains.unpadded.fasta'
				pathOfMafResult = 'mira_mapping/' + args.processName + '_assembly/' + args.processName + '_d_results/' + args.processName + '_out.maf'
				pathOfCafResult = 'mira_mapping/' + args.processName + '_assembly/' + args.processName + '_d_results/' + args.processName + '_out.caf'

			print ''
			print 'Checking results for circularization...'
			#circularizationcheck will return a tuple with (True, start, end)
			fourthStep = circularizationCheck.circularizationCheck(pathOfResult, args.circularSize, args.circularOffSet, blastFolder)
			print ''

			resultFile = args.processName + '.fasta'

			if fourthStep[0] == True:
				print 'Evidences of circularization were found!'
				print 'Sequence is going to be trimmed according to circularization position. \nMAF and CAF files are unaltered.'
				print ''
				with open(resultFile, "w") as outputResult: #create draft file to be checked and annotated
					finalResults = SeqIO.read(open(pathOfResult, 'rU'), "fasta", generic_dna)
					finalResults.seq = finalResults.seq.upper()
					count = SeqIO.write(finalResults[fourthStep[2]:], outputResult, "fasta") #trims according to circularization position
			else:
				print 'Evidences of circularization could not be found, but everyother step was successful.'
				print 'Check results from MIRA Mapping.'
				print ''
				with open(resultFile, "w") as outputResult: #create draft file to be checked and annotated
					finalResults = SeqIO.read(open(pathOfResult, 'rU'), "fasta", generic_dna)
					finalResults.seq = finalResults.seq.upper()
					count = SeqIO.write(finalResults, outputResult, "fasta") #no need to trim, since circularization wasn't found

			pathOfFinalResults = args.processName + '_Final_Results/'
			if not os.path.exists(pathOfFinalResults): os.makedirs(pathOfFinalResults)
			
			#copying results to the final results folder and creating draft file to be checked...
			print '#'*75
			print '## Creating target DNA draft file from mapping/MITObim results...'

			#creating some stat file:
			finalResults = SeqIO.read(open(resultFile, 'rU'), "fasta", generic_dna)
			finalStatsFile = open(pathOfFinalResults + args.processName + '.stats', 'w')

			finalStatsFile.write('Statistics for final sequence:\n\n')
			finalStatsFile.write('Length: ' + str(len(finalResults.seq)) + "\n")
			finalStatsFile.write('GC content: ' + ("{0:.2f}".format(SeqUtils.GC(finalResults.seq))) + '%\n')
			numberOfNs = finalResults.seq.lower().count('n')
			finalStatsFile.write("Length without Ns: " + str(len(finalResults.seq) - numberOfNs) + "\n")
			finalStatsFile.write("Number of Ns: " + str(numberOfNs) + "\n")
			if fourthStep[0] == True:
				finalStatsFile.write("Circularization: Yes\n")
			else:
				finalStatsFile.write("Circularization: No\n")
			finalStatsFile.write("K-mer used: " + str(firstStep[1]) + "\n")

			destFile = pathOfFinalResults + args.processName + '.unordered.fasta'
			shutil.copyfile(resultFile, destFile)
			
			destFile = pathOfFinalResults + args.processName + '.unordered.maf'
			shutil.copyfile(pathOfMafResult, destFile)
			
			destFile = pathOfFinalResults + args.processName + '.unordered.caf'
			shutil.copyfile(pathOfCafResult, destFile)
				
			print '## Final sequence saved to %s' % pathOfFinalResults 

			#from now on, just checking how the build went to output to user and then annotate
			print '## Now running tRNAscan-SE to check the final build...'
			
			if args.skipTrnaScan == True:
				print ''
				print '## --skiptrna is turned on, going to ignore this check and move on...'

			fifthStep = tRNAscanChecker.tRNAscanCheck(resultFile, fourthStep[0], args.skipTrnaScan, args.organismType, args.coveCutOff,
								args.buildBacteria, args.buildArchea) #returns a Assembly object 
												      #with statistics and alignment info
			if args.skipTrnaScan == False:
				print '## %s tRNAs were found.' % len(fifthStep)

			if args.ignoreFirstBuildChecks == False:
				print ''
				print '## Procceding to genomic check and annotation...'
				print ''

				#time to look for genomic features, searching according to the -o flag
				#or with the genbank reference that was given
				organismType = args.organismType
				module_dir = os.path.dirname(__file__)
				module_dir = os.path.abspath(module_dir)
				'''
					1. The Standard Code
					2. The Vertebrate Mitochondrial Code
					3. The Yeast Mitochondrial Code
					4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
					5. The Invertebrate Mitochondrial Code
					6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
					9. The Echinoderm and Flatworm Mitochondrial Code
					10. The Euplotid Nuclear Code
					11. The Bacterial, Archaeal and Plant Plastid Code
					12. The Alternative Yeast Nuclear Code
					13. The Ascidian Mitochondrial Code
					14. The Alternative Flatworm Mitochondrial Code
					16. Chlorophycean Mitochondrial Code
					21. Trematode Mitochondrial Code
					22. Scenedesmus obliquus Mitochondrial Code
					23. Thraustochytrium Mitochondrial Code
					24. Pterobranchia Mitochondrial Code
					25. Candidate Division SR1 and Gracilibacteria Code
				'''
				if args.buildCloroplast == True:
					refSeqFileForGenes = os.path.join(module_dir, 'references/cloroplast.gb')
					#args.organismType = 11
				elif args.organismType == 11: #plant plastid or bacterial or archea
					refSeqFileForGenes = os.path.join(module_dir, 'references/magnolia.gb')
					if args.buildBacteria == True:
						refSeqFileForGenes = os.path.join(module_dir, 'references/bacteria.gb') #some bacterial DNA
					elif args.buildArchea == True:
						refSeqFileForGenes = os.path.join(module_dir, 'references/archea.gb') #some archea DNA
				elif args.organismType == 3: #yeast
					refSeqFileForGenes = os.path.join(module_dir, 'references/yeast.gb')
				elif args.organismType == 5: #invertebrate, insect
					refSeqFileForGenes = os.path.join(module_dir, 'references/beetle.gb')
				elif args.organismType == 6: #ciliate
					refSeqFileForGenes = os.path.join(module_dir, 'references/paramecium.gb')
				elif args.organismType != 1: #human, 2(vertebrate), also the default option
					refSeqFileForGenes = os.path.join(module_dir, 'references/human.gb')

				#we don't need a sequence reference for this check
				if args.refSeqFile != None: #if user gave a reference file, let's consider its features and everything else
					if args.refSeqFile[-6:] != '.fasta':
						refSeqFileForGenes = args.refSeqFile
						usingOwnGenBankReference = True

				#add the genetic feature check to this Assembly object
				fifthStep.checkFeatures = geneChecker.geneCheck(refSeqFileForGenes, resultFile, args.cutoffEquality, usingOwnGenBankReference, blastFolder, organismType = args.organismType)
				presentFeatures = fifthStep.checkFeatures[0]
				importantFeatures = fifthStep.checkFeatures[1]
				numberOfSplits = len(fifthStep.checkFeatures[2])
				print '#'*75
				print ''
				print '#'*25
				print '## Annotating...'

				#Annotation and creation of genbank file down here:
				resultGbFile = pathOfFinalResults + args.processName + '.unordered.gb'
				listOfFeaturesToOutput = []
				listOfFoundTRNAs = []
				for foundFeature in presentFeatures:
					thisFeatureFound = presentFeatures[foundFeature][1]
					#comparing tRNAscan-SE results with this, in case tRNAscan-SE was run
					if "trn" in thisFeatureFound.seq2.lower():
						if args.skipTrnaScan == False:
							for tRNAFound in fifthStep.tRNAs:
							#down here we update the start and end positions of tRNAs found with Needle, with the
							#results outputted by tRNAScan-SE
							#tRNAconvert = guarantees all tRNA names are in tRNA-Phe format
								if 'trna-' + tRNAFound.tRNAtype.lower() == tRNAconvert(thisFeatureFound.seq2.lower()):
									#making sure the frame matches tRNAscan-SE results
									if tRNAFound.tRNAcoordinates[0] > tRNAFound.tRNAcoordinates[1]:
										thisFeatureFound.frame = -1
									#changing start and end according to tRNAscan-SE result
									thisFeatureFound.startBase = min(tRNAFound.tRNAcoordinates[0],
													tRNAFound.tRNAcoordinates[1])
									thisFeatureFound.endBase = max(tRNAFound.tRNAcoordinates[0],
													tRNAFound.tRNAcoordinates[1])
									break

						listOfFoundTRNAs.append(thisFeatureFound.seq2.lower())

					listOfFeaturesToOutput.append(thisFeatureFound)

				#if tRNAscan-SE was run, check the tRNAs it found and input them in the features to output list
				tRNAsFoundWithIntrons = [] #create empty list to hold problematic tRNAs and output them to the stats file
				if args.skipTrnaScan == False:
					for tRNAFound in fifthStep.tRNAs:
						tRNAName = 'trna-' + tRNAFound.tRNAtype.lower()
						if tRNAFound.tRNAintronBegin > 0:
							tRNAsFoundWithIntrons.append(prettyRNAName(tRNAName))
						if tRNAName not in tRNAconvert(listOfFoundTRNAs) and 'trna-sec' not in tRNAName and 'trna-sup' not in tRNAName:
							newTRNAStart = tRNAFound.tRNAcoordinates[0]
							newTRNAEnd = tRNAFound.tRNAcoordinates[1]
							newTRNALen = max(newTRNAStart, newTRNAEnd) - min(newTRNAStart, newTRNAEnd)
							#creating the new tRNA with the Alignment class and using prettyRNAName() to make sure it's
							#named as (for example) tRNA-Phe
							newTRNA = geneChecker.Alignment(prettyRNAName(tRNAName), prettyRNAName(tRNAName), newTRNALen)
							newTRNA.startBase = min(newTRNAStart, newTRNAEnd)
							newTRNA.endBase = max(newTRNAStart, newTRNAEnd)
							if newTRNAStart > newTRNAEnd:
								newTRNA.frame = -1
							else:
								newTRNA.frame = 1

							presentFeatures[prettyRNAName(tRNAName)] = (False, thisFeatureFound, False)

							listOfFeaturesToOutput.append(newTRNA)

				listOfFeaturesToOutput.sort()

				print '#'*25
				print ''
				print '#'*25
				print '#### FINAL RESULTS: ####'
				print '#'*25
				print 'Features found: ' + str(len(presentFeatures)) + ' / ' + str(len(importantFeatures))
				print 'Number of splits founds: ',numberOfSplits
				print 'Split or duplicated genes: ' + ', '.join(fifthStep.checkFeatures[2])
				print 'tRNAs found with introns: ' + str(len(tRNAsFoundWithIntrons))
				finalStatsFile.write('Features found: ' + str(len(presentFeatures)) + ' / ' + str(len(importantFeatures)) + '\n')
				if numberOfSplits > 0:
					finalStatsFile.write('Number of splits founds: ' + str(numberOfSplits) + '\n')
					finalStatsFile.write('Split or duplicated genes: ' + ', '.join(fifthStep.checkFeatures[2]) + '\n')
				if len(tRNAsFoundWithIntrons) > 0:
					finalStatsFile.write('tRNAs found with introns: ' + ', '.join(tRNAsFoundWithIntrons) + '\n')

				print ''
				finalStatsFile.write('\nFeatures not found:\n')
				#printing missing features to the stats file:
				for targetFeature in importantFeatures:
					#targetFeature = id of a feature that we searched for
					#presentFeatures = dict containing the alignment objects of all found features
					if targetFeature.lower().startswith('trn'):
						targetId = prettyRNAName(targetFeature)
					else:
						targetId = targetFeature
					if targetId not in presentFeatures and targetFeature not in presentFeatures:
						print '%s was not found.' % targetFeature
						finalStatsFile.write(targetFeature + '\n')
				print ''

				finalResults = genbankOutput.genbankOutput(resultGbFile, resultFile, listOfFeaturesToOutput,
                                                           buildCloroplast = args.buildCloroplast, dLoopSize = args.dLoopSize)
				finalResults.name = args.processName[0:10]
				finalResults.id = args.processName[0:10]
				
				resultSequinFile = resultGbFile.replace('.gb','.tbl')

				with open(resultGbFile, "w") as outputResult: #create the file!
					count = SeqIO.write(finalResults, outputResult, "genbank")
					createImageOfAnnotation(finalResults, resultGbFile.replace('.gb','.png'))
				with open(resultSequinFile, "w") as outputSeqIn:
					outputSeqIn.write('>Features ' + finalResults.name + '\n')
					for gbkFeature in finalResults.features:
						if gbkFeature.location.strand == 1 or gbkFeature.location.strand == None:
							outputSeqIn.write(str(gbkFeature.location.start + 1) + ' ' + str(gbkFeature.location.end)\
									 + ' ' + str(gbkFeature.type) + '\n\t\t')
						else:
							outputSeqIn.write(str(gbkFeature.location.end) + ' ' + str(gbkFeature.location.start + 1)\
									 + ' ' + str(gbkFeature.type) + '\n\t\t')
						for qualifier in gbkFeature.qualifiers:
							if qualifier == 'product' or qualifier == 'gene':
								outputSeqIn.write(str(qualifier) + ' ' + str(gbkFeature.qualifiers[qualifier]) + '\n')
						outputSeqIn.write('\n')
					print '.tbl (Sequin) file created.'

				if ('TRNF' in presentFeatures) or ('tRNA-Phe' in presentFeatures) or ('trnf' in presentFeatures) \
				or ('trnF' in presentFeatures):
					print 'Creating ordered genbank file (with tRNA-Phe at the start)...'
					resultOrderedGbFile = pathOfFinalResults + args.processName + '.gb'

					if 'TRNF' in presentFeatures:
						lookForPhe = 'TRNF'
					elif 'tRNA-Phe' in presentFeatures:
						lookForPhe = 'tRNA-Phe'
					elif 'trnf' in presentFeatures:
						lookForPhe = 'trnf'
					elif 'trnF' in presentFeatures:
						lookForPhe = 'trnF'

					pheAlignment = presentFeatures[lookForPhe][1]
					pheStart = pheAlignment.startBase
					pheEnd = pheAlignment.endBase
					orderedFinalResults = finalResults[pheStart:] + finalResults[0:pheStart]

					with open(resultOrderedGbFile, "w") as outputResult: #create the file!
						count = SeqIO.write(orderedFinalResults, outputResult, "genbank")
						count = SeqIO.write(orderedFinalResults, resultOrderedGbFile.replace('.gb','.fasta'), "fasta")
						createImageOfAnnotation(orderedFinalResults, resultOrderedGbFile.replace('.gb','.png'))
					with open(resultSequinFile.replace('.unordered.tbl','.tbl'),"w") as outputSeqIn:
						outputSeqIn.write('>Features ' + finalResults.name + '\n')
						for gbkFeature in orderedFinalResults.features:
							if gbkFeature.location.strand == 1 or gbkFeature.location.strand == None:
								outputSeqIn.write(str(gbkFeature.location.start + 1) + ' ' + str(gbkFeature.location.end)\
									 	 + ' ' + str(gbkFeature.type) + '\n\t\t')
							else:
								outputSeqIn.write(str(gbkFeature.location.end) + ' ' + str(gbkFeature.location.start)\
										  + ' ' + str(gbkFeature.type) + '\n\t\t')
							for qualifier in gbkFeature.qualifiers:
								if qualifier == 'product' or qualifier == 'gene':
									outputSeqIn.write(str(qualifier) + ' ' + str(gbkFeature.qualifiers[qualifier]) + '\n')
							outputSeqIn.write('\n')
						print 'Ordered .tbl file created.'

					print 'Annotation done. Genbank file created.'
					print ''

				#If circularization couldn't be found 
				'''
				if fourthStep[0] == False and len(presentFeatures) < len(importantFeatures) and args.refSeqFile != None:
					print "Circularization couldn't be found and some features were missing.\nGoing to try and find contigs for the edges..."
					print ''

					#allHitsButBest = SeqIO.parse(open("all_hits_except_best.fasta", "rU"), "fasta", generic_dna)

					#create edge files here
					edgesStart = finalResults[0:args.edgesToLook]
					edgesStart.name = 'edge_start'
					edgesStart.id = 'edge_start'
					edgesEnd = finalResults[(-1) * args.edgesToLook:]
					edgesEnd.name = 'edge_end'
					edgesEnd.id = 'edge_end'
					edgesList = [edgesStart, edgesEnd]
					edgeWrite = SeqIO.write(edgesList, "final_result_edges.fasta", "fasta")

					print "Formatting database for blast to find edges..."
			
					if blastFolder == 'installed':
						command = "formatdb -i " + "final_result_edges.fasta" + " -p F" #need to formatdb edgefile first
					else:
						command = blastFolder + "/bin/makeblastdb -in final_result_edges.fasta -dbtype nucl" #need to formatdb refseq first
					args = shlex.split(command)
					formatDB = Popen(args, stdout=open(os.devnull, 'wb'))
					formatDB.wait()
	
					print "Running blast against finding_edges file to determine if a hit was built..."
					with open("finding_edges.blast.xml",'w') as blastResultFile:
						if blastFolder == 'installed':
							command = "blastall -p blastn -d final_result_edges.fasta -i all_hits_except_best.fasta -e " + str(blasteVal) + " -m 7" #call BLAST with XML output
						else:
							command = blastFolder + "/bin/blastn -task blastn -db final_result_edges.fasta -query all_hits_except_best.fasta -outfmt 5 -evalue " + str(blasteVal) #call BLAST with XML output
						args = shlex.split(command)
						blastAll = Popen(args, stdout=blastResultFile)
						blastAll.wait()
				'''
			else: #if args.ignoreFirstBuildChecks == True
				print ''
				print '--relaxed is turned on, skipping annotation...\n'
				print ''

			if fourthStep[0] == False:
				print " Warning: Circularization wasn't found, so the end of the sequence might have to be extended."

		finalStatsFile.close()
		print '#'*28
		print 'generalMaker is done running.\nMain result files can be found at %s' % pathOfFinalResults
		print 'Check the .stats file in the main results folder for general informations about your assembly.'
		print '#'*28
		print '\nCheck any warnings outputted by generalMaker and then check your final sequence manually for any small changes you might need to make.'
		#cleaning up!
		#try:
		if args.keepTmpFolders == False:
			#remove mira mapping folder
			shutil.rmtree('mira_mapping/',ignore_errors=True)
			#remove mitobim folders
			for x in xrange(args.mitobimIterations):
				shutil.rmtree('iteration' + str(x + 1) + '/',ignore_errors=True)
			#remove kmer folders
			for x in args.kmers.lower().split(','):
				shutil.rmtree('kmer_' + str(x) + '/',ignore_errors=True)
		if args.ignoreFirstBuildChecks == False:
			os.remove('best_query.fasta.nin')
			os.remove('best_query.fasta.nhr')
			os.remove('best_query.fasta.nsq')
			os.remove('important_features.fasta.pin')
			os.remove('important_features.fasta.phr')
			os.remove('important_features.fasta.psq')
		if usingOwnGenBankReference == True:
			tmpRefSeqStuff = args.refSeqFile + '.fasta'
			os.remove(tmpRefSeqStuff)
		else:
			tmpRefSeqStuff = args.refSeqFile
		if args.ignoreFirstBuildChecks == False:
			if args.refSeqFile is not None:
				os.remove(tmpRefSeqStuff + '.nin')
				os.remove(tmpRefSeqStuff + '.nhr')
				os.remove(tmpRefSeqStuff + '.nsq')
				os.remove(args.processName + '.fasta.nin')
				os.remove(args.processName + '.fasta.nhr')
				os.remove(args.processName + '.fasta.nsq')
		os.remove(args.processName + '.fasta')
		if args.skipTrnaScan == False:
			os.remove('tRNAscan.log')
			os.remove('best_query.trnascan')
			shutil.move(args.processName + '.trnascan', pathOfFinalResults + args.processName + '.trnascan')
		if args.skipMitobim == False:
			os.remove('mitobim.log')
		os.remove('circularization_check.blast.xml')
		if args.ignoreFirstBuildChecks == False:
			os.remove('important_features.fasta')
			os.remove('important_features.cds.fasta')
			if args.refSeqFile is not None:
				os.remove('possible_hits.blast.xml')
			if not os.path.exists('mitomaker_tmp'): os.makedirs('mitomaker_tmp')
			shutil.move('best_query.fasta','mitomaker_tmp/best_query.fasta')
			shutil.move('possible_hits.fasta','mitomaker_tmp/possible_hits.fasta')
			shutil.move('important_features.blast.xml','mitomaker_tmp/important_features.blast.xml')
			shutil.move('important_features.cds.blast.xml','mitomaker_tmp/important_features.cds.blast.xml')
		#except:
		#	print 'Could not clean up temporary files.'
		#end clean up.
	else:
		print "After " + str(len(args.kmers.split(','))) + " runs, target DNA wasn't built. Giving up."
		print "Check reads for quality and/or each step output."
