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
import shlex, os, shutil

def miraMapping(processName, processorsToUse, inputFile, miraTechnology, useNewMira, pathToNewMira, pathToOldMira, pairedEnd, copyKmers, lastKmer):
	#let's create a folder for mira!
	print 'Creating folder for MIRA Mapping on phase 1 best sequence...'
	pathToWork = 'mira_mapping/'
	shutil.rmtree('mira_mapping/',ignore_errors=True)
	#if not os.path.exists(pathToWork): os.makedirs(pathToWork)
	os.makedirs(pathToWork)
	
	#now, let's read the input file and get the reads from there...
	stringOfReadsFiles = ''
	with open(inputFile, 'r') as soapInputFile:
		for line in soapInputFile:
			if '.fastq' in line:
				#listOfReadsFiles.append(line.replace('\n',''))
				stringOfReadsFiles += ' ' + line.replace('\n','').split('=')[-1]

	if useNewMira == False: #use mira 3.4, default behaviour
		#grabbed all input files, let's cat them into one since it's mapping...
		with open(pathToWork + processName + '_in.' + miraTechnology + '.fastq', 'w') as miraInput:
			command = "cat" + stringOfReadsFiles
			args = shlex.split(command)
			createUniqueReadsFile = Popen(args, stdout=miraInput)
			createUniqueReadsFile.wait()
	
		#let's create the backbone file
		destFile = pathToWork + processName + '_backbone_in.fasta'
		shutil.copyfile('best_query.fasta', destFile)
	
		#time to run mira mapping
		print ''
		print 'Now running MIRA Mapping. Logs will be saved in mira_mapping/'
		print ''
		try:
			#no need to create manifest file
			with open(pathToWork + 'mira.log','w') as miraLogFile:
				if lastKmer != 'default' and copyKmers == True:
					command = pathToOldMira + '/bin/mira --project=' + processName + ' --job=mapping,genome,accurate,' + miraTechnology + ' -SK:not=' + str(processorsToUse) + ':bph=' + lastKmer + ' -MI:somrnl=0:sonfs=no -SB:bft=fasta:bbq=22:bsn=' + processName
				else:
					command = pathToOldMira + '/bin/mira --project=' + processName + ' --job=mapping,genome,accurate,' + miraTechnology + ' -SK:not=' + str(processorsToUse) + ' -MI:somrnl=0:sonfs=no -SB:bft=fasta:bbq=22:bsn=' + processName
				args = shlex.split(command)
				miraRun = Popen(args, cwd=pathToWork, stdout=miraLogFile)
				miraRun.wait()
				return True
		except:
			print "Mira mapping with MIRA3.4 failed. Check it's logs for more information."
			return False
	else: #use mira4

		#cat all reads to be used by mitobim, but for the mapping assembly, different readgroups will be created
		with open(pathToWork + processName + '_in.' + miraTechnology + '.fastq', 'w') as miraInput:
			command = "cat" + stringOfReadsFiles
			args = shlex.split(command)
			createUniqueReadsFile = Popen(args, stdout=miraInput)
			createUniqueReadsFile.wait()

		with open(pathToWork + 'mapping.manifest', 'w') as manifestFile:
			manifestFile.write('project = ' + processName + '\n')
			manifestFile.write('job = genome, mapping, accurate\n')
			if lastKmer == 'default' or copyKmers == False:
				manifestFile.write('parameters = -NW:cmrnl=warn:cac=warn:cnfs=warn -GE:not=' + str(processorsToUse) + ' -SB:tor=no\n\n')
			elif copyKmers == True:
				manifestFile.write('parameters = -NW:cmrnl=warn:cac=warn:cnfs=warn -GE:not=' + str(processorsToUse) + ' -SB:tor=no -SK:bph=' + lastKmer + '\n\n')
			manifestFile.write('readgroup\nis_reference\n')

			#let's create the backbone file
			destFile = pathToWork + processName + '_backbone_in.fna'
			shutil.copyfile('best_query.fasta', destFile)

			#write the backbone part of the manifest
			manifestFile.write('data = ' + processName + '_backbone_in.fna\n')
			manifestFile.write('default_qual = 20\n')
			manifestFile.write('technology = text\nstrain = backbone\n\n')

			with open(inputFile, 'r') as soapInputFile:
				listOfInputs = {}
				listOfTechnologies = {}
				listOfQualities = {}
				n = 0
				for soapLine in soapInputFile:
					if '[LIB]' in soapLine: #start of a new library
						n += 1
						if n > 7:
							print "WARNING: MIRA won't work with more than 7 libraries, going to procceed \nwith the first 7 libraries only."
							break
						listOfInputs[n] = []
						listOfTechnologies[n] = miraTechnology
						listOfQualities[n] = None
					if 'q' == soapLine[0] or 'f' == soapLine[0]: #append the data files
						listOfInputs[n].append(soapLine.replace('\n','').split('=')[-1]) #grab only the file part of this line
					elif 'technology' in soapLine:
						listOfTechnologies[n] = soapLine.replace('\n','').split('=')[-1]
					elif 'default_qual' in soapLine:
						listOfQualities[n] = soapLine.replace('\n','').split('=')[-1]
						
				#time to grab the readgroups from the dict
				for readGroup in listOfInputs:
					manifestFile.write('readgroup\n')
					manifestFile.write('data =')
					for readFile in listOfInputs[readGroup]:
						manifestFile.write(' ' + readFile)
					manifestFile.write('\n')
					manifestFile.write('technology = ' + listOfTechnologies[readGroup] + '\n')
					if listOfQualities[readGroup] != None:
						manifestFile.write('default_qual = ' + listOfQualities[readGroup] + '\n')
					manifestFile.write('strain = ' + processName + '_' + str(readGroup) + '\n')
					if len(listOfInputs[readGroup]) > 1:
						if pairedEnd == False: #user did not ask for paired info to be automatically discovered
							manifestFile.write('template_size = unknown infoonly\nsegment_placement = unknown infoonly\n\n')
						else:
							manifestFile.write('autopairing\n\n')
					else:
						manifestFile.write('\n')

		print ''
		print 'Now running MIRA Mapping with MIRA4. Logs will be saved in mira_mapping/'
		print ''

		with open(pathToWork + 'mira.log','w') as miraLogFile:
			try:
				if pathToNewMira.lower() != 'installed':
					command = pathToNewMira + '/bin/mira mapping.manifest'
				else:
					command = 'mira mapping.manifest'
				args = shlex.split(command)
				miraRun = Popen(args, cwd=pathToWork, stdout=miraLogFile)
				miraRun.wait()
				return True
			except:
				print "An error occured while running MIRA4 for mapping assembly. Check it's logs for more information."
				print 'Printing last 10 lines of log file...\n'
				with open(pathToWork + 'mira.log','r') as logFile:
					content = logFile.readlines()
					for n in xrange(-1,-11,-1):
						print content[n]
				print ''
				return False
