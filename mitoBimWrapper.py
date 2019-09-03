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

def mitoBimWrapper(mitobimIterations, processName, miraTechnology, mitobimFolder, readLen, newMira, newMitobimFolder, pathToNewMira, pathToOldMira,
                  kmerUsed):
	'''
	Wrapper for MITObim. Just runs the program.
	'''
	#let's clean up, just in case
	print ''
	print 'Deleting any iterationX/ folders from previous MITObim runs that will be used by this one...'
	print ''
	
	for x in xrange(mitobimIterations):
		shutil.rmtree('iteration' + str(x + 1) + '/',ignore_errors=True)
	
	#now, time to create the mitobim call line...
	#grab the absolute path of both the backbone and the readpool, as mitobim doesn't work well without absolute paths
	print "Checking readpool and maf files for MITObim..."
	readPoolPath = os.path.abspath('mira_mapping/' + processName + '_in.' + miraTechnology + '.fastq')
	mafPath = os.path.abspath('mira_mapping/' + processName + '_assembly/' + processName + '_d_results/'
							+ processName + '_out.maf')
	
	#debugging messages
	#print 'Readpool: ' + readPoolPath
	#print 'Maf: ' + mafPath
	print 'Running MITObim.pl...'

	try:
		#time to run mitobim
		if newMira == True:
			#run with mitoBIM 1.7, if mapping was made with mira4.0
			with open('mitobim.log','w') as mitobimLogFile:
				command = 'perl ' + newMitobimFolder + 'MITObim.pl -start 1 -end ' + str(mitobimIterations) + ' -sample ' + processName + '_1' ' -ref backbone -readpool ' + readPoolPath + ' -kmer ' + str(kmerUsed) + ' -maf ' + mafPath + ' --clean --readlength ' + str(readLen) + ' --mirapath ' + pathToNewMira + '/bin/'
				args = shlex.split(command)
				mitobimCaller = Popen(args, stdout=mitobimLogFile)
				mitobimCaller.wait()
				return True
		else:
			#run with mitoBIM 1.6
			with open('mitobim.log','w') as mitobimLogFile:
				command = 'perl ' + mitobimFolder + 'MITObim.pl -start 1 -end ' + str(mitobimIterations) + ' -strain ' + processName + ' -ref ReferenceStrain -readpool ' + readPoolPath + ' -maf ' + mafPath + ' --clean --readlength ' + str(readLen) + ' --mirapath ' + pathToOldMira + '/bin/'
				args = shlex.split(command)
				mitobimCaller = Popen(args, stdout=mitobimLogFile)
				mitobimCaller.wait()
				return True
	except:
		print ''
		print "MITObim failed. Check it's logs for more information..."
		print 'Printing last 10 lines of log file and aborting.'
		with open('mitobim.log','r') as logFile:
			content = logFile.readlines()
			for n in xrange(-1,-11,-1):
				print content[n]
			print ''
		return False
