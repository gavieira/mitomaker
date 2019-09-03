#!/usr/bin/python
'''
Script made to test if all required programs for generalMaker are working properly.
'''
from subprocess import Popen
import shlex, os
from Bio import SeqIO, SeqFeature, SeqUtils
from Bio.Alphabet import generic_dna, generic_protein

def testMitoMaker():
	#main mitomaker folder
	module_dir = os.path.dirname(__file__)
	module_dir = os.path.abspath(module_dir)
	testcase_dir = os.path.join(module_dir, 'testcase/')

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
				elif configPart == 'blastfolder':
					blastFolder = line.replace('\n','').replace(' ','').split('=')[-1]
				elif configPart == 'trnascanfolder':
					tRNAscanFolder = line.replace('\n','').replace(' ','').split('=')[-1]
					
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

	if blastFolder.lower() == 'default':
		blastFolder = os.path.join(module_dir, 'blast/')
		print 'WARNING: Blast is still set as default folder. Change it in the config file if you encounter problems running this script.\n'

	if tRNAscanFolder.lower() == 'default':
		tRNAscanFolder = os.path.join(module_dir, 'tRNAscan/')
		print 'WARNING: tRNAscan-SE is still set as default folder. Change it in the config file if you encounter problems running this script.\n'

	print 'Program folders:'
	print 'MITObim = %s' % mitobimFolder
	print 'MITObim1.7 = %s' % newMitobimFolder
	print 'MIRA4 folder = %s' % pathToNewMira
	print 'MIRA3 folder = %s' % pathToOldMira
	print 'SOAPtrans folder = %s' % soapTransFolder
	print 'tRNAscan-SE folder = %s' % tRNAscanFolder
	print ''

	print 'Running mitoMaker.py test case. This might take a while...'
	print 'Creating input file...'
	inputReads = testcase_dir + 'test.fastq'

	with open(testcase_dir + 'test.input', 'w') as inputFile:
		inputFile.write('max_rd_len=100\n\n')
		inputFile.write('[LIB]\n')
		inputFile.write('rank=1\n')
		inputFile.write('asm_flag=3\n')
		inputFile.write('q=' + inputReads + '\n')

	try:
		print 'Checking if blast is working properly..'
		
		if blastFolder.lower() == 'installed':
			command = 'blastall --version'
		else:
			command = blastFolder + 'bin/blastn --help'
		
		args = shlex.split(command)
		testRun = Popen(args)
		testRun.wait()
		print 'Ok'
	except:
		print 'Blast+ does not seem to be properly working.'
		return False

	try:
		print 'Checking if MIRA4 is working properly...'
		if pathToNewMira.lower() == 'installed':
			command = 'mira --version'
		else:
			command = '%sbin/mira --version' % pathToNewMira
		args = shlex.split(command)
		miraRun = Popen(args)
		miraRun.wait()
		print 'Ok'
	except:
		print 'MIRA4 does not seem to be properly working...'
		return False

	try:
		print 'Checking if SOAPdenovo-Trans is working properly...'
		command = '%sSOAPdenovo-Trans-31mer --version' % soapTransFolder
		args = shlex.split(command)
		soapRun = Popen(args)
		soapRun.wait()
		print 'Ok'
	except:
		print 'SOAPdenovo-Trans does not seem to be properly working...'
		return False
	
	try:
		command = "./tRNAscan-SE"
		args = shlex.split(command)
		tRNAscanRun = Popen(args, cwd=tRNAscanFolder)
		tRNAscanRun.wait()
	except:
		print 'tRNAscan-SE failed. Try to recompile it in the tRNAScan-SE folder given by default.'
		return False

	print '\nChecking main script now...\n'
	print 'Running with SOAPdenovo-Trans first...\n'
	
	'''
	Run first with SOAP and checking MITObim as well.
	MITObim shouldn't have problems as it is a simple perl script that is provided in generalMaker's folder
	'''
	try:
		command = 'python -u %s/mitoMaker.py -j test -i test.input -r test_reference.gb -k 23,53' % module_dir
		args = shlex.split(command)
		testRun = Popen(args, cwd=testcase_dir)
		testRun.wait()
	except KeyboardInterrupt:
		raise
		return False
	except:
		print '-'*15
		print 'Failure with SOAP test...\nSkipping to next test...'
		print '-'*15
	else:
		print '\nHuge success! Everything worked fine with SOAP...'

	print 'Running with SOAPdenovo now...\n'
	
	try:
		command = 'python -u %s/mitoMaker.py -j test -i test.input -r test_reference.gb -k 31 --forcedenovo' % module_dir
		args = shlex.split(command)
		testRun = Popen(args, cwd=testcase_dir)
		testRun.wait()
	except KeyboardInterrupt:
		raise
		return False
	except:
		print '-'*15
		print 'Failure with SOAP test...\nSkipping to MIRA4 test...'
		print '-'*15
	else:
		print '\nHuge success! Everything worked fine with SOAPdenovo...'
		
	'''
	Run with SPAdes
	'''
	print 'Running with SPAdes now...\n'
	try:
		command = 'python -u %s/mitoMaker.py -j test -i test.input -r test_reference.gb -k 23 --spades' % module_dir
		args = shlex.split(command)
		testRun = Popen(args, cwd=testcase_dir)
		testRun.wait()
	except KeyboardInterrupt:
		raise
		return False
	except:
		print '-'*15
		print 'Failure with SPAdes test...\nSkipping to next test...'
		print '-'*15
	else:
		print '\nHuge success! Everything worked fine with SOAP...'

	print 'Now running with MIRA, to check if it will work properly and the config file is alright.'

	try:
		command = 'python -u %s/mitoMaker.py -j test -i test.input -r test_reference.gb -o 2 -k default,31 --recursivemira --skipmitobim' % module_dir
		args = shlex.split(command)
		testRun = Popen(args, cwd=testcase_dir)
		testRun.wait()
	except KeyboardInterrupt:
		raise
		return False
	except:
		print 'Failure with MIRA test...'
		return False

	print '\nGreat!\nEverything is working fine.'
	

if __name__ == '__main__':
	testMitoMaker()
