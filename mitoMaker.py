#!/usr/bin/python
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
import sys, shlex, os

if __name__ == "__main__":
	module_dir = os.path.dirname(__file__)
	module_dir = os.path.abspath(module_dir)

	flagsToAppend = ''
	if '-o' not in sys.argv:
		flagsToAppend += ' -o 2'
	if '-k' not in sys.argv:
		flagsToAppend += ' -k 23,31,43,53'
	if '--forcedenovo' not in sys.argv:
		flagsToAppend += ' --soaptrans'

	command = 'python -u %s/generalMaker.py %s %s' % (module_dir, ' '.join(sys.argv[1:]), flagsToAppend)
	args = shlex.split(command)
	generalMaker = Popen(args)
	generalMaker.wait()
