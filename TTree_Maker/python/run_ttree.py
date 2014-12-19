#! /usr/bin/env python

##############################################################################################
##########								   Imports  								##########
##############################################################################################
import os
import sys
from DataFormats.FWLite import Events, Handle
import ROOT
from optparse import OptionParser
from ttree_maker import treemaker

##############################################################################################
##########								Parser Options								##########
##############################################################################################

parser = OptionParser()
#Run options
parser.add_option('-i',  '--input', 	  type = 'string', action='store', default='input',	dest='input',		help='Path to input file holding list of files to run on')
parser.add_option('-og', '--on_grid', 	  type = 'string', action='store', default='yes',	dest='on_grid',		help='Changes everything to relative paths if running on the grid, default is "yes"')
parser.add_option('-me', '--max_events',  type = 'int',    action='store', default=-1,		dest='max_events',  help='Maximum number of events to process (default is -1 for "all")')
parser.add_option('-pe', '--print_every', type = 'int',    action='store', default=1000,	dest='print_every', help='Print progress after how many events?')
parser.add_option('-nj', '--n_jobs', 	  type = 'int',    action='store', default=1,		dest='n_jobs',		help='Number of total grid jobs')
parser.add_option('-ij', '--i_job', 	  type = 'int',    action='store', default=0,		dest='i_job',		help='Which job is this in the sequence?')
#Sample options
parser.add_option('-n',  '--name', 			type = 'string', action='store', 			   dest='name', 		 help='Name of sample or process (used to name output files, etc.)')
parser.add_option('-d',  '--data', 			type = 'string', action='store', default='no', dest='data', 		 help='Set to "yes" if file is a data file (MC file assumed)')
parser.add_option('-cs', '--cross_section', type = 'float',  action='store', default=1.0,  dest='cross_section', help='Cross section of process')
parser.add_option('-ne', '--n_events', 		type = 'float',  action='store', default=1.0,  dest='n_events', 	 help='Number of events generated')
(options, args) = parser.parse_args()

##############################################################################################
##########							Set Up Event Loop								##########
##############################################################################################

print 'Opening files for sample '+options.name+' . . .'  
#Build path to input file
input_files_list = ''
if options.on_grid == 'yes' :
	input_files_list += 'tardir/'
else :
	input_files_list += './'
input_files_list += options.input
if not '.' in options.input :
	input_files_list += '.txt'
print '	Using input file '+input_files_list+''
#open input file in read only mode 
input_files_list = open(input_files_list,'r')
files = []
print '	Getting these files: '
#Read files in line by line
for input_file in input_files_list :
	print '		'+input_file.rstrip()+''
	files.append(input_file.rstrip())
events = Events(files)
ntotalevents = events.size()
#Set filename for analyzer from sample name
filename = options.name
analyzer = treemaker(filename, options.data, options.cross_section/options.n_events)

#Counters
real_count = 0
count = 0

##############################################################################################
##########								Main Event Loop								##########
##############################################################################################

print 'Files opened, starting event loop'
for event in events:
	#increment the "real" counter
	++real_count
	#check the grid split
	if ((realCount-1)-options.i_job) % options.n_jobs != 0 :
		continue
	++count
	#check the max events 
	if count == options.max_events+1 :
		'Processed event number '+str(count-1)+', exiting'
		break
	#print progress
	if count % options.print_every == 0 or count == 1:
			print 'Count at '+str(count)+' out of '+str(ntotalevents/options.n_jobs)+', (%.4f%% complete)'%(float(count) / float(ntotalevents/options.n_jobs) * 100.0)
	#analyze event and add to TTree
	num = analyzer.analyze(event)
	#reset analyzer
	analyzer.reset(num)
#clean up after yourself
del analyzer
