#run_ttree.py: helper code to run the tree_maker workhorse code that produces ttrees from nTuples
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

import os
import sys
from DataFormats.FWLite import Events, Handle
import ROOT
from optparse import OptionParser
from ttree_maker import treemaker

##########								Parser Options								##########

parser = OptionParser()
#Run options
parser.add_option('--input', 	  type='string', action='store', default='input', dest='input',	   	  
	help='Path to input file holding list of files to run on')
parser.add_option('--on_grid', 	  type='string', action='store', default='no',	  dest='on_grid',	  
	help='Changes everything to relative paths if running on the grid, default is "no"')
parser.add_option('--event_type', type='string', action='store', default='none',  dest='event_type',  
	help='ttbar event type: "qq_semilep","gg_semilep","dilep", or "had", or "none" for background (default "none")')
parser.add_option('--sideband',	  type='string', action='store', default='no',    dest='sideband',    
	help='Perform analysis in a sideband? "yes" or "no" (default: "no")')
parser.add_option('--leptons', 	  type='string', action='store', default='muons', dest='leptons',	  
	help='Specify "muons" or "electrons" lepton type, default is "muons"')
parser.add_option('--top_type',   type='int', 	 action='store', default=1, 	  dest='top_type',	  
	help='Specify top type "1" (fully merged) or "2" (partially merged), default is "1"')
parser.add_option('--max_events', type='int',    action='store', default=-1,	  dest='max_events',  
	help='Maximum number of events to process (default is -1 for "all")')
parser.add_option('--print_every',type='int',    action='store', default=1000,	  dest='print_every', 
	help='Print progress after how many events?')
parser.add_option('--n_jobs', 	  type='int',    action='store', default=1,		  dest='n_jobs',	  
	help='Number of total grid jobs')
parser.add_option('--i_job', 	  type='int',    action='store', default=0,		  dest='i_job',	   	  
	help='Which job is this in the sequence?')
#Sample options
parser.add_option('--name', 		 type='string', action='store', 			  	dest='name', 		    
	help='Name of sample or process (used to name output files, etc.)')
parser.add_option('--data', 		 type='string', action='store', default='no', 	dest='data', 		    
	help='Set to "yes" if file is a data file (MC file assumed)')
parser.add_option('--generator', 	 type='string', action='store', default='none', dest='generator', 		
	help='Monte Carlo generator for this file (powheg, madgraph, pythia8); default is "none"')
parser.add_option('--cross_section', type='float',  action='store', default=1.0,  	dest='cross_section', 	
	help='Cross section of process')
parser.add_option('--n_events', 	 type='float',  action='store', default=1.0,  	dest='n_events', 	    
	help='Number of events generated')
(options, args) = parser.parse_args()

##########							Set Up Event Loop								##########

#append 'sb' to the name for a sideband file
if options.sideband.lower() == 'yes' :
	options.name = options.name+'_sb'

#append the top type to the name
if options.top_type == 1 :
	options.name = options.name+'_type1'
elif options.top_type == 2 :
	options.name = options.name+'_type2'

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
print 'Using input file '+input_files_list+''
#open input file in read only mode 
input_files_list = open(input_files_list,'r')
files = []
print 'Getting these files: '
#Read files in line by line
for input_file in input_files_list :
	print '	'+input_file.rstrip()+''
	files.append(input_file.rstrip())
events = Events(files)
ntotalevents = events.size()
#Set filename for analyzer from sample name
filename = options.name
if options.n_jobs>1 :
	filename+='_'+str(options.i_job)
filename+='_tree.root'
#Initialize analyzer
analyzer = treemaker(filename, options.data, options.generator, options.event_type, options.sideband, 
	options.leptons, options.top_type, options.cross_section/options.n_events,options.on_grid)

#Counters
real_count = 0
count = 0

##########								Main Event Loop								##########

print 'Files opened, starting event loop'
for event in events:
	#increment the "real" counter
	real_count+=1
	#check the grid split
	if ((real_count-1)-options.i_job) % options.n_jobs != 0 :
		continue
	count+=1
	#check the max events 
	if count == options.max_events+1 :
		print 'Processed event number '+str(count-1)+', exiting'
		break
	#print progress
	if count % options.print_every == 0 or count == 1:
			print ( 'Count at '+str(count)+' out of '+str(ntotalevents/options.n_jobs)+', (%.4f%% complete)'
				%(float(count) / float(ntotalevents/options.n_jobs) * 100.0) )
	#analyze event and add to TTree
	err = analyzer.analyze(event)
	#reset analyzer
	analyzer.reset(err)
#clean up after yourself
del analyzer