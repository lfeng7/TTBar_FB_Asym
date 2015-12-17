#Trigger efficiency calculator originally designed to measure efficiency 
#of the HLT_Ele30_CaloIdVT_TrkIdT_PFNoPUJet100_PFNoPUJet25_v* trigger path
#using SingleMu data events

##########								   Imports  								##########

from ROOT import *
from optparse import OptionParser
from array import array
from math import *
import sys

#Global variables
##Histogram bins
#eta_bins 	  = [-1.*math.pi,math.pi]
#el_pt_bins   = [20.,130.,135.,140.,145.,150.,155.,165.,175.,185.,200.,220.,245.,275.,350.,500.]
#jet1_pt_bins = [100.,700.,1100.]
#jet2_pt_bins = [20.,350.,650.]

##########								Parser Options								##########

parser = OptionParser()
#Run options
parser.add_option('--input', 	  type='string', action='store', default='input', dest='input',	   	  
	help='Path to input file holding list of files to run on')
parser.add_option('--print_every',type='int',    action='store', default=1000,	  dest='print_every', 
	help='Print progress after how many events?')
parser.add_option('--max_events', type='int',    action='store', default=-1,	  dest='max_events',  
	help='Maximum number of events to process (default is -1 for "all")')
parser.add_option('--n_jobs', type='int',    action='store', default=1,	  dest='n_jobs',  
	help='Total number of jobs')
parser.add_option('--i_job', type='int',    action='store', default=0,	  dest='i_job',  
	help='Index of this job')
parser.add_option('--on_grid', type='string',    action='store', default='no',	  dest='on_grid',  
	help='Running on the grid?')
(options, args) = parser.parse_args()


##########									Set up 									##########

print 'Opening files . . . '
#Build path to input file
input_files_list = './'
if options.on_grid == 'yes' :
	input_files_list+='./tardir/'
input_files_list += options.input
if not '.' in options.input :
	input_files_list += '.txt'
print 'Using input file '+input_files_list+''
#open input file in read only mode 
input_files_list = open(input_files_list,'r')
#Set up the chain and add files
chain = TChain('tree')
filenamelist = []
for input_file in input_files_list :
	filenamelist.append(input_file.rstrip())
print 'Chaining files: '
for filename in filenamelist :
	print filename
	chain.Add(filename)
#Get relevant branches
muon1_pt 	  = array('d',[-1.0]);  chain.SetBranchAddress('muon1_pt',muon1_pt)
muon1_eta 	  = array('d',[100.0]); chain.SetBranchAddress('muon1_eta',muon1_eta)
muon1_isLoose = array('I',[2]);  	chain.SetBranchAddress('muon1_isLoose',muon1_isLoose)
muon1_relPt   = array('d',[-1.0]);  chain.SetBranchAddress('muon1_relPt',muon1_relPt)
muon1_dR 	  = array('d',[-1.0]);  chain.SetBranchAddress('muon1_dR',muon1_dR)
muon1_Q 	  = array('i',[0]); 	chain.SetBranchAddress('muon1_Q',muon1_Q)
ele1_pt 	  = array('d',[-1.0]);  chain.SetBranchAddress('ele1_pt',ele1_pt)
ele1_eta 	  = array('d',[100.0]); chain.SetBranchAddress('ele1_eta',ele1_eta)
ele1_isLoose  = array('I',[2]);  	chain.SetBranchAddress('ele1_isLoose',ele1_isLoose)
ele1_relPt 	  = array('d',[-1.0]);  chain.SetBranchAddress('ele1_relPt',ele1_relPt)
ele1_dR 	  = array('d',[-1.0]);  chain.SetBranchAddress('ele1_dR',ele1_dR)
ele1_Q 		  = array('i',[0]); 	chain.SetBranchAddress('ele1_Q',ele1_Q)
muon2_pt 	  = array('d',[-1.0]);  chain.SetBranchAddress('muon2_pt',muon2_pt)
muon2_eta 	  = array('d',[100.0]); chain.SetBranchAddress('muon2_eta',muon2_eta)
muon2_isLoose = array('I',[2]);  	chain.SetBranchAddress('muon2_isLoose',muon2_isLoose)
muon2_relPt   = array('d',[-1.0]);  chain.SetBranchAddress('muon2_relPt',muon2_relPt)
muon2_dR 	  = array('d',[-1.0]);  chain.SetBranchAddress('muon2_dR',muon2_dR)
muon2_Q 	  = array('i',[0]); 	chain.SetBranchAddress('muon2_Q',muon2_Q)
ele2_pt 	  = array('d',[-1.0]);  chain.SetBranchAddress('ele2_pt',ele2_pt)
ele2_eta 	  = array('d',[100.0]); chain.SetBranchAddress('ele2_eta',ele2_eta)
ele2_isLoose  = array('I',[2]);  	chain.SetBranchAddress('ele2_isLoose',ele2_isLoose)
ele2_relPt 	  = array('d',[-1.0]);  chain.SetBranchAddress('ele2_relPt',ele2_relPt)
ele2_dR 	  = array('d',[-1.0]);  chain.SetBranchAddress('ele2_dR',ele2_dR)
ele2_Q 		  = array('i',[0]); 	chain.SetBranchAddress('ele2_Q',ele2_Q)
el_trigger 	  = array('I',[2]); 	chain.SetBranchAddress('el_trigger',el_trigger)
#Set up output file
filename = 'electron_trigger_efficiency'
if options.n_jobs != 1 :
	filename+='_'+str(options.i_job)
filename+='.root'
outfile  = TFile(filename,'recreate')
#Set up output TTree
tree 	  = TTree('tree','tree')
el_pt 	  = array('d',[-1.0]); tree.Branch('el_pt',el_pt,'el_pt/D')
el_eta 	  = array('d',[100.]); tree.Branch('el_eta',el_eta,'el_eta/D')
pass_trig = array('I',[2]); tree.Branch('pass_trig',pass_trig,'pass_trig/i')
#Counters
realcount = 0
count = 0
##########								Main Event Loop								##########
nEntries = chain.GetEntries()
for entry in range(nEntries) :
	#increment the 'real' counter
	realcount+=1
	#check the grid split
	if ((realcount-1)-options.i_job) % options.n_jobs != 0 :
		continue
	#check the max events
	count+=1
	if count == options.max_events+1 :
		print 'Processed event number '+str(count-1)+', exiting'
		break
	#print progress
	if count % options.print_every == 0 or count == 1:
		print 'Count at '+str(count)+' out of '+str(nEntries/options.n_jobs)+', (%.4f%% complete)'%(float(count) / float(nEntries/options.n_jobs) * 100.0)
	chain.GetEntry(entry)
	cuts = []
	#muon1 selection
	cuts.append(muon1_pt[0]>40. and abs(muon1_eta[0])<2.4)
	cuts.append(muon1_isLoose[0]==1)
	cuts.append(muon1_relPt[0]>25. or muon1_dR[0]>0.5)
	#muon2 rejection
	cuts.append(muon2_pt[0]<40. or abs(muon2_eta[0])>2.4 or muon2_isLoose[0]!=1 or (muon2_relPt[0]<25. and muon2_dR[0]<0.5))
	#electron1 selection
	cuts.append(ele1_pt[0]>40. and abs(ele1_eta[0])<2.4)
	cuts.append(ele1_isLoose[0]==1)
	cuts.append(ele1_relPt[0]>25. or ele1_dR[0]>0.5)
	#electron2 rejection
	cuts.append(ele2_pt[0]<25. or abs(ele2_eta[0])>2.4 or ele2_isLoose[0]!=1 or (ele2_relPt[0]<25. and ele2_dR[0]<0.5))
	#Require electron and muon to have opposite charges
	cuts.append(ele1_Q[0]!=muon1_Q[0])
	#check all cuts
	if cuts.count(False) > 0 :
		continue
	el_pt[0] = ele1_pt[0]
	el_eta[0] = ele1_eta[0]
	pass_trig[0] = el_trigger[0]
	tree.Fill()

#Write the tree, close the file
outfile.cd()
tree.Write()
outfile.Close()
	
