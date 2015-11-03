#Trigger efficiency calculator originally designed to measure efficiency 
#of the HLT_Ele30_CaloIdVT_TrkIdT_PFNoPUJet100_PFNoPUJet25_v* trigger path

##########								   Imports  								##########

from ROOT import *
from DataFormats.FWLite import Events, Handle
from optparse import OptionParser
from array import array
from math import *
import sys

#Global variables
##Histogram bins
#eta_bins 	 = [-1.0*math.pi,-2.1,-1.6,-1.2,-0.9,-0.6,-0.3,-0.2,0.2,0.3,0.6,0.9,1.2,1.6,2.1,math.pi]
#el_pt_bins 	 = [0.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,90.,100.,125.,150.,175.,200.,250.,300.,400.,550.]
#jet1_pt_bins = [50.,100.,125.,150.,175.,200.,250.,300.,350.,400.,450.,500.,550.,600.,700.,800.,1200.]
#jet2_pt_bins = [0.,25.,50.,75.,100.,125.,150.,175.,200.,250.,275.,300.,350.,400.,500.,700.]
#Trigger paths
REF_TRIG_PATH = 'HLT_Ele27_WP80_v'
#REF_TRIG_PATH = 'HLT_Ele8_CaloIdT_TrkIdVL_v'
#REF_TRIG_PATH = 'HLT_Ele8_CaloIdT_TrkIdVL_Jet30_v3'
TRIG_PATH = 'HLT_Ele30_CaloIdVT_TrkIdT_PFNoPUJet100_PFNoPUJet25_v'

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

##################################  Handles and Labels  ##################################

vector_of_4vecs = 'vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >'
#Trigger Bit Information
trigHandle = Handle('edm::TriggerResults'); 	trigLabel = ('TriggerResults','','HLT')
#muons
muLabel = ('jhuMuonPFlowLoose','muonLoose'); 	muHandle = Handle(vector_of_4vecs)
#electrons
elLabel = ('jhuElePFlowLoose','electronLoose'); elHandle = Handle(vector_of_4vecs)
#Jets
jetLabel = ('jhuAk5','AK5'); 				jetHandle = Handle(vector_of_4vecs)

##########							Set Up Event Loop								##########

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
files = []
print 'Getting these files: '
#Read files in line by line
for input_file in input_files_list :
	print '	'+input_file.rstrip()+''
	files.append(input_file.rstrip())
events = Events(files)
ntotalevents = events.size()
#Set up output file
filename = TRIG_PATH+'_efficiency'
if options.n_jobs != 1 :
	filename+='_'+str(options.i_job)
filename+='.root'
outfile  = TFile(filename,'recreate')
#Set up output TTree
tree 	  = TTree('tree','tree')
el_pt 	  = array('d',[-1.0]); tree.Branch('el_pt',el_pt,'el_pt/D')
el_eta 	  = array('d',[100.]); tree.Branch('el_eta',el_eta,'el_eta/D')
jet1_pt   = array('d',[-1.0]); tree.Branch('jet1_pt',jet1_pt,'jet1_pt/D')
jet2_pt   = array('d',[-1.0]); tree.Branch('jet2_pt',jet2_pt,'jet2_pt/D')
pass_trig = array('I',[2]); tree.Branch('pass_trig',pass_trig,'pass_trig/i')
pass_ref  = array('I',[2]); tree.Branch('pass_ref',pass_ref,'pass_ref/i')
#Counters
realcount = 0
count = 0

##########								Main Event Loop								##########

print 'Files opened, starting event loop'
for event in events:
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
		print 'Count at '+str(count)+' out of '+str(ntotalevents/options.n_jobs)+', (%.4f%% complete)'%(float(count) / float(ntotalevents/options.n_jobs) * 100.0)

	#Read Trigger information
	event.getByLabel(trigLabel,trigHandle)
	if not trigHandle.isValid() :
		print 'trigger handle is invalid, fool, crashing now, sorry : ('
		break
	trigResults = trigHandle.product()
	trigNames = event.object().triggerNames(trigResults)
	pass_trig[0] = 2
	pass_ref[0] = 2
	for i in range(trigResults.size()) :
		s = str(trigNames.triggerName(i))
		if s.startswith(TRIG_PATH) :
			if trigResults.accept(i) :
				pass_trig[0] = 1
			else :
				pass_trig[0] = 0
			if pass_ref[0] != 2 :
				break
		elif s.startswith(REF_TRIG_PATH) :
			if trigResults.accept(i) :
				pass_ref[0] = 1
			else :
				pass_ref[0] = 0
			if pass_trig[0] != 2 :
				break
	#make sure it passed at least one trigger
	if pass_trig[0] != 1 and pass_ref[0] != 1 :
		continue

	#get all the info from the event
	#jets
	event.getByLabel(jetLabel,jetHandle)
	if not jetHandle.isValid() :
		print 'jet handle is invalid, crashing : ('
		break
	jetVars = jetHandle.product()
	jets = []
	for jetVar in jetVars :
		newJet = TLorentzVector(); newJet.SetPtEtaPhiM(jetVar.Pt(),jetVar.Eta(),jetVar.Phi(),jetVar.M())
		jets.append(newJet)
	jets.sort(key = lambda x: x.Pt(),reverse=True)
	#muons
	event.getByLabel(muLabel,muHandle)
	if not muHandle.isValid() :
		print 'muon handle is invalid, crashing : ('
		break
	muonVars = muHandle.product()
	muons = []
	for muonVar in muonVars :
		newMuon = TLorentzVector(); newMuon.SetPtEtaPhiM(muonVar.Pt(),muonVar.Eta(),muonVar.Phi(),muonVar.M())
		muons.append(newMuon)
	muons.sort(key = lambda x: x.Pt(),reverse=True)
	#electrons
	event.getByLabel(elLabel,elHandle)
	if not elHandle.isValid() :
		print 'electron handle is invalid, crashing : ('
		break
	elVars = elHandle.product()
	els = []
	for elVar in elVars :
		newEl = TLorentzVector(); newEl.SetPtEtaPhiM(elVar.Pt(),elVar.Eta(),elVar.Phi(),elVar.M())
		els.append(newEl)
	els.sort(key = lambda x: x.Pt(),reverse=True)

	#write into the TTree
	if len(els) > 0 :
		el_pt[0]   = els[0].Pt()
		el_eta[0]  = els[0].Eta()
	if len(jets) > 0 :
		jet1_pt[0] = jets[0].Pt()
	if len(jets) > 1 :
		jet2_pt[0] = jets[1].Pt()
	tree.Fill()

#Write the tree, close the file
outfile.cd()
tree.Write()
outfile.Close()
	
