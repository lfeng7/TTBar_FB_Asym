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
ele_pt_bins   = array('d',[40.,45.,50.,55.,60.,70.,80.,100.,150.,200.,300.])
eta_bins 	  = array('d',[-2.5,-1.5,-1.,-0.5,-0.25,0.,0.25,0.5,1.,1.5,2.5])
jet1_pt_bins = array('d',[150.,175.,200.,250.,300.,400.,500.,1000.])
jet2_pt_bins = array('d',[50.,100.,150.,250.,350.,600.])
n_ele_pt_bins = len(ele_pt_bins)-1
n_eta_bins = len(eta_bins)-1
n_jet1_pt_bins = len(jet1_pt_bins)-1
n_jet2_pt_bins = len(jet2_pt_bins)-1
#constants
LUMINOSITY = 19748.

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
#Physics objects
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
mu_trigger 	  = array('I',[2]); 	chain.SetBranchAddress('mu_trigger',mu_trigger)
el_trigger 	  = array('I',[2]); 	chain.SetBranchAddress('el_trigger',el_trigger)
lepb_pt 	  = array('d',[-1.0]);  chain.SetBranchAddress('lepb_pt', lepb_pt)
hadt_pt 	  = array('d',[-1.0]);  chain.SetBranchAddress('hadt_pt', hadt_pt)
hadt_M 		  = array('d',[-1.0]);  chain.SetBranchAddress('hadt_M', hadt_M)
lepW_pt 	  = array('d',[-1.0]);  chain.SetBranchAddress('lepW_pt', lepW_pt)
lept_M 		  = array('d',[-1.0]);  chain.SetBranchAddress('lept_M', lept_M)
#reweighting factors
weight 	  = array('d',[1.0]);  chain.SetBranchAddress('weight', weight)
sf_pileup = array('d',[1.0]);  chain.SetBranchAddress('sf_pileup', sf_pileup)
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
jet1_pt   = array('d',[-1.0]); tree.Branch('jet1_pt',jet1_pt,'jet1_pt/D')
jet2_pt   = array('d',[-1.0]); tree.Branch('jet2_pt',jet2_pt,'jet2_pt/D')
pass_trig = array('I',[2]); tree.Branch('pass_trig',pass_trig,'pass_trig/i')
pass_id   = array('I',[2]); tree.Branch('pass_id',pass_id,'pass_id/i')
#Set up histograms and efficiency graphs
histos_and_graphs = []
ele_pt_all   = TH1D('ele_pt_all','electron p_{T} for all events; p_{T} (GeV)',n_ele_pt_bins,ele_pt_bins); histos_and_graphs.append(ele_pt_all)
ele_eta_all  = TH1D('ele_eta_all','electron #eta for all events; #eta',n_eta_bins,eta_bins); histos_and_graphs.append(ele_eta_all)
jet1_pt_all  = TH1D('jet1_pt_all','jet 1 p_{T} for all events; p_{T} (GeV)',n_jet1_pt_bins,jet1_pt_bins); histos_and_graphs.append(jet1_pt_all)
jet2_pt_all  = TH1D('jet2_pt_all','jet 2 p_{T} for all events; p_{T} (GeV)',n_jet2_pt_bins,jet2_pt_bins); histos_and_graphs.append(jet2_pt_all)
ele_pt_pass  = TH1D('ele_pt_pass','electron p_{T} for passing events; p_{T} (GeV)',n_ele_pt_bins,ele_pt_bins); histos_and_graphs.append(ele_pt_pass)
ele_eta_pass = TH1D('ele_eta_pass','electron #eta for passing events; #eta',n_eta_bins,eta_bins); histos_and_graphs.append(ele_eta_pass)
jet1_pt_pass = TH1D('jet1_pt_pass','jet 1 p_{T} for passing events; p_{T} (GeV)',n_jet1_pt_bins,jet1_pt_bins); histos_and_graphs.append(jet1_pt_pass)
jet2_pt_pass = TH1D('jet2_pt_pass','jet 2 p_{T} for passing events; p_{T} (GeV)',n_jet2_pt_bins,jet2_pt_bins); histos_and_graphs.append(jet2_pt_pass)
ele_pt_x = array('d',n_ele_pt_bins*[0.])
ele_pt_xe = array('d',n_ele_pt_bins*[0.])
ele_pt_y = array('d',n_ele_pt_bins*[0.])
ele_pt_ye = array('d',n_ele_pt_bins*[0.])
ele_eta_x = array('d',n_eta_bins*[0.])
ele_eta_xe = array('d',n_eta_bins*[0.])
ele_eta_y = array('d',n_eta_bins*[0.])
ele_eta_ye = array('d',n_eta_bins*[0.])
jet1_pt_x = array('d',n_jet1_pt_bins*[0.])
jet1_pt_xe = array('d',n_jet1_pt_bins*[0.])
jet1_pt_y = array('d',n_jet1_pt_bins*[0.])
jet1_pt_ye = array('d',n_jet1_pt_bins*[0.])
jet2_pt_x = array('d',n_jet2_pt_bins*[0.])
jet2_pt_xe = array('d',n_jet2_pt_bins*[0.])
jet2_pt_y = array('d',n_jet2_pt_bins*[0.])
jet2_pt_ye = array('d',n_jet2_pt_bins*[0.])
ele_pt_gr    = TGraphErrors(n_ele_pt_bins,ele_pt_x,ele_pt_y,ele_pt_xe,ele_pt_ye); histos_and_graphs.append(ele_pt_gr)
ele_pt_gr.SetName('ele_pt_gr'); ele_pt_gr.SetTitle('Probe efficiency vs. electron p_{T}'); 
ele_pt_gr.GetXaxis().SetName('electron p_{T} (GeV)'); ele_pt_gr.GetYaxis().SetName('Probe efficiency')
ele_eta_gr   = TGraphErrors(n_eta_bins,ele_eta_x,ele_eta_y,ele_eta_xe,ele_eta_ye); histos_and_graphs.append(ele_eta_gr)
ele_eta_gr.SetName('ele_eta_gr'); ele_eta_gr.SetTitle('Probe efficiency vs. electron #eta'); 
ele_eta_gr.GetXaxis().SetName('#eta'); ele_eta_gr.GetYaxis().SetName('Probe efficiency')
jet1_pt_gr   = TGraphErrors(n_jet1_pt_bins,jet1_pt_x,jet1_pt_y,jet1_pt_xe,jet1_pt_ye); histos_and_graphs.append(jet1_pt_gr)
jet1_pt_gr.SetName('jet1_pt_gr'); jet1_pt_gr.SetTitle('Probe efficiency vs.leading jet p_{T}'); 
jet1_pt_gr.GetXaxis().SetName('leading jet p_{T} (GeV)'); jet1_pt_gr.GetYaxis().SetName('Probe efficiency')
jet2_pt_gr   = TGraphErrors(n_jet2_pt_bins,jet2_pt_x,jet2_pt_y,jet2_pt_xe,jet2_pt_ye); histos_and_graphs.append(jet2_pt_gr)
jet2_pt_gr.SetName('jet2_pt_gr'); jet2_pt_gr.SetTitle('Probe efficiency vs. subleading jet p_{T}'); 
jet2_pt_gr.GetXaxis().SetName('subleading jet p_{T} (GeV)'); jet2_pt_gr.GetYaxis().SetName('Probe efficiency')
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
	#muon1 selection and trigger
	cuts.append(mu_trigger[0]==1)
	cuts.append(muon1_pt[0]>ele1_pt[0])
	cuts.append(muon1_pt[0]>40. and abs(muon1_eta[0])<2.4)
	cuts.append(muon1_isLoose[0]==1)
	cuts.append(muon1_relPt[0]>25. or muon1_dR[0]>0.5)
	#muon2 rejection
#	cuts.append(muon2_pt[0]<40. or abs(muon2_eta[0])>2.4 or muon2_isLoose[0]!=1 or (muon2_relPt[0]<25. and muon2_dR[0]<0.5))
	#electron1 selection
	cuts.append(ele1_pt[0]>40. and abs(ele1_eta[0])<2.4)
#	cuts.append(ele1_isLoose[0]==1)
	cuts.append(ele1_relPt[0]>25. or ele1_dR[0]>0.5)
	#electron2 rejection
#	cuts.append(ele2_pt[0]<25. or abs(ele2_eta[0])>2.4 or ele2_isLoose[0]!=1 or (ele2_relPt[0]<25. and ele2_dR[0]<0.5))
	#Require electron and muon to have opposite charges
	cuts.append(ele1_Q[0]!=muon1_Q[0])
	#Other leptonic side cuts
	cuts.append(lepW_pt[0]>50.)
	cuts.append(lept_M[0]>140. and lept_M[0]<250.)
#	#Hadronic side preselection
#	cuts.append(hadt_pt[0]>300. and hadt_M[0]>100.)
	cuts.append(hadt_pt[0]>50. and lepb_pt[0]>50. and max(hadt_pt[0],lepb_pt[0])>150.)
	#check all cuts
	if cuts.count(False) > 0 :
		continue
	#fill the tree
	el_pt[0] = ele1_pt[0]
	el_eta[0] = ele1_eta[0]
	if lepb_pt[0]>hadt_pt[0] :
		jet1_pt[0] = lepb_pt[0]
		jet2_pt[0] = hadt_pt[0]
	else :
		jet2_pt[0] = lepb_pt[0]
		jet1_pt[0] = hadt_pt[0]
	pass_trig[0] = el_trigger[0]
	pass_id[0] = ele1_isLoose[0]
	tree.Fill()
	#fill the histograms
	ele_pt_all.Fill(ele1_pt[0])
	ele_eta_all.Fill(ele1_eta[0])
	jet1_pt_all.Fill(jet1_pt[0])
	jet2_pt_all.Fill(jet2_pt[0])
	if pass_trig[0] == 1 and pass_id[0] == 1 :
#	if pass_trig[0]==1 :
		ele_pt_pass.Fill(ele1_pt[0])
		ele_eta_pass.Fill(ele1_eta[0])
		jet1_pt_pass.Fill(jet1_pt[0])
		jet2_pt_pass.Fill(jet2_pt[0])
#Make the graph y-values and errors
for i in range(n_ele_pt_bins) :
	x_value = (ele_pt_bins[i+1]+ele_pt_bins[i])/2
	x_err   = (ele_pt_bins[i+1]-ele_pt_bins[i])/2
	passing_events = ele_pt_pass.GetBinContent(ele_pt_pass.FindBin(x_value))
	all_events = ele_pt_all.GetBinContent(ele_pt_all.FindBin(x_value))
	y_value = 0.; y_err = 0.
	if all_events > 0. :
		y_value = passing_events/all_events
		y_err = y_value*sqrt(1./passing_events+1./all_events)
	ele_pt_gr.SetPoint(i,x_value,y_value)
	ele_pt_gr.SetPointError(i,x_err,y_err)
for i in range(n_eta_bins) :
	x_value = (eta_bins[i+1]+eta_bins[i])/2
	x_err   = (eta_bins[i+1]-eta_bins[i])/2
	passing_events = ele_eta_pass.GetBinContent(ele_eta_pass.FindBin(x_value))
	all_events = ele_eta_all.GetBinContent(ele_eta_all.FindBin(x_value))
	y_value = 0.; y_err = 0.
	if all_events > 0. :
		y_value = passing_events/all_events
		y_err = y_value*sqrt(1./passing_events+1./all_events)
	ele_eta_gr.SetPoint(i,x_value,y_value)
	ele_eta_gr.SetPointError(i,x_err,y_err)
for i in range(n_jet1_pt_bins) :
	x_value = (jet1_pt_bins[i+1]+jet1_pt_bins[i])/2
	x_err   = (jet1_pt_bins[i+1]-jet1_pt_bins[i])/2
	passing_events = jet1_pt_pass.GetBinContent(jet1_pt_pass.FindBin(x_value))
	all_events = jet1_pt_all.GetBinContent(jet1_pt_all.FindBin(x_value))
	y_value = 0.; y_err = 0.
	if all_events > 0. :
		y_value = passing_events/all_events
		y_err = y_value*sqrt(1./passing_events+1./all_events)
	jet1_pt_gr.SetPoint(i,x_value,y_value)
	jet1_pt_gr.SetPointError(i,x_err,y_err)
for i in range(n_jet2_pt_bins) :
	x_value = (jet2_pt_bins[i+1]+jet2_pt_bins[i])/2
	x_err   = (jet2_pt_bins[i+1]-jet2_pt_bins[i])/2
	passing_events = jet2_pt_pass.GetBinContent(jet2_pt_pass.FindBin(x_value))
	all_events = jet2_pt_all.GetBinContent(jet2_pt_all.FindBin(x_value))
	y_value = 0.; y_err = 0.
	if all_events > 0. :
		y_value = passing_events/all_events
		y_err = y_value*sqrt(1./passing_events+1./all_events)
	jet2_pt_gr.SetPoint(i,x_value,y_value)
	jet2_pt_gr.SetPointError(i,x_err,y_err)

#Write the tree, histograms, and graphs; close the file
outfile.cd()
tree.Write()
for thing in histos_and_graphs :
	thing.Write()
outfile.Close()
	
