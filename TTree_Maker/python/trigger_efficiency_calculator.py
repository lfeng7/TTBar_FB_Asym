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
ele_pt_bins   = array('d',[40.,50.,60.,80.,100.,150.,200.,300.])
eta_bins 	  = array('d',[-2.5,-1.,-0.5,0.,0.5,1.,2.5])
ele_pt_bins_2D   = array('d',[40.,50.,60.,80.,100.,150.,300.])
eta_bins_2D 	 = array('d',[0.,0.75,1.4,1.6,2.5])
jet1_pt_bins = array('d',[150.,200.,300.,500.,1000.])
jet2_pt_bins = array('d',[50.,100.,150.,200.,500.])
n_ele_pt_bins = len(ele_pt_bins)-1
n_eta_bins = len(eta_bins)-1
n_ele_pt_bins_2D = len(ele_pt_bins_2D)-1
n_eta_bins_2D = len(eta_bins_2D)-1
n_jet1_pt_bins = len(jet1_pt_bins)-1
n_jet2_pt_bins = len(jet2_pt_bins)-1
#constants
LUMINOSITY = 19748.
CSVL = 0.244
CSVM = 0.679
CSVT = 0.898

##########								Parser Options								##########

parser = OptionParser()
#Run options
parser.add_option('--input', 	  type='string', action='store', default='input', dest='input',	   	  
	help='Path to input file holding list of files to run on')
parser.add_option('--print_every',type='int',    action='store', default=1000,	  dest='print_every', 
	help='Print progress after how many events?')
parser.add_option('--max_events', type='int',    action='store', default=-1,	  dest='max_events',  
	help='Maximum number of events to process (default is -1 for "all")')
(options, args) = parser.parse_args()

##########							Filegroup Class 							##########

class filegroup : 
	#docstring
	"""filegroup class"""
	
	#__init__function
	def __init__(self,name,chain,outputfile,mode) :
		self.name = name
		self.chain = chain
		#Get relevant branches
		#Physics objects
		muon1_pt 	  = array('d',[-1.0]);  self.chain.SetBranchAddress('muon1_pt',muon1_pt)
		muon1_eta 	  = array('d',[100.0]); self.chain.SetBranchAddress('muon1_eta',muon1_eta)
		muon1_isLoose = array('I',[2]);  	self.chain.SetBranchAddress('muon1_isLoose',muon1_isLoose)
		muon1_relPt   = array('d',[-1.0]);  self.chain.SetBranchAddress('muon1_relPt',muon1_relPt)
		muon1_dR 	  = array('d',[-1.0]);  self.chain.SetBranchAddress('muon1_dR',muon1_dR)
		muon1_Q 	  = array('i',[0]); 	self.chain.SetBranchAddress('muon1_Q',muon1_Q)
		ele1_pt 	  = array('d',[-1.0]);  self.chain.SetBranchAddress('ele1_pt',ele1_pt)
		ele1_eta 	  = array('d',[100.0]); self.chain.SetBranchAddress('ele1_eta',ele1_eta)
		ele1_phi 	  = array('d',[-1.0]);  self.chain.SetBranchAddress('ele1_phi',ele1_phi)
		ele1_M 		  = array('d',[-1.0]);  self.chain.SetBranchAddress('ele1_M',ele1_M)
		ele1_isLoose  = array('I',[2]);  	self.chain.SetBranchAddress('ele1_isLoose',ele1_isLoose)
		ele1_relPt 	  = array('d',[-1.0]);  self.chain.SetBranchAddress('ele1_relPt',ele1_relPt)
		ele1_dR 	  = array('d',[-1.0]);  self.chain.SetBranchAddress('ele1_dR',ele1_dR)
		ele1_Q 		  = array('i',[0]); 	self.chain.SetBranchAddress('ele1_Q',ele1_Q)
		muon2_pt 	  = array('d',[-1.0]);  self.chain.SetBranchAddress('muon2_pt',muon2_pt)
		muon2_eta 	  = array('d',[100.0]); self.chain.SetBranchAddress('muon2_eta',muon2_eta)
		muon2_isLoose = array('I',[2]);  	self.chain.SetBranchAddress('muon2_isLoose',muon2_isLoose)
		muon2_relPt   = array('d',[-1.0]);  self.chain.SetBranchAddress('muon2_relPt',muon2_relPt)
		muon2_dR 	  = array('d',[-1.0]);  self.chain.SetBranchAddress('muon2_dR',muon2_dR)
		muon2_Q 	  = array('i',[0]); 	self.chain.SetBranchAddress('muon2_Q',muon2_Q)
		ele2_pt 	  = array('d',[-1.0]);  self.chain.SetBranchAddress('ele2_pt',ele2_pt)
		ele2_eta 	  = array('d',[100.0]); self.chain.SetBranchAddress('ele2_eta',ele2_eta)
		ele2_phi 	  = array('d',[-1.0]);  self.chain.SetBranchAddress('ele2_phi',ele2_phi)
		ele2_M 		  = array('d',[-1.0]);  self.chain.SetBranchAddress('ele2_M',ele2_M)
		ele2_isLoose  = array('I',[2]);  	self.chain.SetBranchAddress('ele2_isLoose',ele2_isLoose)
		ele2_relPt 	  = array('d',[-1.0]);  self.chain.SetBranchAddress('ele2_relPt',ele2_relPt)
		ele2_dR 	  = array('d',[-1.0]);  self.chain.SetBranchAddress('ele2_dR',ele2_dR)
		ele2_Q 		  = array('i',[0]); 	self.chain.SetBranchAddress('ele2_Q',ele2_Q)
		mu_trigger 	  = array('I',[2]); 	self.chain.SetBranchAddress('mu_trigger',mu_trigger)
		el_trigger 	  = array('I',[2]); 	self.chain.SetBranchAddress('el_trigger',el_trigger)
		lepb_pt 	  = array('d',[-1.0]);  self.chain.SetBranchAddress('lepb_pt', lepb_pt)
		lepb_csv 	  = array('d',[-1.0]);  self.chain.SetBranchAddress('lepb_csv', lepb_csv)
		hadt_pt 	  = array('d',[-1.0]);  self.chain.SetBranchAddress('hadt_pt', hadt_pt)
		hadt_csv 	  = array('d',[-1.0]);  self.chain.SetBranchAddress('hadt_csv', hadt_csv)
		hadt_M 		  = array('d',[-1.0]);  self.chain.SetBranchAddress('hadt_M', hadt_M)
		lepW_pt 	  = array('d',[-1.0]);  self.chain.SetBranchAddress('lepW_pt', lepW_pt)
		lept_M 		  = array('d',[-1.0]);  self.chain.SetBranchAddress('lept_M', lept_M)
		#reweighting factors
		weight 	  = array('d',[1.0]);  chain.SetBranchAddress('weight', weight)
		sf_pileup = array('d',[1.0]);  chain.SetBranchAddress('sf_pileup', sf_pileup)
		sf_top_pT = array('d',[1.0]);  chain.SetBranchAddress('sf_top_pT', sf_top_pT)
		#Set up output TTree
		self.tree = TTree('tree','tree')
		el_pt 	  = array('d',[-1.0]); self.tree.Branch('el_pt',el_pt,'el_pt/D')
		el_eta 	  = array('d',[100.]); self.tree.Branch('el_eta',el_eta,'el_eta/D')
		jet1_pt   = array('d',[-1.0]); self.tree.Branch('jet1_pt',jet1_pt,'jet1_pt/D')
		jet2_pt   = array('d',[-1.0]); self.tree.Branch('jet2_pt',jet2_pt,'jet2_pt/D')
		pass_trig = array('I',[2]); self.tree.Branch('pass_trig',pass_trig,'pass_trig/i')
		pass_id   = array('I',[2]); self.tree.Branch('pass_id',pass_id,'pass_id/i')
		#Set up histograms and efficiency graphs
		self.histos_and_graphs = []
		self.ele_pt_all   = TH1D(self.name+'_ele_pt_all',self.name+'electron p_{T} for all events; p_{T} (GeV)',n_ele_pt_bins,ele_pt_bins) 
		self.histos_and_graphs.append(self.ele_pt_all)
		self.ele_eta_all  = TH1D(self.name+'_ele_eta_all',self.name+'electron #eta for all events; #eta',n_eta_bins,eta_bins) 
		self.histos_and_graphs.append(self.ele_eta_all)
		self.jet1_pt_all  = TH1D(self.name+'_jet1_pt_all',self.name+'jet 1 p_{T} for all events; p_{T} (GeV)',n_jet1_pt_bins,jet1_pt_bins) 
		self.histos_and_graphs.append(self.jet1_pt_all)
		self.jet2_pt_all  = TH1D(self.name+'_jet2_pt_all',self.name+'jet 2 p_{T} for all events; p_{T} (GeV)',n_jet2_pt_bins,jet2_pt_bins) 
		self.histos_and_graphs.append(self.jet2_pt_all)
		self.ele_pt_pass  = TH1D(self.name+'_ele_pt_pass',self.name+'electron p_{T} for passing events; p_{T} (GeV)',n_ele_pt_bins,ele_pt_bins) 
		self.histos_and_graphs.append(self.ele_pt_pass)
		self.ele_eta_pass = TH1D(self.name+'_ele_eta_pass',self.name+'electron #eta for passing events; #eta',n_eta_bins,eta_bins) 
		self.histos_and_graphs.append(self.ele_eta_pass)
		self.jet1_pt_pass = TH1D(self.name+'_jet1_pt_pass',self.name+'jet 1 p_{T} for passing events; p_{T} (GeV)',n_jet1_pt_bins,jet1_pt_bins) 
		self.histos_and_graphs.append(self.jet1_pt_pass)
		self.jet2_pt_pass = TH1D(self.name+'_jet2_pt_pass',self.name+'jet 2 p_{T} for passing events; p_{T} (GeV)',n_jet2_pt_bins,jet2_pt_bins) 
		self.histos_and_graphs.append(self.jet2_pt_pass)
		self.histo_2d_all  = TH2D(self.name+'_histo_2D_all','',n_ele_pt_bins_2D,ele_pt_bins_2D,n_eta_bins_2D,eta_bins_2D); self.histos_and_graphs.append(self.histo_2d_all)
		self.histo_2d_pass = TH2D(self.name+'_histo_2D_pass','',n_ele_pt_bins_2D,ele_pt_bins_2D,n_eta_bins_2D,eta_bins_2D); self.histos_and_graphs.append(self.histo_2d_pass)
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
		self.ele_pt_gr    = TGraphErrors(n_ele_pt_bins,ele_pt_x,ele_pt_y,ele_pt_xe,ele_pt_ye); self.histos_and_graphs.append(self.ele_pt_gr)
		self.ele_pt_gr.SetName(self.name+'ele_pt_gr'); self.ele_pt_gr.SetTitle(self.name+' probe efficiency vs. electron p_{T}'); 
		self.ele_pt_gr.GetXaxis().SetName('electron p_{T} (GeV)'); self.ele_pt_gr.GetYaxis().SetName('Probe efficiency')
		self.ele_eta_gr   = TGraphErrors(n_eta_bins,ele_eta_x,ele_eta_y,ele_eta_xe,ele_eta_ye); self.histos_and_graphs.append(self.ele_eta_gr)
		self.ele_eta_gr.SetName(self.name+'ele_eta_gr'); self.ele_eta_gr.SetTitle(self.name+' probe efficiency vs. electron #eta'); 
		self.ele_eta_gr.GetXaxis().SetName('#eta'); self.ele_eta_gr.GetYaxis().SetName('Probe efficiency')
		self.jet1_pt_gr   = TGraphErrors(n_jet1_pt_bins,jet1_pt_x,jet1_pt_y,jet1_pt_xe,jet1_pt_ye); self.histos_and_graphs.append(self.jet1_pt_gr)
		self.jet1_pt_gr.SetName(self.name+'jet1_pt_gr'); self.jet1_pt_gr.SetTitle(self.name+' probe efficiency vs.leading jet p_{T}'); 
		self.jet1_pt_gr.GetXaxis().SetName('leading jet p_{T} (GeV)'); self.jet1_pt_gr.GetYaxis().SetName('Probe efficiency')
		self.jet2_pt_gr   = TGraphErrors(n_jet2_pt_bins,jet2_pt_x,jet2_pt_y,jet2_pt_xe,jet2_pt_ye); self.histos_and_graphs.append(self.jet2_pt_gr)
		self.jet2_pt_gr.SetName(self.name+'jet2_pt_gr'); self.jet2_pt_gr.SetTitle(self.name+' probe efficiency vs. subleading jet p_{T}'); 
		self.jet2_pt_gr.GetXaxis().SetName('subleading jet p_{T} (GeV)'); self.jet2_pt_gr.GetYaxis().SetName('Probe efficiency')
		#Counter
		count = 0
		##########								Main Event Loop								##########
		print 'Filling trees for '+self.name+'. . .'
		nEntries = chain.GetEntries()
		for entry in range(nEntries) :
			#check the max events
			count+=1
			if count == options.max_events+1 :
				print 'Processed event number '+str(count-1)+', exiting'
				break
			#print progress
			if count % options.print_every == 0 or count == 1:
				print 'Count at '+str(count)+' out of '+str(nEntries)+', (%.4f%% complete)'%(float(count) / float(nEntries) * 100.0)
			chain.GetEntry(entry)
			cuts = []
			cuts.append(self.name.lower().find('data')!=-1 or weight[0]!=1.0)
			#cuts for trigger analysis
			if m == 'trigger' :
				#muon1 selection and trigger
				cuts.append(mu_trigger[0]==1)
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
				cuts.append(ele1_Q[0]+muon1_Q[0]==0)
				#Other leptonic side cuts
				cuts.append(lepW_pt[0]>50.)
				cuts.append(lept_M[0]>140. and lept_M[0]<250.)
				#Hadronic side cuts
				cuts.append(hadt_pt[0]>50. and lepb_pt[0]>50. and max(hadt_pt[0],lepb_pt[0])>150.)
			#cuts for ID analysis
			elif m == 'ID' :
				#electron 1 and/or 2 selection
				ele1_tag_selec = ele1_pt[0]>40. and abs(ele1_eta[0])<2.4 and ele1_isLoose[0]==1 and (ele1_relPt[0]>25. or ele1_dR[0]>0.5)
				ele1_probe_selec = ele1_pt[0]>40. and abs(ele1_eta[0])<2.4 and (ele1_relPt[0]>25. or ele1_dR[0]>0.5)
				ele2_tag_selec = ele2_pt[0]>40. and abs(ele2_eta[0])<2.4 and ele2_isLoose[0]==1 and (ele2_relPt[0]>25. or ele2_dR[0]>0.5)
				ele2_probe_selec = ele2_pt[0]>40. and abs(ele2_eta[0])<2.4 and (ele2_relPt[0]>25. or ele2_dR[0]>0.5)
				cuts.append((ele1_tag_selec and ele2_probe_selec) or (ele2_tag_selec and ele1_probe_selec))
				#muon 1 and 2 rejection
				cuts.append(muon1_pt[0]<40. or abs(muon1_eta[0])>2.4 or muon1_isLoose[0]!=1 or (muon1_relPt[0]<25. and muon1_dR[0]<0.5))
				cuts.append(muon2_pt[0]<40. or abs(muon2_eta[0])>2.4 or muon2_isLoose[0]!=1 or (muon2_relPt[0]<25. and muon2_dR[0]<0.5))
				#require two electrons to have opposite charges
				cuts.append(ele1_Q[0]+ele2_Q[0]==0)
				#make the dielectron fourvector and cut on its mass
				e1_fvec = TLorentzVector()
				e2_fvec = TLorentzVector()
				e1_fvec.SetPtEtaPhiM(ele1_pt[0],ele1_eta[0],ele1_phi[0],ele1_M[0])
				e2_fvec.SetPtEtaPhiM(ele2_pt[0],ele2_eta[0],ele2_phi[0],ele2_M[0])
				ee_fvec = e1_fvec+e2_fvec
				cuts.append(ee_fvec.M()>12.)
				cuts.append(ee_fvec.M()<76. or ee_fvec.M()>106.)
				#other leptonic side cuts
				cuts.append(lepW_pt[0]>50.)
				#cuts.append(lept_M[0]>140. and lept_M[0]<250.)
				#other jet requirements
				cuts.append(hadt_pt[0]>50. and lepb_pt[0]>50. and max(hadt_pt[0],lepb_pt[0])>150.)
				cuts.append(lepb_csv[0]>CSVL or hadt_csv[0]>CSVL)
				#if the tagged electron is ele2, switch some variables around
				if ele2_tag_selec and ele1_probe_selec and not (ele1_tag_selec and ele2_probe_selec) :
					ele2_pt[0] = ele1_pt[0]
					ele2_eta[0] = ele1_eta[0]
					ele2_Q[0] = ele1_Q[0]
					ele2_isLoose[0] = ele1_isLoose[0]

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
			self.tree.Fill()
			#fill the histograms
			# weight should be (245.8/21560109.)*LUMINOSITY*sf_pileup[0] for ttbar MC
			ev_weight = 1.0
			if self.name.lower().find('mc')!=-1 :
				ev_weight = weight[0]*LUMINOSITY*sf_pileup[0]*sf_top_pT[0]
			elpt = ele1_pt[0]
			eleta = ele1_eta[0]
			if m=='ID' :
				elpt = ele2_pt[0]
				eleta = ele2_eta[0]
			self.ele_pt_all.Fill(elpt,ev_weight)
			self.ele_eta_all.Fill(eleta,ev_weight)
			self.jet1_pt_all.Fill(jet1_pt[0],ev_weight)
			self.jet2_pt_all.Fill(jet2_pt[0],ev_weight)
			self.histo_2d_all.Fill(elpt,abs(eleta),ev_weight)
			if (m=='trigger' and pass_trig[0]==1) or (m=='ID' and ele2_isLoose[0]==1) :
				self.ele_pt_pass.Fill(elpt,ev_weight)
				self.ele_eta_pass.Fill(eleta,ev_weight)
				self.jet1_pt_pass.Fill(jet1_pt[0],ev_weight)
				self.jet2_pt_pass.Fill(jet2_pt[0],ev_weight)
				self.histo_2d_pass.Fill(elpt,abs(eleta),ev_weight)
		print 'Done'
		#Make the graph y-values and errors
		for i in range(n_ele_pt_bins) :
			x_value = (ele_pt_bins[i+1]+ele_pt_bins[i])/2
			x_err   = (ele_pt_bins[i+1]-ele_pt_bins[i])/2
			passing_events = self.ele_pt_pass.GetBinContent(self.ele_pt_pass.FindBin(x_value))
			all_events = self.ele_pt_all.GetBinContent(self.ele_pt_all.FindBin(x_value))
			y_value = 0.; y_err = 0.
			if all_events > 0. :
				y_value = passing_events/all_events
				if passing_events>0. :
					y_err = y_value*sqrt(1./passing_events+1./all_events)
			self.ele_pt_gr.SetPoint(i,x_value,y_value)
			self.ele_pt_gr.SetPointError(i,x_err,y_err)
		for i in range(n_eta_bins) :
			x_value = (eta_bins[i+1]+eta_bins[i])/2
			x_err   = (eta_bins[i+1]-eta_bins[i])/2
			passing_events = self.ele_eta_pass.GetBinContent(self.ele_eta_pass.FindBin(x_value))
			all_events = self.ele_eta_all.GetBinContent(self.ele_eta_all.FindBin(x_value))
			y_value = 0.; y_err = 0.
			if all_events > 0. :
				y_value = passing_events/all_events
				if passing_events>0. :
					y_err = y_value*sqrt(1./passing_events+1./all_events)
			self.ele_eta_gr.SetPoint(i,x_value,y_value)
			self.ele_eta_gr.SetPointError(i,x_err,y_err)
		for i in range(n_jet1_pt_bins) :
			x_value = (jet1_pt_bins[i+1]+jet1_pt_bins[i])/2
			x_err   = (jet1_pt_bins[i+1]-jet1_pt_bins[i])/2
			passing_events = self.jet1_pt_pass.GetBinContent(self.jet1_pt_pass.FindBin(x_value))
			all_events = self.jet1_pt_all.GetBinContent(self.jet1_pt_all.FindBin(x_value))
			y_value = 0.; y_err = 0.
			if all_events > 0. :
				y_value = passing_events/all_events
				if passing_events>0. :
					y_err = y_value*sqrt(1./passing_events+1./all_events)
			self.jet1_pt_gr.SetPoint(i,x_value,y_value)
			self.jet1_pt_gr.SetPointError(i,x_err,y_err)
		for i in range(n_jet2_pt_bins) :
			x_value = (jet2_pt_bins[i+1]+jet2_pt_bins[i])/2
			x_err   = (jet2_pt_bins[i+1]-jet2_pt_bins[i])/2
			passing_events = self.jet2_pt_pass.GetBinContent(self.jet2_pt_pass.FindBin(x_value))
			all_events = self.jet2_pt_all.GetBinContent(self.jet2_pt_all.FindBin(x_value))
			y_value = 0.; y_err = 0.
			if all_events > 0. :
				y_value = passing_events/all_events
				if passing_events>0. :
					y_err = y_value*sqrt(1./passing_events+1./all_events)
			self.jet2_pt_gr.SetPoint(i,x_value,y_value)
			self.jet2_pt_gr.SetPointError(i,x_err,y_err)
		#Fit the graphs with constants, then save the fit functions as graphs
		self.ele_pt_const = TF1('ele_pt_const','[0]',ele_pt_bins[0],ele_pt_bins[n_ele_pt_bins])
		self.ele_eta_const = TF1('ele_eta_const','[0]',eta_bins[0],eta_bins[n_eta_bins])
		self.jet1_pt_const = TF1('jet1_pt_const','[0]',jet1_pt_bins[0],jet1_pt_bins[n_jet1_pt_bins])
		self.jet2_pt_const = TF1('jet2_pt_const','[0]',jet2_pt_bins[0],jet2_pt_bins[n_jet2_pt_bins])
		self.ele_pt_gr.Fit('ele_pt_const')
		self.ele_eta_gr.Fit('ele_eta_const')
		self.jet1_pt_gr.Fit('jet1_pt_const')
		self.jet2_pt_gr.Fit('jet2_pt_const')
		self.ele_pt_fit_value = self.ele_pt_const.GetParameter(0)
		self.ele_pt_fit_err   = self.ele_pt_const.GetParError(0)
		self.ele_eta_fit_value = self.ele_eta_const.GetParameter(0)
		self.ele_eta_fit_err   = self.ele_eta_const.GetParError(0)
		self.jet1_pt_fit_value = self.jet1_pt_const.GetParameter(0)
		self.jet1_pt_fit_err   = self.jet1_pt_const.GetParError(0)
		self.jet2_pt_fit_value = self.jet2_pt_const.GetParameter(0)
		self.jet2_pt_fit_err   = self.jet2_pt_const.GetParError(0)
		self.ele_pt_fit_gr     = TGraphErrors(n_ele_pt_bins,ele_pt_x,ele_pt_y,ele_pt_xe,ele_pt_ye)
		self.ele_eta_fit_gr    = TGraphErrors(n_eta_bins,ele_eta_x,ele_eta_y,ele_eta_xe,ele_eta_ye)
		self.jet1_pt_fit_gr     = TGraphErrors(n_jet1_pt_bins,jet1_pt_x,jet1_pt_y,jet1_pt_xe,jet1_pt_ye)
		self.jet2_pt_fit_gr     = TGraphErrors(n_jet2_pt_bins,jet2_pt_x,jet2_pt_y,jet2_pt_xe,jet2_pt_ye)
		for i in range(n_ele_pt_bins) :
			x_value = (ele_pt_bins[i+1]+ele_pt_bins[i])/2
			x_err   = (ele_pt_bins[i+1]-ele_pt_bins[i])/2
			self.ele_pt_fit_gr.SetPoint(i,x_value,self.ele_pt_fit_value)
			self.ele_pt_fit_gr.SetPointError(i,x_err,self.ele_pt_fit_err)
		for i in range(n_eta_bins) :
			x_value = (eta_bins[i+1]+eta_bins[i])/2
			x_err   = (eta_bins[i+1]-eta_bins[i])/2
			self.ele_eta_fit_gr.SetPoint(i,x_value,self.ele_eta_fit_value)
			self.ele_eta_fit_gr.SetPointError(i,x_err,self.ele_eta_fit_err)
		for i in range(n_jet1_pt_bins) :
			x_value = (jet1_pt_bins[i+1]+jet1_pt_bins[i])/2
			x_err   = (jet1_pt_bins[i+1]-jet1_pt_bins[i])/2
			self.jet1_pt_fit_gr.SetPoint(i,x_value,self.jet1_pt_fit_value)
			self.jet1_pt_fit_gr.SetPointError(i,x_err,self.jet1_pt_fit_err)
		for i in range(n_jet2_pt_bins) :
			x_value = (jet2_pt_bins[i+1]+jet2_pt_bins[i])/2
			x_err   = (jet2_pt_bins[i+1]-jet2_pt_bins[i])/2
			self.jet2_pt_fit_gr.SetPoint(i,x_value,self.jet2_pt_fit_value)
			self.jet2_pt_fit_gr.SetPointError(i,x_err,self.jet2_pt_fit_err)
		#Write the tree, histograms, and graphs; close the file
		outputfile.cd()
		self.tree.Write()
		for thing in self.histos_and_graphs :
			thing.Write()	

##########									Set up 									##########

print 'Opening files . . . '
#Build path to input file
input_files_list = './'
input_files_list += options.input
if not '.' in options.input :
	input_files_list += '.txt'
print 'Using input file '+input_files_list+''
#figure out whether we're going for trigger or ID efficiency
m = 'trigger'
if input_files_list.find('ee_skim') != -1 :
	m = 'ID'
#open input file in read only mode 
input_files_list = open(input_files_list,'r')
#Set up the chains and add files
data_chain = TChain('tree')
mc_chain = TChain('tree')
filenamelist = []
for input_file in input_files_list :
	filenamelist.append(input_file.rstrip())
print 'Chaining files: '
for filename in filenamelist :
	if filename.startswith('#') :
		continue
	if filename.find('SingleMu')!=-1 or filename.find('SingleEl')!=-1 :
		print 'Adding file '+filename+' to data chain'
		data_chain.Add(filename)
	else :
		print 'Adding file '+filename+' to MC chain'
		mc_chain.Add(filename)
#Set up output file
filename = 'electron_'
if m=='trigger':
	filename+='trigger_'
elif m=='ID':
	filename+='ID_'
filename+='efficiency.root'
outfile  = TFile(filename,'recreate')
#Set up file group objects and make histos/graphs
data_group = filegroup('data',data_chain,outfile,m)
mc_group = filegroup('mc',mc_chain,outfile,m)
#calculate numbers
data_events_selected = data_group.ele_pt_all.Integral()
data_events_passing  = data_group.ele_pt_pass.Integral()
mc_events_selected = mc_group.ele_pt_all.Integral()
mc_events_passing  = mc_group.ele_pt_pass.Integral()
data_eff = data_events_passing/data_events_selected
data_eff_err = data_eff*sqrt(1./data_events_passing+1./data_events_selected)
mc_eff = mc_events_passing/mc_events_selected
mc_eff_err = mc_eff*sqrt(1./mc_events_passing+1./mc_events_selected)
data_mc_sf = data_eff/mc_eff
data_mc_sf_err = data_mc_sf*sqrt((data_eff_err/data_eff)*(data_eff_err/data_eff)+(mc_eff_err/mc_eff)*(mc_eff_err/mc_eff))
print 'DATA EVENTS SELECTED = '+str(data_events_selected)
print 'DATA EVENTS PASSING = '+str(data_events_passing) 
print 'DATA EFFICIENCY = '+str(data_eff)+' +/- '+str(data_eff_err)
print 'MC EVENTS SELECTED = '+str(mc_events_selected)
print 'MC EVENTS PASSING = '+str(mc_events_passing)
print 'MC EFFICIENCY = '+str(mc_eff)+' +/- '+str(mc_eff_err)
print 'DATA/MC SCALEFACTOR = '+str(data_mc_sf)+' +/- '+str(data_mc_sf_err)
#make pretty plots
canvs = []
el_pt_canv = TCanvas('el_pt_canv','el_pt_canv',1200,900);	canvs.append(el_pt_canv)
el_eta_canv = TCanvas('el_eta_canv','el_eta_canv',1200,900);	canvs.append(el_eta_canv)
jet1_pt_canv = TCanvas('jet1_pt_canv','jet1_pt_canv',1200,900);	canvs.append(jet1_pt_canv)
jet2_pt_canv = TCanvas('jet2_pt_canv','jet2_pt_canv',1200,900);	canvs.append(jet2_pt_canv)
data_graphs = []
data_graphs.append(data_group.ele_pt_gr)
data_graphs.append(data_group.ele_eta_gr)
data_graphs.append(data_group.jet1_pt_gr)
data_graphs.append(data_group.jet2_pt_gr)
mc_graphs = []
mc_graphs.append(mc_group.ele_pt_gr)
mc_graphs.append(mc_group.ele_eta_gr)
mc_graphs.append(mc_group.jet1_pt_gr)
mc_graphs.append(mc_group.jet2_pt_gr)
data_eff_lines = []
data_eff_lines.append(TLine(ele_pt_bins[0],data_eff,ele_pt_bins[n_ele_pt_bins],data_eff))
data_eff_lines.append(TLine(eta_bins[0],data_eff,eta_bins[n_eta_bins],data_eff))
data_eff_lines.append(TLine(jet1_pt_bins[0],data_eff,jet1_pt_bins[n_jet1_pt_bins],data_eff))
data_eff_lines.append(TLine(jet2_pt_bins[0],data_eff,jet2_pt_bins[n_jet2_pt_bins],data_eff))
mc_eff_lines = []
mc_eff_lines.append(TLine(ele_pt_bins[0],mc_eff,ele_pt_bins[n_ele_pt_bins],mc_eff))
mc_eff_lines.append(TLine(eta_bins[0],mc_eff,eta_bins[n_eta_bins],mc_eff))
mc_eff_lines.append(TLine(jet1_pt_bins[0],mc_eff,jet1_pt_bins[n_jet1_pt_bins],mc_eff))
mc_eff_lines.append(TLine(jet2_pt_bins[0],mc_eff,jet2_pt_bins[n_jet2_pt_bins],mc_eff))
mc_fit_graphs = []
mc_fit_graphs.append(mc_group.ele_pt_fit_gr)
mc_fit_graphs.append(mc_group.ele_eta_fit_gr)
mc_fit_graphs.append(mc_group.jet1_pt_fit_gr)
mc_fit_graphs.append(mc_group.jet2_pt_fit_gr)
data_fit_graphs = []
data_fit_graphs.append(data_group.ele_pt_fit_gr)
data_fit_graphs.append(data_group.ele_eta_fit_gr)
data_fit_graphs.append(data_group.jet1_pt_fit_gr)
data_fit_graphs.append(data_group.jet2_pt_fit_gr)
for i in range(len(canvs)) :
	canvs[i].cd()
	data_graphs[i].SetMarkerStyle(21)
	data_graphs[i].SetMarkerColor(kRed)
	data_graphs[i].SetFillStyle(3004)
	data_graphs[i].SetFillColor(kRed)
	mc_graphs[i].SetMarkerStyle(21)
	mc_graphs[i].SetMarkerColor(kBlue)
	mc_graphs[i].SetFillStyle(3005)
	mc_graphs[i].SetFillColor(kBlue)
	data_fit_graphs[i].SetLineWidth(3)
	data_fit_graphs[i].SetLineStyle(2)
	data_fit_graphs[i].SetLineColor(kRed)
	data_fit_graphs[i].SetFillColor(kRed)
	data_fit_graphs[i].SetFillStyle(3003)
	mc_fit_graphs[i].SetLineWidth(3)
	mc_fit_graphs[i].SetLineStyle(2)
	mc_fit_graphs[i].SetLineColor(kBlue)
	mc_fit_graphs[i].SetFillColor(kBlue)
	mc_fit_graphs[i].SetFillStyle(3003)
	data_eff_lines[i].SetLineWidth(3)
	data_eff_lines[i].SetLineColor(kRed)
	mc_eff_lines[i].SetLineWidth(3)
	mc_eff_lines[i].SetLineColor(kBlue)
	if m=='trigger' :
		data_graphs[i].GetYaxis().SetRangeUser(0.5,1.3)
	elif m=='ID' :
		data_graphs[i].GetYaxis().SetRangeUser(0.0,1.0)
	data_graphs[i].SetTitle(data_graphs[i].GetTitle().lstrip(data_group.name))
	leg = TLegend(0.2,0.7,0.4,0.85)
	leg.AddEntry(data_graphs[i],'data','PEF')
	#leg.AddEntry(data_fit_graphs[i],'data fit','LF')
	#leg.AddEntry(data_eff_lines[i],'overall data','L')
	leg.AddEntry(mc_graphs[i],'MC','PEF')
	#leg.AddEntry(mc_fit_graphs[i],'MC fit','LF')
	#leg.AddEntry(mc_eff_lines[i],'overall MC','L')
	data_graphs[i].Draw('APE3')
	mc_graphs[i].Draw('PE3 SAME')
	#data_eff_lines[i].Draw('SAME')
	#mc_eff_lines[i].Draw('SAME')
	#data_fit_graphs[i].Draw('LE3 SAME')
	#mc_fit_graphs[i].Draw('LE3 SAME')
	leg.Draw('SAME')
	outfile.cd()
	canvs[i].Write()
#make the 2D scalefactor plot
sf_plot = TH2D('sf_plot','Data/MC efficiency scalefactors over electron phase space; p_{T} (GeV); |#eta|',n_ele_pt_bins_2D,ele_pt_bins_2D,n_eta_bins_2D,eta_bins_2D)
for i in range(n_ele_pt_bins_2D) :
	for j in range(n_eta_bins_2D) :
		xvalue = (ele_pt_bins_2D[i+1]+ele_pt_bins_2D[i])/2.
		yvalue = (eta_bins_2D[j+1]+eta_bins_2D[j])/2.
		thisbin = mc_group.histo_2d_all.FindFixBin(xvalue,yvalue)
		mcall   = mc_group.histo_2d_all.GetBinContent(thisbin)
		mcpass  = mc_group.histo_2d_pass.GetBinContent(thisbin)
		dataall   = data_group.histo_2d_all.GetBinContent(thisbin)
		datapass  = data_group.histo_2d_pass.GetBinContent(thisbin)
		content = 1.0; err = 1.0
		if datapass!=0. and dataall!=0. and mcpass!=0. and mcall!=0. :
			content = (datapass/dataall)/(mcpass/mcall)
			err = content*sqrt(1./datapass+1./dataall+1./mcpass+1./mcall)
		sf_plot.SetBinContent(thisbin,content)
		sf_plot.SetBinError(thisbin,err)
sf_plot_canv = TCanvas('sf_plot_canv','sf_plot_canv',1200,900);
sf_plot_canv.cd()
sf_plot.Draw("COLZ")
outfile.cd()
sf_plot_canv.Write()
sf_plot.Write()



outfile.Close()