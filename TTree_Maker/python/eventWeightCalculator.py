#eventWeightCalculator contains functions that return scalefactors for a bunch of different corrections
#currently not really implemented because we don't know what these should be for 13TeV data!
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

import ROOT
from math import *
import pickle

#global variables
MUON_ID_EFF_PKL_FILENAME   = 'MuonEfficiencies_Run2012ReReco_53X.pkl'
ELECTRON_ID_EFF_HISTO_FILENAME = 'electrons_scale_factors.root'
ELECTRON_ID_EFF_HISTO_NAME = 'electronsDATAMCratio_FO_ID'
MUON_TRIG_EFF_PKL_FILENAME = 'SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.pkl' 
DATA_PU_FILENAME = 'data_pileup_distribution.root'
MC_PU_FILENAME = 'dumped_Powheg_TT.root'

##########							   MC_corrector Class 								##########

class MC_corrector :
	##################################		#__doc__		##################################
	"""MC corrector class; calculates scalefactors and reweights for an event"""

	##################################  #__init__ function  ##################################
	def __init__(self,generator,eventType,lepType,onGrid) :
		self.generator = generator
		self.event_type = eventType
		self.lep_type = lepType
		self.on_grid = onGrid.lower()
		self.__open_pileup_files__()
		self.__open_lep_eff_files__()

	def __open_pileup_files__(self) :
		prepend = ''
		if self.on_grid == 'yes' :
			prepend+='./tardir/'
		else :
			prepend+='../../other_input_files/'
		#note that these files should definitely be recalculated once I decide on a trigger, etc.
		data_pu_file = ROOT.TFile(prepend+DATA_PU_FILENAME)
		MC_pu_file = ROOT.TFile(prepend+MC_PU_FILENAME)
		self.data_pu_dist = ROOT.TH1D()
		self.MC_pu_dist   = ROOT.TH1D()
		data_dist = data_pu_file.Get('pileup')
		MC_dist = MC_pu_file.Get('pileup')
		self.data_pu_dist = data_dist.Clone()
		self.MC_pu_dist = MC_dist.Clone()
		self.data_pu_dist.SetDirectory(0)
		self.MC_pu_dist.SetDirectory(0)
		self.data_pu_dist.Scale(1.0/self.data_pu_dist.Integral())
		self.MC_pu_dist.Scale(1.0/self.MC_pu_dist.Integral())

	def __open_lep_eff_files__(self) :
		prepend = ''
		if self.on_grid == 'yes' :
			prepend+='./tardir/'
		else :
			prepend+='../../other_input_files/'
		if self.lep_type == 1 :
			print 'WARNING: ELECTRON TRIGGER EFFICIENCY SCALE FACTOR CORRECTIONS NOT YET IMPLEMENTED!'
			self.lepton_id_histo_file = ROOT.TFile(prepend+ELECTRON_ID_EFF_HISTO_FILENAME)
			self.lepton_id_histo = self.lepton_id_histo_file.Get(ELECTRON_ID_EFF_HISTO_NAME)
			return
		elif self.lep_type == 0 :
			self.lepton_id_file   = open(prepend+MUON_ID_EFF_PKL_FILENAME)
			self.lepton_trig_file = open(prepend+MUON_TRIG_EFF_PKL_FILENAME)
			self.lepton_id_dict   = pickle.load(self.lepton_id_file)
			self.lepton_trig_dict = pickle.load(self.lepton_trig_file)

	def __del__(self) :
		if self.lep_type == 0 :
			self.lepton_id_file.close()
			self.lepton_trig_file.close()
		return #DEBUG

	#top pT reweighting, for Powheg and madgraph semileptonic/dileptonic events
	#takes in the MC truth t and tbar vectors
	#returns the reweighting factor
	#8TeV Powheg and Madgraph semi/dilep events: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
	def getToppT_reweight(self,t_vec,tbar_vec) :
		if (self.generator != 'powheg' and self.generator != 'madgraph') or self.event_type>2 : 
			return 1.0
		a = 0.; b = 0.
		if self.event_type<2 :
			a = 0.159; b = -0.00141
		elif self.event_type == 2 :
			a = 0.148; b = -0.00129
		sf_top = exp(a+b*t_vec.Pt()); sf_antitop = exp(a+b*tbar_vec.Pt())
		return sqrt(sf_top*sf_antitop)
		#return 1.0 #DEBUG RETURN
	
	#pileup reweighting factor
	#takes in the generated and reconstructed number of interactions
	#returns the scalefactor
	def getpileup_reweight(self,MCpileup) :
		data_pu = self.data_pu_dist.GetBinContent(self.data_pu_dist.FindFixBin(1.0*MCpileup))
		mc_pu 	= self.MC_pu_dist.GetBinContent(self.MC_pu_dist.FindFixBin(1.0*MCpileup))
		return data_pu/mc_pu
		#return 1.0 #DEBUG RETURN
	
	#lepton ID efficiency scale factor
	#takes in the measured event pileup and lepton pT and eta
	#returns a three-tuple of (scalefactor, sf_low, sf_hi)
	def getID_eff(self,pileup,meas_lep_pt,meas_lep_eta) :
		#muons
		if self.lep_type == 0 :
			vtxsubdict = self.lepton_id_dict['Tight']['vtxpt20-500']
			nextvtxkey = getNextKey(vtxsubdict.keys(),pileup) 
			sf_vtx = vtxsubdict[nextvtxkey]['data/mc']['efficiency_ratio']
			sf_vtx_low = vtxsubdict[nextvtxkey]['data/mc']['err_low']
			sf_vtx_hi  = vtxsubdict[nextvtxkey]['data/mc']['err_hi']
			etasubdict = self.lepton_id_dict['Tight']['etapt20-500']
			nextetakey = getNextKey(etasubdict.keys(),meas_lep_eta)
			sf_eta = etasubdict[nextetakey]['data/mc']['efficiency_ratio']
			sf_eta_low = etasubdict[nextetakey]['data/mc']['err_low']
			sf_eta_hi  = etasubdict[nextetakey]['data/mc']['err_hi']
			if abs(meas_lep_eta) < 0.9 :
				ptsubdict = self.lepton_id_dict['Tight']['ptabseta<0.9']
			elif abs(meas_lep_eta) < 1.2 :
				ptsubdict = self.lepton_id_dict['Tight']['ptabseta0.9-1.2']
			elif abs(meas_lep_eta) < 2.1 :
				ptsubdict = self.lepton_id_dict['Tight']['ptabseta1.2-2.1']
			else :
				ptsubdict = self.lepton_id_dict['Tight']['ptabseta2.1-2.4']
			nextptkey = getNextKey(ptsubdict.keys(),meas_lep_pt)
			sf_pt = ptsubdict[nextptkey]['data/mc']['efficiency_ratio']
			sf_pt_low = ptsubdict[nextptkey]['data/mc']['err_low']
			sf_pt_hi  = ptsubdict[nextptkey]['data/mc']['err_hi']
			sf = sf_vtx*sf_eta*sf_pt
			sf_low = sf*(1.0-sqrt((sf_vtx_low/sf_vtx)**2+(sf_eta_low/sf_eta)**2+(sf_pt_low/sf_pt)**2))
			sf_hi  = sf*(1.0+sqrt((sf_vtx_hi/sf_vtx)**2+(sf_eta_hi/sf_eta)**2+(sf_pt_hi/sf_pt)**2))
		#electrons
		#correction taken right from file listed under "Triggering MVA" in 
		#"Recommended Working Points with Scale Factors" on
		#https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentification
		elif self.lep_type == 1 :
			if meas_lep_pt>200. :
				meas_lep_pt = 200.
			if abs(meas_lep_eta)>2.5 :
				meas_lep_eta = 2.5
			bin = self.lepton_id_histo.FindFixBin(abs(meas_lep_eta),meas_lep_pt)
			sf = self.lepton_id_histo.GetBinContent(bin)
			sf_err = self.lepton_id_histo.GetBinError(bin)
			sf_low = sf-sf_err; sf_hi = sf+sf_err
		return sf,sf_low,sf_hi
		#return 1.0 #DEBUG RETURN
	
	#trigger efficiency scale factor with lepton pT-based trigger
	#takes in the measured event pileup and lepton pT and eta
	#returns a three-tuple of (scalefactor, err_low, err_hi)
	def gettrig_eff(self,pileup,meas_lep_pt,meas_lep_eta) :
		#muons
		if self.lep_type == 0 :
			vtxsubdict = self.lepton_trig_dict['Mu40']['TightID']['VTX']
			nextvtxkey = getNextKey(vtxsubdict.keys(),pileup) 
			sf_vtx = vtxsubdict[nextvtxkey]['data/mc']['efficiency_ratio']
			sf_vtx_low = vtxsubdict[nextvtxkey]['data/mc']['err_low']
			sf_vtx_hi  = vtxsubdict[nextvtxkey]['data/mc']['err_hi']
			etasubdict = self.lepton_trig_dict['Mu40']['TightID']['ETA']
			nextetakey = getNextKey(etasubdict.keys(),meas_lep_eta)
			sf_eta = etasubdict[nextetakey]['data/mc']['efficiency_ratio']
			sf_eta_low = etasubdict[nextetakey]['data/mc']['err_low']
			sf_eta_hi  = etasubdict[nextetakey]['data/mc']['err_hi']
			if abs(meas_lep_eta) < 0.9 :
				ptsubdict = self.lepton_trig_dict['Mu40']['TightID']['PT_ABSETA_Barrel_0to0p9']
			elif abs(meas_lep_eta) < 1.2 :
				ptsubdict = self.lepton_trig_dict['Mu40']['TightID']['PT_ABSETA_Transition_0p9to1p2']
			else :
				ptsubdict = self.lepton_trig_dict['Mu40']['TightID']['PT_ABSETA_Endcaps_1p2to2p1']
			nextptkey = getNextKey(ptsubdict.keys(),meas_lep_pt)
			sf_pt = ptsubdict[nextptkey]['data/mc']['efficiency_ratio']
			sf_pt_low = ptsubdict[nextptkey]['data/mc']['err_low']
			sf_pt_hi  = ptsubdict[nextptkey]['data/mc']['err_hi']
			sf = sf_vtx*sf_eta*sf_pt
			sf_low = sf*(1.0-sqrt((sf_vtx_low/sf_vtx)**2+(sf_eta_low/sf_eta)**2+(sf_pt_low/sf_pt)**2))
			sf_hi  = sf*(1.0+sqrt((sf_vtx_hi/sf_vtx)**2+(sf_eta_hi/sf_eta)**2+(sf_pt_hi/sf_pt)**2))
		#electrons
		elif self.lep_type == 1 :
			#DEBUG RETURN, THIS IS NOT A REAL THING YET!!!!
			sf = 1.0; sf_low = 1.0; sf_hi = 1.0
		return sf,sf_low,sf_hi
		#return 1.0 #DEBUG RETURN

#getNextKey finds the next key to use in dereferencing the dictionary
def getNextKey(keyList,num) :
	if num < float(keyList[0].split('_')[0]) :
		return keyList[0]
	elif num >= float(keyList[len(keyList)-1].split('_')[1]) :
		return keyList[len(keyList)-1]
	else :
		for key in keyList :
			if num >= float(key.split('_')[0]) and num < float(key.split('_')[1]) :
				return key