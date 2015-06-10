
#TTree Maker workhorse code
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym


#Global variables
#Error codes:
ERR_NONE = 0 		   #	0 = no Error
ERR_INVALID_INIT = 1   #	1 = invalid initialization options
ERR_INVALID_HANDLE = 2 #	2 = invalid Handle for event
#Cutflow codes:
CUTFLOW_FULL_EVENT = 0
CUTFLOW_NOT_ENOUGH_JETS = 1
CUTFLOW_NO_LEPTON = 2
CUTFLOW_KINEMATIC_FIT_FAIL = 3
#Beam energy
SQRT_S=13000.0
BEAM_ENERGY=SQRT_S/2.0

##########								   Imports  								##########

import ROOT
from DataFormats.FWLite import Events, Handle
from array import array
from math import *
import sys
from muon import muon
from electron import electron
from jet import jet
from eventTypeHelper import eventTypeCheck, findInitialQuark, findMCTops
from lepHelper import muonCuts, electronCuts
from metHelper import setupMET
from jetHelper import selectJets
from ttbarReconstructor import reconstruct
from angleReconstructor import getObservables, getMCObservables
from eventWeightCalculator import *

##########							   Treemaker Class 								##########

class treemaker :
	##################################		#__doc__		##################################
	"""treemaker class; calculates and outputs all TTree variables for an event"""

	##################################  Handles and Labels  ##################################
	#MC GenEvent info
	genHandle = Handle('vector<reco::GenParticle>'); genLabel  = ('prunedGenParticles','')
	vector_of_4vecs = 'vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >'
	#muons
	muHandles = []; muLabels = []
	muLabels.append(('jhuMuonPFlowLoose','muonLoose')); 	   muHandles.append(Handle(vector_of_4vecs))
	muLabels.append(('jhuMuonPFlowLoose','muonLoosecharge'));  muHandles.append(Handle('vector<int>'))
	muLabels.append(('jhuMuonPFlowLoose','muonLooseistight')); muHandles.append(Handle('vector<unsigned int>'))
	muLabels.append(('jhuMuonPFlowLoose','muonLooseisloose')); muHandles.append(Handle('vector<unsigned int>'))
	#electrons
	elHandles = []; elLabels = []
	elLabels.append(('jhuElePFlowLoose','electronLoose')); 		  elHandles.append(Handle(vector_of_4vecs))
	elLabels.append(('jhuElePFlowLoose','electronLoosecharge'));  elHandles.append(Handle('vector<int>'))
	elLabels.append(('jhuElePFlowLoose','electronLooseistight')); elHandles.append(Handle('vector<unsigned int>'))
	elLabels.append(('jhuElePFlowLoose','electronLooseisloose')); elHandles.append(Handle('vector<unsigned int>')) 
	#MET
	metHandles = []; metLabels = []
	metLabels.append(('jhuGen','metpt'));  metHandles.append(Handle('double'))
	metLabels.append(('jhuGen','metphi')); metHandles.append(Handle('double'))
	#Jets
	jetHandles = [];	jetLabels = []
	jetLabels.append(('jhuCa8pp','PrunedCA8')); 			 jetHandles.append(Handle(vector_of_4vecs))
	jetLabels.append(('jhuCa8','UnprunedCA8')); 			 jetHandles.append(Handle(vector_of_4vecs))
	jetLabels.append(('jhuCa8','UnprunedCA8tau1')); 		 jetHandles.append(Handle('vector<double>'))
	jetLabels.append(('jhuCa8','UnprunedCA8tau2')); 		 jetHandles.append(Handle('vector<double>'))
	jetLabels.append(('jhuCa8','UnprunedCA8tau3')); 		 jetHandles.append(Handle('vector<double>'))
	jetLabels.append(('jhuCa8pp','PrunedCA8csv')); 			 jetHandles.append(Handle('vector<double>'))
	jetLabels.append(('jhuCa8pp','PrunedCA8PartonFlavour')); jetHandles.append(Handle('vector<int>'))
	#pileup
	pileupLabel = ('jhuGen','npv'); 	  pileupHandle 	= Handle('unsigned int')
	MCpileupLabel = ('jhuGen','npvTrue'); MCpileupHandle = Handle('unsigned int')

	##################################  ANALYZE FUNCTION  ##################################
	def analyze(self,event) :
		self.ERR_CODE = ERR_NONE; self.cutflow[0] = CUTFLOW_FULL_EVENT
		#keep track of whether event has been cut
		keepEvent = True
		#event type split
		if self.is_data == 0 :
			#GenParticles
			event.getByLabel(self.genLabel,self.genHandle)
			if not self.genHandle.isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			GenParticles = self.genHandle.product()
			#pythia8 nTuple genParticles
			genPartVars = []
			for i in range(len(self.genPartHandles)) :
				event.getByLabel(self.genPartLabels[i],self.genPartHandles[i])
				if not self.genPartHandles[i].isValid() :
					self.ERR_CODE = ERR_INVALID_HANDLE
					return self.ERR_CODE
				genPartVars.append(self.genPartHandles[i].product())
			if self.event_type != 4 :
				keepEvent,add_twice = eventTypeCheck(self.MC_generator,GenParticles,genPartVars,self.event_type) 
							#abovefunction in eventTypeHelper.py
				if add_twice :
					self.addTwice[0] = 1
		if not keepEvent :
			return self.ERR_CODE
		#Mother particle and MC truth top assignment
		if self.is_data == 0 : #MC truth values only relevant for semileptonic qqbar->ttbar
			q_vec 	 = findInitialQuark(self.MC_generator,GenParticles,genPartVars) #function in eventTypeHelper.py
			qbar_vec = ROOT.TLorentzVector(q_vec.X(),q_vec.Y(),-1.0*q_vec.Z(),q_vec.E())
			MCt_vec, MCtbar_vec = findMCTops(self.MC_generator,GenParticles) #function in eventTypeHelper.py
		else : #if we don't have the MC truth information, we have to assign which is which later when we do the boost
			q_vec 	 = ROOT.TLorentzVector(1.0,0.0,sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY)
			qbar_vec = ROOT.TLorentzVector(1.0,0.0,-1.0*sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY)
			MCt_vec    = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
			MCtbar_vec = ROOT.TLorentzVector(-1.0,0.0,0.0,1.0)
		self.__fillMC__(q_vec,qbar_vec,MCt_vec,MCtbar_vec)

		#get all the info from the event
		#MET
		metVars = []
		for i in range(len(self.metHandles)) :
			event.getByLabel(self.metLabels[i],self.metHandles[i])
			if not self.metHandles[i].isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			metVars.append(self.metHandles[i].product())
		met = ROOT.TLorentzVector(); met.SetPtEtaPhiM(metVars[0][0], 0., metVars[0][1], 0.)
		self.__fillMET__(met)
		#jets
		jets = []; jetVars = []
		for i in range(len(self.jetHandles)) :
			event.getByLabel(self.jetLabels[i],self.jetHandles[i])
			if not self.jetHandles[i].isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			jetVars.append(self.jetHandles[i].product())
		for i in range(len(jetVars[0])) :
			flavor = -1
			if len(jetVars>6) :
				flavor = jetVars[6][i]
			newJet = jet(jetVars[0][i],jetVars[1],jetVars[2],jetVars[3],jetVars[4],jetVars[5][i],flavor)
			jets.append(newJet)
		#separate the jets into the top and b candidates
		jets = selectJets(jets)
		self.__fillJets__(jets)
		if len(jets)<2 :
			return self.__closeout__(CUTFLOW_NOT_ENOUGH_JETS)
		self.sf_btag_eff[0] = jets[1].btagSF 
		self.sf_btag_eff_low[0] = jets[1].btagSFlow 
		self.sf_btag_eff_hi[0] = jets[1].btagSFhigh
		#muons
		muons = []; muVars = []
		for i in range(len(self.muHandles)) :
			event.getByLabel(self.muLabels[i],self.muHandles[i])
			if not self.muHandles[i].isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			muVars.append(self.muHandles[i].product())
		for i in range(len(muVars[0])) :
			newMuon = muon(muVars[0][i],muVars[1][i],muVars[2][i],muVars[3][i],jets)
			muons.append(newMuon)
		muons.sort(key = lambda x: x.vec.Pt())
		self.__fillMuons__(muons)
		if self.lep_type == 0 and len(muons)<1 :
			return self.__closeout__(CUTFLOW_NO_LEPTON)
		#electrons
		electrons = []; elVars = []
		for i in range(len(self.elHandles)) :
			event.getByLabel(self.elLabels[i],self.elHandles[i])
			if not self.elHandles[i].isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			elVars.append(self.elHandles[i].product())
		for i in range(len(elVars[0])) :
			newElectron = electron(elVars[0][i],elVars[1][i],elVars[2][i],elVars[3][i],met,jets)
			electrons.append(newElectron)
		electrons.sort(key = lambda x: x.vec.Pt())
		self.__fillElectrons__(electrons)
		if self.lep_type == 1 and len(electrons)<1 :
			return self.__closeout__(CUTFLOW_NO_LEPTON)
		#pileup
		event.getByLabel(self.pileupLabel,self.pileupHandle)
		if not self.pileupHandle.isValid() :
			self.ERR_CODE = ERR_INVALID_HANDLE
			return self.ERR_CODE
		pileup = self.pileupHandle.product()[0]
		MCpileup = 0
		if self.is_data == 0 :
			event.getByLabel(self.MCpileupLabel,self.MCpileupHandle)
			if not self.MCpileupHandle.isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			MCpileup = self.MCpileupHandle.product()[0]
		self.pileup[0] = int(pileup); self.MC_pileup[0] = int(MCpileup)

#		#lepton selection, put nice leptons ahead of high-pT leptons
#		muons = muonCuts(muons); electrons = electronCuts(electrons)
#		self.__fillMuons__(muons); self.__fillElectrons__(electrons)
		
		#neutrino handling and setup for fit
		met1_vec, met2_vec = setupMET(lep_vec,metVars) #function in metHelper.py
		self.met_pt[0], self.met_eta[0] = met1_vec.Pt(),  met1_vec.Eta() 
		self.met_phi[0], self.met_M[0]  = met1_vec.Phi(), met1_vec.M()
		self.nFits[0] = 2
		if met1_vec.Pz() == met2_vec.Pz() :
			self.nFits[0] = 1
		#fill the rest of the leptonic fourvectors
		if self.lep_type==0 :
			self.__fillLepSide__(muons[0].vec,met1_vec,jets[1].vec)
		elif self.lep_type==1 :
			self.__fillLepSide__(electrons[0].vec,met1_vec,jets[1].vec)

		#event reconstruction with kinematic fit
		scaledlep = ROOT.TLorentzVector(); scaledmet = ROOT.TLorentzVector() 
		scaledlepb = ROOT.TLorentzVector(); scaledhadt = ROOT.TLorentzVector()
		if self.lep_type == 0 :
			scaledlep, scaledmet, scaledlepb, scaledhadt, self.chi2[0] = reconstruct(muons[0].vec,met1_vec,met2_vec,jets) 
		elif self.lep_type == 1 :
			scaledlep, scaledmet, scaledlepb, scaledhadt, self.chi2[0] = reconstruct(electrons[0].vec,met1_vec,met2_vec,jets) 
		#above function in ttbarReconstructor.py

		#fill the TTree with the scaled fourvector variables
		self.__fillScaledFourVecs__(scaledlep,scaledmet,scaledlepb,scaledhadt)

		#reconstruct the observables using both the scaled and unscaled vectors
		if self.lep_type == 0 :
			self.cstar[0], self.x_F[0], self.M[0] = getObservables(muons[0].vec+met1_vec+jets[1].vec,jets[0].vec,self.Q_l[0]) 
		elif self.lep_type == 1 :
			self.cstar[0], self.x_F[0], self.M[0] = getObservables(electrons[0].vec+met1_vec+jets[1].vec,jets[0].vec,self.Q_l[0]) 
		( self.cstar_scaled[0], self.x_F_scaled[0], 
			self.M_scaled[0] ) = getObservables(scaledlep+scaledmet+scaledlepb,scaledhadt,self.Q_l[0]) 
		#above function in angleReconstructor.py

		#MC Truth observable and reweighting calculation
		if self.is_data==0 :
			if self.event_type!=4 :
				( self.cstar_MC[0],self.x_F_MC[0],self.M_MC[0],
					self.w_a[0],self.w_s_xi[0],self.w_a_xi[0],
					self.w_s_delta[0],self.w_a_delta[0],
					self.w_a_opp[0],self.w_s_xi_opp[0],self.w_a_xi_opp[0],
					self.w_s_delta_opp[0],self.w_a_delta_opp[0] ) = getMCObservables(q_vec,qbar_vec,MCt_vec,MCtbar_vec) 
			#scale factor and reweighting calculations
			#8TeV numbers
			self.sf_top_pT[0] = self.corrector.getToppT_reweight(MCt_vec,MCtbar_vec)
			self.sf_pileup[0] = self.corrector.getpileup_reweight(MCpileup)
			( self.sf_lep_ID[0], self.sf_lep_ID_low[0], 
				self.sf_lep_ID_hi[0] ) = self.corrector.getID_eff(pileup,meas_lep_pt,meas_lep_eta)
			( self.sf_trig_eff[0], self.sf_trig_eff_low[0], 
				self.sf_trig_eff_hi[0] ) = self.corrector.gettrig_eff(pileup,meas_lep_pt,meas_lep_eta)
		self.__closeout__(0) #yay! A successful event!

	##################################  #__init__ function  ##################################
	def __init__(self,fileName,isData,generator,eventType,lepType,reweight,onGrid) :
		self.ERR_CODE = ERR_NONE
		#handle input options
		self.__handleInput__(fileName,isData,generator,eventType,lepType,reweight)
		#book TTree
		self.__book__()
		#Set Monte Carlo reweighter
		self.corrector = MC_corrector(self.MC_generator,self.event_type,self.lep_type,onGrid)
	
	##################################   #__handleInput__   ##################################   
	#############################   __init__ helper function   ###############################
	def __handleInput__(self,fileName,isData,generator,eventType,lepType,reweight) :
		self.ERR_CODE = ERR_NONE
		#output file name
		self.file_name = fileName
		#data or MC?
		if isData == 'no' :
			self.is_data = 0
		elif isData == 'yes' :
			self.is_data = 1
		else :
			print 'ERROR: cannot determine if inputted file is a data or MC file!'
			print '	options.data = '+isData+''
			self.ERR_CODE = ERR_INVALID_INIT
		#MC generator?
		gen = generator.lower()
		if gen == 'powheg' or gen == 'madgraph' or gen == 'pythia8' :
			self.MC_generator = gen
		elif gen == 'mg5' :
			self.MC_generator = 'madgraph'
		elif gen == 'none' :
			self.MC_generator = 'none'
			self.jetLabels.pop(); self.jetHandles.pop()
		#event type?
		if eventType == 'qq_semilep' or eventType == 'qqbar_semilep' or eventType == 'qq' or eventType == 'qqbar' :
			print 'only SEMILEPTONIC QQBAR EVENTS will be analyzed from this file'
			self.event_type = 0
		elif eventType == 'gg_semilep' or eventType == 'gg' :
			print 'only SEMILEPTONIC GG (qg,qiqbarj,etc.) EVENTS will be analyzed from this file'
			self.event_type = 1
		elif eventType == 'dilep' or eventType == 'dileptonic' or eventType == 'di' or eventType == 'fulllep' :
			print 'only DILEPTONIC EVENTS will be analyzed from this file'
			self.event_type = 2
		elif eventType == 'had' or eventType == 'hadronic' :
			print 'only HADRONIC EVENTS will be analyzed from this file'
			self.event_type = 3
		elif eventType == 'none' :
			print 'ALL event types will be analyzed from this file'
			self.event_type = 4
		else :
			print 'ERROR: unrecognized event type specification! Cannot run analysis!'
			print '	options.event_type = '+eventType+''
			self.ERR_CODE = ERR_INVALID_INIT
		#lepton type?
		if lepType == 'muons' or lepType == 'muon' or lepType == 'mu' :
			print 'File will be analyzed using MUON selection'
			self.lep_type = 0
		elif (lepType == 'electrons' or lepType == 'electron' or lepType == 'ele') :
			print 'File will be analyzed using ELECTRON selection'
			self.lep_type = 1
		else :
			print 'ERROR: cannot determine if muon or electron analysis is being performed!'
			print '	lepType = '+lepType+''
			self.ERR_CODE = ERR_INVALID_INIT
		#event scaling?
		self.w = reweight

	##################################   __book__ function  ##################################
	#############################     (init helper function)   ###############################
	def __book__(self) :
		#Define output file and tree
		self.f = ROOT.TFile(self.file_name,'recreate')
		self.f.cd()
		self.tree = ROOT.TTree('tree','tree')
		#list of branch variables with initial values
		self.initial_branches = []
		#Add branches to TTree
		#cutflow
		self.cutflow = array('I',[0]); self.addBranch('cutflow',self.cutflow,'i',0)
		#weights (renormalization, scale factors, analysis)
		self.weight    = array('d',[self.w]); self.addBranch('weight',self.weight,'D',self.w)
		self.w_a 	   	   = array('d',[1.0]); self.addBranch('w_a',		  self.w_a,	   	   	 'D',1.0)
		self.w_s_xi    	   = array('d',[1.0]); self.addBranch('w_s_xi',	      self.w_s_xi,   	 'D',1.0)
		self.w_a_xi    	   = array('d',[1.0]); self.addBranch('w_a_xi',	      self.w_a_xi,   	 'D',1.0)
		self.w_s_delta 	   = array('d',[1.0]); self.addBranch('w_s_delta',    self.w_s_delta,	 'D',1.0)
		self.w_a_delta 	   = array('d',[1.0]); self.addBranch('w_a_delta',    self.w_a_delta,	 'D',1.0)
		self.w_a_opp 	   = array('d',[1.0]); self.addBranch('w_a_opp',	  self.w_a_opp,	 	 'D',1.0)
		self.w_s_xi_opp    = array('d',[1.0]); self.addBranch('w_s_xi_opp',	  self.w_s_xi_opp,   'D',1.0)
		self.w_a_xi_opp    = array('d',[1.0]); self.addBranch('w_a_xi_opp',	  self.w_a_xi_opp,   'D',1.0)
		self.w_s_delta_opp = array('d',[1.0]); self.addBranch('w_s_delta_opp',self.w_s_delta_opp,'D',1.0)
		self.w_a_delta_opp = array('d',[1.0]); self.addBranch('w_a_delta_opp',self.w_a_delta_opp,'D',1.0)
		self.sf_top_pT 		 = array('d',[1.0]); self.addBranch('sf_top_pT', 	   self.sf_top_pT, 		 'D',1.0)
		self.sf_btag_eff 	 = array('d',[1.0]); self.addBranch('sf_btag_eff', 	   self.sf_btag_eff, 	 'D',1.0)
		self.sf_btag_eff_low = array('d',[1.0]); self.addBranch('sf_btag_eff_low', self.sf_btag_eff_low, 'D',1.0)
		self.sf_btag_eff_hi  = array('d',[1.0]); self.addBranch('sf_btag_eff_hi',  self.sf_btag_eff_hi,  'D',1.0)
		self.sf_pileup 		 = array('d',[1.0]); self.addBranch('sf_pileup', 	   self.sf_pileup, 		 'D',1.0)
		self.sf_lep_ID 		 = array('d',[1.0]); self.addBranch('sf_lep_ID', 	   self.sf_lep_ID, 		 'D',1.0)
		self.sf_lep_ID_low 	 = array('d',[1.0]); self.addBranch('sf_lep_ID_low',   self.sf_lep_ID_low, 	 'D',1.0)
		self.sf_lep_ID_hi 	 = array('d',[1.0]); self.addBranch('sf_lep_ID_hi',    self.sf_lep_ID_hi, 	 'D',1.0)
		self.sf_trig_eff 	 = array('d',[1.0]); self.addBranch('sf_trig_eff', 	   self.sf_trig_eff, 	 'D',1.0)
		self.sf_trig_eff_low = array('d',[1.0]); self.addBranch('sf_trig_eff_low', self.sf_trig_eff_low, 'D',1.0)
		self.sf_trig_eff_hi  = array('d',[1.0]); self.addBranch('sf_trig_eff_hi',  self.sf_trig_eff_hi,  'D',1.0)
		#error codes
		self.error_code = array('I',[self.ERR_CODE]); self.addBranch('error_code',self.error_code,'i',self.ERR_CODE)
		#lepton charge
		self.Q_l = array('i',[0]); self.addBranch('Q_l',self.Q_l,'I',0)
		#kinematic fit chi2
		self.chi2 = array('d',[0.0]); self.addBranch('chi2',self.chi2,'D',0.0)
		#number of kinematic fits performed
		self.nFits = array('I',[0]); self.addBranch('nFits',self.nFits,'i',0)
		#whether or not this event should be added twice and have its weight halved based on whether its initial state
		#was symmetric (this will only be nonzero for qqbar and some gg events)
		self.addTwice = array('I',[2]); self.addBranch('addTwice',self.addTwice,'i',2)
		#initial quark vector
		self.q_pt  = array('d',[-1.0]);  self.addBranch('q_pt', self.q_pt, 'D',-1.0)
		self.q_eta = array('d',[100.0]); self.addBranch('q_eta',self.q_eta,'D',100.0)
		self.q_phi = array('d',[100.0]); self.addBranch('q_phi',self.q_phi,'D',100.0)
		self.q_M   = array('d',[-1.0]);  self.addBranch('q_M',  self.q_M,  'D',-1.0)
		#initial antiquark vector
		self.qbar_pt  = array('d',[-1.0]);  self.addBranch('qbar_pt', self.qbar_pt, 'D',-1.0)
		self.qbar_eta = array('d',[100.0]); self.addBranch('qbar_eta',self.qbar_eta,'D',100.0)
		self.qbar_phi = array('d',[100.0]); self.addBranch('qbar_phi',self.qbar_phi,'D',100.0)
		self.qbar_M   = array('d',[-1.0]);  self.addBranch('qbar_M',  self.qbar_M,  'D',-1.0)
		#MC top vector
		self.MCt_pt  = array('d',[-1.0]);  self.addBranch('MCt_pt', self.MCt_pt, 'D',-1.0)
		self.MCt_eta = array('d',[100.0]); self.addBranch('MCt_eta',self.MCt_eta,'D',100.0)
		self.MCt_phi = array('d',[100.0]); self.addBranch('MCt_phi',self.MCt_phi,'D',100.0)
		self.MCt_M   = array('d',[-1.0]);  self.addBranch('MCt_M',  self.MCt_M,  'D',-1.0)
		#MC antitop vector
		self.MCtbar_pt  = array('d',[-1.0]);  self.addBranch('MCtbar_pt', self.MCtbar_pt, 'D',-1.0)
		self.MCtbar_eta = array('d',[100.0]); self.addBranch('MCtbar_eta',self.MCtbar_eta,'D',100.0)
		self.MCtbar_phi = array('d',[100.0]); self.addBranch('MCtbar_phi',self.MCtbar_phi,'D',100.0)
		self.MCtbar_M   = array('d',[-1.0]);  self.addBranch('MCtbar_M',  self.MCtbar_M,  'D',-1.0)
		#muon 1
		self.muon1_pt  = array('d',[-1.0]);  self.addBranch('muon1_pt', self.muon1_pt, 'D',-1.0)
		self.muon1_eta = array('d',[100.0]); self.addBranch('muon1_eta',self.muon1_eta,'D',100.0)
		self.muon1_phi = array('d',[100.0]); self.addBranch('muon1_phi',self.muon1_phi,'D',100.0)
		self.muon1_M   = array('d',[-1.0]);  self.addBranch('muon1_M',  self.muon1_M,  'D',-1.0)
		self.muon1_Q 	   = array('i',[0]); 	self.addBranch('muon1_Q',self.muon1_Q,'I',0)
		self.muon1_isTight = array('I',[2]); 	self.addBranch('muon1_isTight',self.muon1_isTight,'i',2)
		self.muon1_isLoose = array('I',[2]); 	self.addBranch('muon1_isLoose',self.muon1_isLoose,'i',2)
		self.muon1_relPt   = array('d',[-1.0]); self.addBranch('muon1_relPt',self.muon1_relPt,'D',-1.0)
		self.muon1_dR 	   = array('d',[-1.0]); self.addBranch('muon1_dR',self.muon1_relPt,'D',-1.0)
		#muon 2
		self.muon2_pt  = array('d',[-1.0]);  self.addBranch('muon2_pt', self.muon2_pt, 'D',-1.0)
		self.muon2_eta = array('d',[100.0]); self.addBranch('muon2_eta',self.muon2_eta,'D',100.0)
		self.muon2_phi = array('d',[100.0]); self.addBranch('muon2_phi',self.muon2_phi,'D',100.0)
		self.muon2_M   = array('d',[-1.0]);  self.addBranch('muon2_M',  self.muon2_M,  'D',-1.0)
		self.muon2_Q 	   = array('i',[0]); 	self.addBranch('muon2_Q',self.muon2_Q,'I',0)
		self.muon2_isTight = array('I',[2]); 	self.addBranch('muon2_isTight',self.muon2_isTight,'i',2)
		self.muon2_isLoose = array('I',[2]); 	self.addBranch('muon2_isLoose',self.muon2_isLoose,'i',2)
		self.muon2_relPt   = array('d',[-1.0]); self.addBranch('muon2_relPt',self.muon2_relPt,'D',-1.0)
		self.muon2_dR 	   = array('d',[-1.0]); self.addBranch('muon2_dR',self.muon2_relPt,'D',-1.0)
		#electron 1
		self.ele1_pt  = array('d',[-1.0]);  self.addBranch('ele1_pt', self.ele1_pt, 'D',-1.0)
		self.ele1_eta = array('d',[100.0]); self.addBranch('ele1_eta',self.ele1_eta,'D',100.0)
		self.ele1_phi = array('d',[100.0]); self.addBranch('ele1_phi',self.ele1_phi,'D',100.0)
		self.ele1_M   = array('d',[-1.0]);  self.addBranch('ele1_M',  self.ele1_M,  'D',-1.0)
		self.ele1_Q 	   	  = array('i',[0]); 	 self.addBranch('ele1_Q',self.ele1_Q,'I',0)
		self.ele1_isTight 	  = array('I',[2]); 	 self.addBranch('ele1_isTight',self.ele1_isTight,'i',2)
		self.ele1_isLoose 	  = array('I',[2]); 	 self.addBranch('ele1_isLoose',self.ele1_isLoose,'i',2)
		self.ele1_relPt 	  = array('d',[-1.0]);   self.addBranch('ele1_relPt',self.ele1_relPt,'D',-1.0)
		self.ele1_dR 		  = array('d',[-1.0]);   self.addBranch('ele1_dR',self.ele1_relPt,'D',-1.0)
		self.ele1_tri_el_val  = array('d',[-100.0]); self.addBranch('ele1_dR',self.ele1_relPt,'D',-100.0)
		self.ele1_tri_jet_val = array('d',[-100.0]); self.addBranch('ele1_dR',self.ele1_relPt,'D',-100.0)
		self.ele1_tri_cut_val = array('d',[-100.0]); self.addBranch('ele1_dR',self.ele1_relPt,'D',-100.0)
		#electron 2
		self.ele2_pt  = array('d',[-1.0]);  self.addBranch('ele2_pt', self.ele2_pt, 'D',-1.0)
		self.ele2_eta = array('d',[100.0]); self.addBranch('ele2_eta',self.ele2_eta,'D',100.0)
		self.ele2_phi = array('d',[100.0]); self.addBranch('ele2_phi',self.ele2_phi,'D',100.0)
		self.ele2_M   = array('d',[-1.0]);  self.addBranch('ele2_M',  self.ele2_M,  'D',-1.0)
		self.ele2_Q 	   	  = array('i',[0]); 	 self.addBranch('ele2_Q',self.ele2_Q,'I',0)
		self.ele2_isTight 	  = array('I',[2]); 	 self.addBranch('ele2_isTight',self.ele2_isTight,'i',2)
		self.ele2_isLoose 	  = array('I',[2]); 	 self.addBranch('ele2_isLoose',self.ele2_isLoose,'i',2)
		self.ele2_relPt 	  = array('d',[-1.0]);   self.addBranch('ele2_relPt',self.ele2_relPt,'D',-1.0)
		self.ele2_dR 		  = array('d',[-1.0]);   self.addBranch('ele2_dR',self.ele2_relPt,'D',-1.0)
		self.ele2_tri_el_val  = array('d',[-100.0]); self.addBranch('ele2_dR',self.ele2_relPt,'D',-100.0)
		self.ele2_tri_jet_val = array('d',[-100.0]); self.addBranch('ele2_dR',self.ele2_relPt,'D',-100.0)
		self.ele2_tri_cut_val = array('d',[-100.0]); self.addBranch('ele2_dR',self.ele2_relPt,'D',-100.0)
		#neutrino
		self.met_pt  = array('d',[-1.0]);  self.addBranch('met_pt', self.met_pt, 'D',-1.0)
		self.met_eta = array('d',[100.0]); self.addBranch('met_eta',self.met_eta,'D',100.0)
		self.met_phi = array('d',[100.0]); self.addBranch('met_phi',self.met_phi,'D',100.0)
		self.met_M   = array('d',[-1.0]);  self.addBranch('met_M',  self.met_M,  'D',-1.0)
		#leptonic W
		self.lepW_pt  = array('d',[-1.0]);  self.addBranch('lepW_pt', self.lepW_pt, 'D',-1.0)
		self.lepW_eta = array('d',[100.0]); self.addBranch('lepW_eta',self.lepW_eta,'D',100.0)
		self.lepW_phi = array('d',[100.0]); self.addBranch('lepW_phi',self.lepW_phi,'D',100.0)
		self.lepW_M   = array('d',[-1.0]);  self.addBranch('lepW_M',  self.lepW_M,  'D',-1.0)
		#leptonic b
		self.lepb_pt  = array('d',[-1.0]);  self.addBranch('lepb_pt', self.lepb_pt, 'D',-1.0)
		self.lepb_eta = array('d',[100.0]); self.addBranch('lepb_eta',self.lepb_eta,'D',100.0)
		self.lepb_phi = array('d',[100.0]); self.addBranch('lepb_phi',self.lepb_phi,'D',100.0)
		self.lepb_M   = array('d',[-1.0]);  self.addBranch('lepb_M',  self.lepb_M,  'D',-1.0)
		self.lepb_tau32  = array('d',[-1.0]); self.addBranch('lepb_tau32',  self.lepb_tau32,  'D',-1.0)
		self.lepb_tau21  = array('d',[-1.0]); self.addBranch('lepb_tau21',  self.lepb_tau21,  'D',-1.0)
		self.lepb_csv 	 = array('d',[-1.0]); self.addBranch('lepb_csv',  self.lepb_csv,  'D',-1.0)
		self.lepb_flavor = array('i',[0]); 	  self.addBranch('lepb_flavor',self.lepb_flavor,'I',0)
		#leptonic top
		self.lept_pt  = array('d',[-1.0]);  self.addBranch('lept_pt', self.lept_pt, 'D',-1.0)
		self.lept_eta = array('d',[100.0]); self.addBranch('lept_eta',self.lept_eta,'D',100.0)
		self.lept_phi = array('d',[100.0]); self.addBranch('lept_phi',self.lept_phi,'D',100.0)
		self.lept_M   = array('d',[-1.0]);  self.addBranch('lept_M',  self.lept_M,  'D',-1.0)
		#hadronic top
		self.hadt_pt  = array('d',[-1.0]);  self.addBranch('hadt_pt', self.hadt_pt, 'D',-1.0)
		self.hadt_eta = array('d',[100.0]); self.addBranch('hadt_eta',self.hadt_eta,'D',100.0)
		self.hadt_phi = array('d',[100.0]); self.addBranch('hadt_phi',self.hadt_phi,'D',100.0)
		self.hadt_M   = array('d',[-1.0]);  self.addBranch('hadt_M',  self.hadt_M,  'D',-1.0)
		self.hadt_tau32  = array('d',[-1.0]); self.addBranch('hadt_tau32',  self.hadt_tau32,  'D',-1.0)
		self.hadt_tau21  = array('d',[-1.0]); self.addBranch('hadt_tau21',  self.hadt_tau21,  'D',-1.0)
		self.hadt_csv 	 = array('d',[-1.0]); self.addBranch('hadt_csv',  self.hadt_csv,  'D',-1.0)
		self.hadt_flavor = array('i',[0]); 	  self.addBranch('hadt_flavor',self.hadt_flavor,'I',0)
		#rescaled fourvectors
		self.scaled_lep_pt  = array('d',[-1.0]);  self.addBranch('scaled_lep_pt', self.scaled_lep_pt, 'D',-1.0)
		self.scaled_lep_eta = array('d',[100.0]); self.addBranch('scaled_lep_eta',self.scaled_lep_eta,'D',100.0)
		self.scaled_lep_phi = array('d',[100.0]); self.addBranch('scaled_lep_phi',self.scaled_lep_phi,'D',100.0)
		self.scaled_lep_M   = array('d',[-1.0]);  self.addBranch('scaled_lep_M',  self.scaled_lep_M,  'D',-1.0)
		self.scaled_met_pt  = array('d',[-1.0]);  self.addBranch('scaled_met_pt', self.scaled_met_pt, 'D',-1.0)
		self.scaled_met_eta = array('d',[100.0]); self.addBranch('scaled_met_eta',self.scaled_met_eta,'D',100.0)
		self.scaled_met_phi = array('d',[100.0]); self.addBranch('scaled_met_phi',self.scaled_met_phi,'D',100.0)
		self.scaled_met_M   = array('d',[-1.0]);  self.addBranch('scaled_met_M',  self.scaled_met_M,  'D',-1.0)
		self.scaled_lepW_pt  = array('d',[-1.0]);  self.addBranch('scaled_lepW_pt', self.scaled_lepW_pt, 'D',-1.0)
		self.scaled_lepW_eta = array('d',[100.0]); self.addBranch('scaled_lepW_eta',self.scaled_lepW_eta,'D',100.0)
		self.scaled_lepW_phi = array('d',[100.0]); self.addBranch('scaled_lepW_phi',self.scaled_lepW_phi,'D',100.0)
		self.scaled_lepW_M   = array('d',[-1.0]);  self.addBranch('scaled_lepW_M',  self.scaled_lepW_M,  'D',-1.0)
		self.scaled_lepb_pt  = array('d',[-1.0]);  self.addBranch('scaled_lepb_pt', self.scaled_lepb_pt, 'D',-1.0)
		self.scaled_lepb_eta = array('d',[100.0]); self.addBranch('scaled_lepb_eta',self.scaled_lepb_eta,'D',100.0)
		self.scaled_lepb_phi = array('d',[100.0]); self.addBranch('scaled_lepb_phi',self.scaled_lepb_phi,'D',100.0)
		self.scaled_lepb_M   = array('d',[-1.0]);  self.addBranch('scaled_lepb_M',  self.scaled_lepb_M,  'D',-1.0)
		self.scaled_lept_pt  = array('d',[-1.0]);  self.addBranch('scaled_lept_pt', self.scaled_lept_pt, 'D',-1.0)
		self.scaled_lept_eta = array('d',[100.0]); self.addBranch('scaled_lept_eta',self.scaled_lept_eta,'D',100.0)
		self.scaled_lept_phi = array('d',[100.0]); self.addBranch('scaled_lept_phi',self.scaled_lept_phi,'D',100.0)
		self.scaled_lept_M   = array('d',[-1.0]);  self.addBranch('scaled_lept_M',  self.scaled_lept_M,  'D',-1.0)
		self.scaled_hadt_pt  = array('d',[-1.0]);  self.addBranch('scaled_hadt_pt', self.scaled_hadt_pt, 'D',-1.0)
		self.scaled_hadt_eta = array('d',[100.0]); self.addBranch('scaled_hadt_eta',self.scaled_hadt_eta,'D',100.0)
		self.scaled_hadt_phi = array('d',[100.0]); self.addBranch('scaled_hadt_phi',self.scaled_hadt_phi,'D',100.0)
		self.scaled_hadt_M   = array('d',[-1.0]);  self.addBranch('scaled_hadt_M',  self.scaled_hadt_M,  'D',-1.0)
		#cosine(theta)
		self.cstar 		  = array('d',[100.0]); self.addBranch('cstar', 	  self.cstar, 		'D',100.0)
		self.cstar_scaled = array('d',[100.0]); self.addBranch('cstar_scaled',self.cstar_scaled,'D',100.0)
		self.cstar_MC 	  = array('d',[100.0]); self.addBranch('cstar_MC',	  self.cstar_MC, 	'D',100.0)
		#Feynman x
		self.x_F 		= array('d',[100.0]); self.addBranch('x_F', 	  self.x_F,   	  'D',100.0)
		self.x_F_scaled = array('d',[100.0]); self.addBranch('x_F_scaled',self.x_F_scaled,'D',100.0)
		self.x_F_MC 	= array('d',[100.0]); self.addBranch('x_F_MC', 	  self.x_F_MC,	  'D',100.0)
		#ttbar invariant mass
		self.M 		  = array('d',[-1.0]); self.addBranch('M', 		 self.M, 	   'D',-1.0)
		self.M_scaled = array('d',[-1.0]); self.addBranch('M_scaled',self.M_scaled,'D',-1.0)
		self.M_MC 	  = array('d',[-1.0]); self.addBranch('M_MC', 	 self.M_MC,    'D',-1.0)
		#pileup
		self.pileup = array('i',[0]); self.addBranch('pileup',self.pileup,'I',0)
		self.MC_pileup = array('i',[0]); self.addBranch('MC_pileup',self.MC_pileup,'I',0)

	##################################   reset function   ##################################
	#########  sets all relevant values back to zero to get ready for next event  ##########
	def reset(self,err) :
		if err != ERR_NONE :
			self.ERR_CODE = err
		for branch in self.initial_branches :
			branch[0][0] = branch[1]

	################################   addBranch function  #################################
	def addBranch(self,name,var,vartype,ini_val) :
		self.tree.Branch(name,var,name+'/'+vartype)
		self.initial_branches.append((var,ini_val))

	##############################   tree filling functions  ###############################		
	def __fillMC__(qvec,qbarvec,MCtvec,MCtbarvec) :
		self.q_pt[0], 		self.q_eta[0] 	   = qvec.Pt(), 	  qvec.Eta()
		self.q_phi[0], 		self.q_M[0] 	   = qvec.Phi(), 	  qvec.M()
		self.qbar_pt[0], 	self.qbar_eta[0]   = qbarvec.Pt(),    qbarvec.Eta()
		self.qbar_phi[0], 	self.qbar_M[0] 	   = qbarvec.Phi(),   qbarvec.M()
		self.MCt_pt[0], 	self.MCt_eta[0]    = MCtvec.Pt(), 	  MCtvec.Eta()
		self.MCt_phi[0], 	self.MCt_M[0] 	   = MCtvec.Phi(),    MCtvec.M()
		self.MCtbar_pt[0], 	self.MCtbar_eta[0] = MCtbarvec.Pt(),  MCtbarvec.Eta()
		self.MCtbar_phi[0], self.MCtbar_M[0]   = MCtbarvec.Phi(), MCtbarvec.M()

	def __fillMET__(self,metvec) :
		self.met_pt[0]  = metvec.Pt(); self.met_eta[0] = metvec.Eta()
		self.met_phi[0] = metvec.Phi(); self.met_M[0] 	= metvec.M()

	def __fillJets__(self,jetList) :
		if len(jetList) > 0 :
			self.hadt_pt  = jetList[0].vec.Pt(); self.hadt_eta = jetList[0].vec.Eta()
			self.hadt_phi = jetList[0].vec.Phi(); self.hadt_M   = jetList[0].vec.M()
			self.hadt_tau32  = jetList[0].tau32; self.hadt_tau21  = jetList[0].tau21
			self.hadt_csv 	 = jetList[0].csv; self.hadt_flavor = jetList[0].flavor
		if len(jetList) > 1 :
			self.lepb_pt  = jetList[1].vec.Pt(); self.lepb_eta = jetList[1].vec.Eta()
			self.lepb_phi = jetList[1].vec.Phi(); self.lepb_M   = jetList[1].vec.M()
			self.lepb_tau32  = jetList[1].tau32; self.lepb_tau21  = jetList[1].tau21
			self.lepb_csv 	 = jetList[1].csv; self.lepb_flavor = jetList[1].flavor

	def __fillMuons__(self,mulist) :
		if len(mulist) > 0 :
			self.muon1_pt  = mulist[0].vec.Pt(); self.muon1_eta = mulist[0].vec.Eta()
			self.muon1_phi = mulist[0].vec.Phi(); self.muon1_M   = mulist[0].vec.M()
			self.muon1_Q 	   = mulist[0].charge; self.muon1_isTight = mulist[0].isTight
			self.muon1_isLoose = mulist[0].isLoose; self.muon1_relPt   = mulist[0].relPt
			self.muon1_dR 	   = mulist[0].dR
			if self.lep_type == 0 :
				self.Q_l[0] = mulist[0].charge
		if len(mulist) > 1 :
			self.muon2_pt  = mulist[1].vec.Pt(); self.muon2_eta = mulist[1].vec.Eta()
			self.muon2_phi = mulist[1].vec.Phi(); self.muon2_M   = mulist[1].vec.M()
			self.muon2_Q 	   = mulist[1].charge; self.muon2_isTight = mulist[1].isTight
			self.muon2_isLoose = mulist[1].isLoose; self.muon2_relPt   = mulist[1].relPt
			self.muon2_dR 	   = mulist[1].dR

	def __fillElectrons__(self,ellist) :
		if len(ellist) > 0 :
			self.ele1_pt  = ellist[0].vec.Pt(); self.ele1_eta = ellist[0].vec.Eta()
			self.ele1_phi = ellist[0].vec.Phi(); self.ele1_M   = ellist[0].vec.M()
			self.ele1_Q 	   = ellist[0].charge; self.ele1_isTight = ellist[0].isTight
			self.ele1_isLoose = ellist[0].isLoose; self.ele1_relPt   = ellist[0].relPt
			self.ele1_dR 	   = ellist[0].dR
			self.ele1_tri_el_val  = ellist[0].triangle_el_val
			self.ele1_tri_jet_val = ellist[0].triangle_jet_val
			self.ele1_tri_cut_val = ellist[0].triangle_cut_val
			if self.lep_type == 1 :
				self.Q_l[0] = ellist[0].charge
		if len(ellist) > 1 :
			self.ele2_pt  = ellist[1].vec.Pt(); self.ele2_eta = ellist[1].vec.Eta()
			self.ele2_phi = ellist[1].vec.Phi(); self.ele2_M   = ellist[1].vec.M()
			self.ele2_Q 	   = ellist[1].charge; self.ele2_isTight = ellist[1].isTight
			self.ele2_isLoose = ellist[1].isLoose; self.ele2_relPt   = ellist[1].relPt
			self.ele2_dR 	   = ellist[1].dR
			self.ele2_tri_el_val  = ellist[1].triangle_el_val
			self.ele2_tri_jet_val = ellist[1].triangle_jet_val
			self.ele2_tri_cut_val = ellist[1].triangle_cut_val

	def __fillLepSide__(self,lepton,met,lepb) :
		lepW = lepton+met
		self.lepW_pt[0] 	= lepW.Pt(); 	self.lepW_eta[0] = lepW.Eta()
		self.lepW_phi[0] = lepW.Phi(); 	self.lepW_M[0] 	= lepW.M()
		lept = lepW+lepb
		self.lept_pt[0] 	= lept.Pt(); 	self.lept_eta[0] = lept.Eta()
		self.lept_phi[0] = lept.Phi(); 	self.lept_M[0] 	= lept.M()

	def __fillScaledFourVecs__(self,lepton,met,lepb,hadt) :
		self.scaled_lep_pt[0] 	= lepton.Pt(); 	self.scaled_lep_eta[0] = lepton.Eta()
		self.scaled_lep_phi[0] = lepton.Phi(); self.scaled_lep_M[0] 	= lepton.M()
		self.scaled_met_pt[0] 	= met.Pt(); 	self.scaled_met_eta[0] = met.Eta()
		self.scaled_met_phi[0] = met.Phi(); 	self.scaled_met_M[0] 	= met.M()
		lepW = lepton+met
		self.scaled_lepW_pt[0] 	= lepW.Pt(); 	self.scaled_lepW_eta[0] = lepW.Eta()
		self.scaled_lepW_phi[0] = lepW.Phi(); 	self.scaled_lepW_M[0] 	= lepW.M()
		self.scaled_lepb_pt[0]  = lepb.Pt();  self.scaled_lepb_eta[0] = lepb.Eta()
		self.scaled_lepb_phi[0] = lepb.Phi(); self.scaled_lepb_M[0] 	= lepb.M()
		lept = lepW+lepb
		self.scaled_lept_pt[0] 	= lept.Pt(); 	self.scaled_lept_eta[0] = lept.Eta()
		self.scaled_lept_phi[0] = lept.Phi(); 	self.scaled_lept_M[0] 	= lept.M()
		self.scaled_hadt_pt[0]  = hadt.Pt();  self.scaled_hadt_eta[0] = hadt.Eta()
		self.scaled_hadt_phi[0] = hadt.Phi(); self.scaled_hadt_M[0] 	= hadt.M()

	########## function to close out the event, called before kicking back to runner #########
	def __closeout__(self,cut_flow) :
		#update cutflow
		self.cutflow[0] = cut_flow
		#update error code
		self.error_code[0] = int(self.ERR_CODE)
		#fill ttree
		self.tree.Fill()
		#return error code
		return self.ERR_CODE

	################################## __del__ function  ###################################
	def __del__(self) :
		self.f.cd()
		self.f.Write()
		self.f.Close()