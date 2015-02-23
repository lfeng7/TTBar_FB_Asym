
#TTree Maker workhorse code for use with the DFMO framework instead of the B2GAnaFW
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym


#Global variables
#Error codes:
ERR_NONE = 0 		   #	0 = no Error
ERR_INVALID_INIT = 1   #	1 = invalid initialization options
ERR_INVALID_HANDLE = 2 #	2 = invalid Handle for event
#Beam energy
SQRT_S=13000.0
BEAM_ENERGY=SQRT_S/2.0

##########								   Imports  								##########

import ROOT
from DataFormats.FWLite import Events, Handle
from array import array
from math import *
import sys
from eventTypeHelper import eventTypeCheck, findInitialQuark, findMCTops
from lepHelper import leptonCuts, getLeptonFourVec#, leptonCleaningKludge
from metHelper import metCut, setupMET
from jetHelper import selectJets
from ttbarReconstructor import reconstruct
from angleReconstructor import getObservables, getMCObservables
from eventWeightCalculator import *

##########							   Treemaker Class 								##########

class treemaker_DFMO :
	##################################		#__doc__		##################################
	"""treemaker class; calculates and outputs all TTree variables for an event"""

	##################################  13TeV Handles and Labels  ##################################
	#MC GenEvent info
	genHandle = Handle('vector<reco::GenParticle>'); genLabel  = ('prunedGenParticles','')
	#muons
	vector_of_4vecs = 'vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >'
	muHandles = []; muLabels = []
	muLabels.append(('jhuMuonPFlowLoose','muonLoose')); 	   muHandles.append(Handle(vector_of_4vecs))
	muLabels.append(('muons','muEta')); 					   muHandles.append(Handle('vector<float>')) #dummies
	muLabels.append(('muons','muPhi')); 					   muHandles.append(Handle('vector<float>')) #dummies
	muLabels.append(('muons','muMass')); 					   muHandles.append(Handle('vector<float>')) #dummies
	muLabels.append(('jhuMuonPFlowLoose','muonLoosecharge'));  muHandles.append(Handle('vector<int>'))
	muLabels.append(('muons','muSumChargedHadronPt')); 		   muHandles.append(Handle('vector<float>')) #dummies
	muLabels.append(('muons','muSumNeutralHadronPt')); 		   muHandles.append(Handle('vector<float>')) #dummies
	muLabels.append(('muons','muSumPhotonPt')); 			   muHandles.append(Handle('vector<float>')) #dummies
	muLabels.append(('muons','muSumPUPt')); 				   muHandles.append(Handle('vector<float>')) #dummies
	muLabels.append(('jhuMuonPFlowLoose','muonLooseistight')); muHandles.append(Handle('vector<unsigned int>'))
	muLabels.append(('jhuMuonPFlowLoose','muonLooseisloose')); muHandles.append(Handle('vector<unsigned int>'))
	muLabels.append(('muons','muGenMuonE')); 				   muHandles.append(Handle('vector<float>')) #dummies
	#electrons
	elHandles = []; elLabels = []
	elLabels.append(('jhuElePFlowLoose','electronLoose')); 		  elHandles.append(Handle(vector_of_4vecs))
	elLabels.append(('electrons','elEta')); 					  elHandles.append(Handle('vector<float>')) #dummies
	elLabels.append(('electrons','elPhi')); 					  elHandles.append(Handle('vector<float>')) #dummies
	elLabels.append(('electrons','elMass')); 					  elHandles.append(Handle('vector<float>')) #dummies
	elLabels.append(('jhuElePFlowLoose','electronLoosecharge'));  elHandles.append(Handle('vector<int>'))
	elLabels.append(('jhuElePFlowLoose','electronLooseiso')); 	  elHandles.append(Handle('vector<double>'))
	elLabels.append(('jhuElePFlowLoose','electronLooseistight')); elHandles.append(Handle('vector<unsigned int>'))
	elLabels.append(('jhuElePFlowLoose','electronLooseisloose')); elHandles.append(Handle('vector<unsigned int>')) 
	#MET
	metHandles = []; metLabels = []
	metLabels.append(('jhuGen','metpt'));  metHandles.append(Handle('double'))
	metLabels.append(('jhuGen','metphi')); metHandles.append(Handle('double'))
	#AK4 Jets, really AK5 jets
	jetHandles_AK4 = [];	jetLabels_AK4 = []
	jetLabels_AK4.append(('jhuAk5','AK5')); 			 jetHandles_AK4.append(Handle(vector_of_4vecs))
	jetLabels_AK4.append(('jetsAK4','jetAK4Eta')); 		 jetHandles_AK4.append(Handle('vector<float>')) #dummies
	jetLabels_AK4.append(('jetsAK4','jetAK4Phi')); 		 jetHandles_AK4.append(Handle('vector<float>')) #dummies
	jetLabels_AK4.append(('jetsAK4','jetAK4Mass')); 	 jetHandles_AK4.append(Handle('vector<float>')) #dummies
	jetLabels_AK4.append(('jhuAk5','AK5csv')); 			 jetHandles_AK4.append(Handle('vector<double>'))
	jetLabels_AK4.append(('jhuAk5','AK5PartonFlavour')); jetHandles_AK4.append(Handle('vector<int>'))
	#AK8 Jets, really CA8 jets
	jetHandles_AK8 = [];	jetLabels_AK8 = []
	jetLabels_AK8.append(('jhuCa8','UnprunedCA8')); 			 jetHandles_AK8.append(Handle(vector_of_4vecs))
	jetLabels_AK8.append(('patjets','patjetEta')); 				 jetHandles_AK8.append(Handle('vector<float>')) #dummies
	jetLabels_AK8.append(('patjets','patjetPhi')); 				 jetHandles_AK8.append(Handle('vector<float>')) #dummies
	jetLabels_AK8.append(('patjets','patjetMass')); 			 jetHandles_AK8.append(Handle('vector<float>')) #dummies
	jetLabels_AK8.append(('jhuCa8','UnprunedCA8csv')); 			 jetHandles_AK8.append(Handle('vector<double>'))
	jetLabels_AK8.append(('jhuCa8','UnprunedCA8tau1')); 		 jetHandles_AK8.append(Handle('vector<double>'))
	jetLabels_AK8.append(('jhuCa8','UnprunedCA8tau2')); 		 jetHandles_AK8.append(Handle('vector<double>'))
	jetLabels_AK8.append(('jhuCa8','UnprunedCA8tau3')); 		 jetHandles_AK8.append(Handle('vector<double>'))
	jetLabels_AK8.append(('jhuCa8','UnprunedCA8PartonFlavour')); jetHandles_AK8.append(Handle('vector<int>'))
	#pythia8 nTuple GenParticles
	genPartHandles = [];	genPartLabels = []
	#genPartLabels.append(('genPart','genPartPt')); 	   genPartHandles.append(Handle('vector<float>'))
	#genPartLabels.append(('genPart','genPartEta'));    genPartHandles.append(Handle('vector<float>'))
	#genPartLabels.append(('genPart','genPartPhi'));    genPartHandles.append(Handle('vector<float>'))
	#genPartLabels.append(('genPart','genPartMass'));   genPartHandles.append(Handle('vector<float>'))
	#genPartLabels.append(('genPart','genPartID')); 	   genPartHandles.append(Handle('vector<float>'))
	#genPartLabels.append(('genPart','genPartMomID'));  genPartHandles.append(Handle('vector<float>'))
	#genPartLabels.append(('genPart','genPartStatus')); genPartHandles.append(Handle('vector<float>'))
	#pileup
	pileupLabel = ('jhuGen','npv'); 	  pileupHandle 	= Handle('unsigned int')
	MCpileupLabel = ('jhuGen','npvTrue'); MCpileupHandle = Handle('unsigned int')

	##################################  ANALYZE FUNCTION  ##################################
	def analyze(self,event) :
		self.ERR_CODE = ERR_NONE
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
		#Mother particle (and MC truth top) assignment
		if self.is_data == 0 : #MC truth values only relevant for semileptonic qqbar->ttbar
			q_vec 	 = findInitialQuark(self.MC_generator,GenParticles,genPartVars) #function in eventTypeHelper.py
			qbar_vec = ROOT.TLorentzVector(q_vec.X(),q_vec.Y(),-1.0*q_vec.Z(),q_vec.E())
			MCt_vec, MCtbar_vec = findMCTops(self.MC_generator,GenParticles) #function in eventTypeHelper.py
		else : #if we don't have the MC truth information, we have to assign which is which later when we do the boost
			q_vec 	 = ROOT.TLorentzVector(1.0,0.0,sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY)
			qbar_vec = ROOT.TLorentzVector(1.0,0.0,-1.0*sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY)
			MCt_vec    = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
			MCtbar_vec = ROOT.TLorentzVector(-1.0,0.0,0.0,1.0)
		self.q_pt[0], 		self.q_eta[0] 	   = q_vec.Pt(), 	   q_vec.Eta()
		self.q_phi[0], 		self.q_M[0] 	   = q_vec.Phi(), 	   q_vec.M()
		self.qbar_pt[0], 	self.qbar_eta[0]   = qbar_vec.Pt(),    qbar_vec.Eta()
		self.qbar_phi[0], 	self.qbar_M[0] 	   = qbar_vec.Phi(),   qbar_vec.M()
		self.MCt_pt[0], 	self.MCt_eta[0]    = MCt_vec.Pt(), 	   MCt_vec.Eta()
		self.MCt_phi[0], 	self.MCt_M[0] 	   = MCt_vec.Phi(),    MCt_vec.M()
		self.MCtbar_pt[0], 	self.MCtbar_eta[0] = MCtbar_vec.Pt(),  MCtbar_vec.Eta()
		self.MCtbar_phi[0], self.MCtbar_M[0]   = MCtbar_vec.Phi(), MCtbar_vec.M()
		#get all the info from the event
		#leptons
		muVars = [];	elVars = []
		muVars_dummy = [];	elVars_dummy = []
		for i in range(len(self.muHandles)) :
			if i == 0 or i == 4 or i == 9 :
				event.getByLabel(self.muLabels[i],self.muHandles[i])
				if not self.muHandles[i].isValid() :
					self.ERR_CODE = ERR_INVALID_HANDLE
					return self.ERR_CODE
				muVars_dummy.append(self.muHandles[i].product())
			muVars.append([])
			for j in range(len(muVars_dummy[0])) :
				muVec = muVars_dummy[0][j]
				if i == 0 :	app = muVec.Pt();
				elif i == 1 :	app = muVec.Eta();
				elif i == 2 :	app = muVec.Phi();
				elif i == 3 :	app = muVec.M();
				elif i == 4 :	app = muVars_dummy[1][j];
				elif i == 9 :	app = muVars_dummy[2][j];
				else :	app = 1.0;
				muVars[i].append(app)
		for i in range(len(self.elHandles)) :
			if i == 0 or i == 4 or i == 5 or i == 6 :
				event.getByLabel(self.elLabels[i],self.elHandles[i])
				if not self.elHandles[i].isValid() :
					self.ERR_CODE = ERR_INVALID_HANDLE
					return self.ERR_CODE
				elVars_dummy.append(self.elHandles[i].product())
			elVars.append([])
			for j in range(len(elVars_dummy[0])) :
				elVec = elVars_dummy[0][j]
				if i == 0 :	app = elVec.Pt();
				elif i == 1 :	app = elVec.Eta();
				elif i == 2 :	app = elVec.Phi();
				elif i == 3 :	app = elVec.M();
				elif i == 4 :	app = elVars_dummy[1][j];
				elif i == 5 :	app = elVars_dummy[2][j];
				elif i == 6 :	app = elVars_dummy[3][j];
				else :	app = 1.0;
				elVars[i].append(app)
		#MET
		metVars = []
		for i in range(len(self.metHandles)) :
			event.getByLabel(self.metLabels[i],self.metHandles[i])
			if not self.metHandles[i].isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			metVars.append(self.metHandles[i].product())
		#AK4 Jets
		jetVars_AK4 = []
		jetVars_AK4_dummy = []
		for i in range(len(self.jetHandles_AK4)) :
			if i == 0 or i == 4 or i == 5 :
				event.getByLabel(self.jetLabels_AK4[i],self.jetHandles_AK4[i])
				if not self.jetHandles_AK4[i].isValid() :
					self.ERR_CODE = ERR_INVALID_HANDLE
					return self.ERR_CODE
				jetVars_AK4_dummy.append(self.jetHandles_AK4[i].product())
			jetVars_AK4.append([])
			for j in range(len(jetVars_AK4_dummy[0])) :
				jetVec = jetVars_AK4_dummy[0][j]
				if i == 0 :		app = jetVec.Pt();
				elif i == 1 :	app = jetVec.Eta();
				elif i == 2 :	app = jetVec.Phi();
				elif i == 3 :	app = jetVec.M();
				elif i == 4 :	app = jetVars_AK4_dummy[1][j];
				elif i == 5 :	app = jetVars_AK4_dummy[2][j];
				else :	app = 1.0;
				jetVars_AK4[i].append(app)
		#AK8 Jets
		jetVars_AK8 = []
		jetVars_AK8_dummy = []
		for i in range(len(self.jetHandles_AK8)) :
			if i == 0 or i == 4 or i == 5 or i == 6 or i == 7 or i == 8 :
				event.getByLabel(self.jetLabels_AK8[i],self.jetHandles_AK8[i])
				if not self.jetHandles_AK8[i].isValid() :
					self.ERR_CODE = ERR_INVALID_HANDLE
					return self.ERR_CODE
				jetVars_AK8_dummy.append(self.jetHandles_AK8[i].product())
			jetVars_AK8.append([])
			for j in range(len(jetVars_AK8_dummy[0])) :
				jetVec = jetVars_AK8_dummy[0][j]
				if i == 0 :	app = jetVec.Pt();
				elif i == 1 :	app = jetVec.Eta();
				elif i == 2 :	app = jetVec.Phi();
				elif i == 3 :	app = jetVec.M();
				elif i == 4 :	app = jetVars_AK8_dummy[1][j];
				elif i == 5 :	app = jetVars_AK8_dummy[2][j];
				elif i == 6 :	app = jetVars_AK8_dummy[3][j];
				elif i == 7 :	app = jetVars_AK8_dummy[4][j];
				elif i == 8 :	app = jetVars_AK8_dummy[5][j];
				else :	app = 1.0;
				jetVars_AK8[i].append(app)
		#kludge-y lepton cleaning
		#jetVars_AK4, jetVars_AK8 = leptonCleaningKludge(muVars,elVars,jetVars_AK4,jetVars_AK8)
		#pileup
		event.getByLabel(self.pileupLabel,self.pileupHandle)
		if not self.pileupHandle.isValid() :
			self.ERR_CODE = ERR_INVALID_HANDLE
			return self.ERR_CODE
		pileup = self.pileupHandle.product()[0]
		if self.is_data == 0 :
			event.getByLabel(self.MCpileupLabel,self.MCpileupHandle)
			if not self.MCpileupHandle.isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			MCpileup = self.MCpileupHandle.product()[0]
		else :
			MCpileup = 0
#		print 'pileup = '+str(pileup)+' MCpileup = '+str(MCpileup)+'' #DEBUG
		self.pileup[0] = int(pileup); self.MC_pileup[0] = int(MCpileup)
		#met cleaning
		met_cut = metCut(metVars,self.MET_control_plots) #function in metHelper.py
		if met_cut!=0 :
			return self.__closeout__(-1*met_cut)
		#lepton selection
		#get the index of the ONE valid lepton OR the cutflow failpoint
		lepIndex = leptonCuts(self.lep_type,self.side_band,muVars,elVars,metVars,jetVars_AK4,self.lepton_control_plots) 
		#above function in lepHelper.py
		if lepIndex < 0 :
			return self.__closeout__(-1*lepIndex)
		#set the fourvector of the lepton
		lep_vec, self.Q_l[0] = getLeptonFourVec(self.lep_type,muVars,elVars,lepIndex) #function in lepHelper.py
		self.lep_pt[0],  self.lep_eta[0] = lep_vec.Pt(),  lep_vec.Eta()
		self.lep_phi[0], self.lep_M[0]   = lep_vec.Phi(), lep_vec.M()
		meas_lep_pt = lep_vec.Pt(); meas_lep_eta = lep_vec.Eta()
		#neutrino handling and setup for fit
		met1_vec, met2_vec = setupMET(lep_vec,metVars) #function in metHelper.py
		self.met_pt[0], self.met_eta[0] = met1_vec.Pt(),  met1_vec.Eta() 
		self.met_phi[0], self.met_M[0]  = met1_vec.Phi(), met1_vec.M()
		if met1_vec.Pz() == met2_vec.Pz() :
			self.nFits[0] = 2
		else :
			self.nFits[0] = 1
		#jet selection
		( jetTuples, self.sf_btag_eff[0], self.sf_btag_eff_low[0], 
			self.sf_btag_eff_hi[0] ) = selectJets(self.is_data,self.top_type,lep_vec,met1_vec,met2_vec,jetVars_AK4,jetVars_AK8,self.jet_control_plots) 
		#above function in jetHelper.py
		if len(jetTuples)==1 :
			return self.__closeout__(-1*jetTuples[0])
		#event reconstruction
		lep_vec, met_vec, jetTuples, self.chi2[0] = reconstruct(lep_vec,met1_vec,met2_vec,jetTuples) 
		#above function in ttbarReconstructor.py
		#fill the TTree with the fourvector variables, and angle and differential cross section variable reconstruction
		if self.top_type == 1 :
			self.__fillFourVecsType1__(lep_vec,met_vec,jetTuples[0][0],jetTuples[1][0]) 		
			( self.cstar[0], self.x_F[0], 
			self.M[0] ) = getObservables(lep_vec+met_vec+jetTuples[0][0],jetTuples[1][0],self.Q_l[0]) 
			#above function in angleReconstructor.py
		elif self.top_type == 2 :
			self.__fillFourVecsType2__(lep_vec,met_vec,jetTuples[0][0],jetTuples[1][0],jetTuples[2][0])
			( self.cstar[0], self.x_F[0], 
			self.M[0] ) = getObservables(lep_vec+met_vec+jetTuples[0][0],jetTuples[1][0]+jetTuples[2][0],self.Q_l[0]) 
			#above function in angleReconstructor.py
		#MC Truth observable and reweighting calculation
		if self.is_data==0 :
			if self.event_type!=none :
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
	def __init__(self,fileName,isData,generator,eventType,sideband,lepType,topType,reweight,onGrid) :
		self.ERR_CODE = ERR_NONE
		#handle input options
		self.__handleInput__(fileName,isData,generator,eventType,sideband,lepType,topType,reweight)
		#book TTree
		self.__book__()
		#Set Monte Carlo reweighter
		self.corrector = MC_corrector(self.MC_generator,self.event_type,self.lep_type,onGrid)
	
	##################################   #__handleInput__   ##################################   
	#############################   __init__ helper function   ###############################
	def __handleInput__(self,fileName,isData,generator,eventType,sideband,lepType,topType,reweight) :
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
		if generator.lower() == 'powheg' or generator.lower() == 'madgraph' or generator.lower() == 'pythia8' :
			self.MC_generator = generator.lower()
		elif generator.lower() == 'mg5' :
			self.MC_generator = 'madgraph'
		elif generator.lower() == 'none' :
			self.MC_generator = 'none'
			self.jetLabels_AK4.pop(); self.jetHandles_AK4.pop(); self.jetLabels_AK8.pop(); self.jetHandles_AK8.pop()
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
		#sideband?
		if sideband == 'no' :
			print 'Analysis will be performed in SIGNAL REGION'
			self.side_band = 0
		elif sideband == 'yes' :
			print 'Analysis will be performed in ISOLATION SIDEBAND REGION'
			self.side_band = 1
		else :
			print 'ERROR: Invalid sidebanding option! Cannot run analysis!'
			print '	options.sideband = '+sideband+''
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
		#top type?
		self.top_type = topType
		if topType == 1 :
			print 'File will be analyzed using TYPE 1 (FULLY MERGED) TOPS'
		elif topType == 2 :
			print 'File will be analyzed using TYPE 2 (PARTIALLY MERGED) TOPS'
		else :
			print 'ERROR: cannot determine if top type 1 or top type 2 analysis is being performed!'
			print '	topType = '+topType+''
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
		self.addTwice = array('I',[0]); self.addBranch('addTwice',self.addTwice,'i',0)
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
		#lepton   vector
		self.lep_pt  = array('d',[-1.0]);  self.addBranch('lep_pt', self.lep_pt, 'D',-1.0)
		self.lep_eta = array('d',[100.0]); self.addBranch('lep_eta',self.lep_eta,'D',100.0)
		self.lep_phi = array('d',[100.0]); self.addBranch('lep_phi',self.lep_phi,'D',100.0)
		self.lep_M   = array('d',[-1.0]);  self.addBranch('lep_M',  self.lep_M,  'D',-1.0)
		#neutrino vector
		self.met_pt  = array('d',[-1.0]);  self.addBranch('met_pt', self.met_pt, 'D',-1.0)
		self.met_eta = array('d',[100.0]); self.addBranch('met_eta',self.met_eta,'D',100.0)
		self.met_phi = array('d',[100.0]); self.addBranch('met_phi',self.met_phi,'D',100.0)
		self.met_M   = array('d',[-1.0]);  self.addBranch('met_M',  self.met_M,  'D',-1.0)
		#leptonic W vector
		self.lepW_pt  = array('d',[-1.0]);  self.addBranch('lepW_pt', self.lepW_pt, 'D',-1.0)
		self.lepW_eta = array('d',[100.0]); self.addBranch('lepW_eta',self.lepW_eta,'D',100.0)
		self.lepW_phi = array('d',[100.0]); self.addBranch('lepW_phi',self.lepW_phi,'D',100.0)
		self.lepW_M   = array('d',[-1.0]);  self.addBranch('lepW_M',  self.lepW_M,  'D',-1.0)
		#leptonic b vector
		self.lepb_pt  = array('d',[-1.0]);  self.addBranch('lepb_pt', self.lepb_pt, 'D',-1.0)
		self.lepb_eta = array('d',[100.0]); self.addBranch('lepb_eta',self.lepb_eta,'D',100.0)
		self.lepb_phi = array('d',[100.0]); self.addBranch('lepb_phi',self.lepb_phi,'D',100.0)
		self.lepb_M   = array('d',[-1.0]);  self.addBranch('lepb_M',  self.lepb_M,  'D',-1.0)
		#hadronic W vector
		self.hadW_pt  = array('d',[-1.0]);  self.addBranch('hadW_pt', self.hadW_pt, 'D',-1.0)
		self.hadW_eta = array('d',[100.0]); self.addBranch('hadW_eta',self.hadW_eta,'D',100.0)
		self.hadW_phi = array('d',[100.0]); self.addBranch('hadW_phi',self.hadW_phi,'D',100.0)
		self.hadW_M   = array('d',[-1.0]);  self.addBranch('hadW_M',  self.hadW_M,  'D',-1.0)
		#hadronic b vector
		self.hadb_pt  = array('d',[-1.0]);  self.addBranch('hadb_pt', self.hadb_pt, 'D',-1.0)
		self.hadb_eta = array('d',[100.0]); self.addBranch('hadb_eta',self.hadb_eta,'D',100.0)
		self.hadb_phi = array('d',[100.0]); self.addBranch('hadb_phi',self.hadb_phi,'D',100.0)
		self.hadb_M   = array('d',[-1.0]);  self.addBranch('hadb_M',  self.hadb_M,  'D',-1.0)
		#leptonic top vector
		self.lept_pt  = array('d',[-1.0]);  self.addBranch('lept_pt', self.lept_pt, 'D',-1.0)
		self.lept_eta = array('d',[100.0]); self.addBranch('lept_eta',self.lept_eta,'D',100.0)
		self.lept_phi = array('d',[100.0]); self.addBranch('lept_phi',self.lept_phi,'D',100.0)
		self.lept_M   = array('d',[-1.0]);  self.addBranch('lept_M',  self.lept_M,  'D',-1.0)
		#hadronic top vector
		self.hadt_pt  = array('d',[-1.0]);  self.addBranch('hadt_pt', self.hadt_pt, 'D',-1.0)
		self.hadt_eta = array('d',[100.0]); self.addBranch('hadt_eta',self.hadt_eta,'D',100.0)
		self.hadt_phi = array('d',[100.0]); self.addBranch('hadt_phi',self.hadt_phi,'D',100.0)
		self.hadt_M   = array('d',[-1.0]);  self.addBranch('hadt_M',  self.hadt_M,  'D',-1.0)
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
		#cosine(theta)
		self.cstar     = array('d',[100.0]); self.addBranch('cstar',   self.cstar,   'D',100.0)
		self.cstar_MC  = array('d',[100.0]); self.addBranch('cstar_MC',self.cstar_MC,'D',100.0)
		#Feynman x
		self.x_F  	 = array('d',[100.0]); self.addBranch('x_F',   self.x_F,   'D',100.0)
		self.x_F_MC  = array('d',[100.0]); self.addBranch('x_F_MC',self.x_F_MC,'D',100.0)
		#ttbar invariant mass
		self.M     = array('d',[-1.0]); self.addBranch('M',   self.M,   'D',-1.0)
		self.M_MC  = array('d',[-1.0]); self.addBranch('M_MC',self.M_MC,'D',-1.0)
		#pileup
		self.pileup = array('i',[0]); self.addBranch('pileup',self.pileup,'I',0)
		self.MC_pileup = array('i',[0]); self.addBranch('MC_pileup',self.MC_pileup,'I',0)

		self.control_plots_folder = self.f.mkdir('control_plots_folder','control_plots')
		self.control_plots_folder.cd()
		#control plots for optimising cuts
		self.MET_control_plots = []; self.lepton_control_plots = []; self.jet_control_plots = []
		#MET
		self.MET_control_plots.append(ROOT.TH1F('met_pt','p_{T} of MET; GeV',60,0.,300.))
		#leptons
		self.lepton_control_plots.append(ROOT.TH1F('lep1_pt','p_{T} of first lepton; p_{T} (GeV)',60,0.,300.))
		self.lepton_control_plots.append(ROOT.TH1F('lep1_eta','#eta of first lepton; #eta',60,-3.0,3.0))
		self.lepton_control_plots.append(ROOT.TH1F('lep2_pt','p_{T} of second lepton; p_{T} (GeV)',60,0.,300.))
		self.lepton_control_plots.append(ROOT.TH1F('lep2_eta','#eta of second lepton; #eta',60,-3.0,3.0))
		self.lepton_control_plots.append(ROOT.TH1F('other_lep1_pt',
			'p_{T} of first "other" lepton; p_{T} (GeV)',60,0.,300.))
		self.lepton_control_plots.append(ROOT.TH1F('other_lep1_eta','#eta of first "other" lepton; #eta',60,-3.0,3.0))
		self.lepton_control_plots.append(ROOT.TH1F('other_lep2_pt',
			'p_{T} of second "other" lepton; p_{T} (GeV)',60,0.,300.))
		self.lepton_control_plots.append(ROOT.TH1F('other_lep2_eta','#eta of second "other" lepton; #eta',60,-3.0,3.0))
		self.lepton_control_plots.append(ROOT.TH2F('lep_2D_cut',
			'selected lepton #Delta R and p_{T,rel} from nearest jet; #Delta R; p_{T,rel} (GeV)',75,0.0,1.5,75,0.0,75.))
		self.lepton_control_plots.append(ROOT.TH2F('triangle_cut_el',
			'|#Delta #phi(e^{-},E^{miss}_{T})-1.5| vs. E^{miss}_{T}; E^{miss}_{T} (GeV); |#Delta #phi(e^{-},E^{miss}_{T})-1.5|',
			60,0,300,50,0.,5.))
		self.lepton_control_plots.append(ROOT.TH2F('triangle_cut_jet',
			'|#Delta #phi(jet,E^{miss}_{T})-1.5| vs. E^{miss}_{T}; E^{miss}_{T} (GeV); |#Delta #phi(jet,E^{miss}_{T})-1.5|',
			60,0,300,50,0.,5.))
		#jets
		self.jet_control_plots.append(ROOT.TH1F('lep_bjet_pT','p_{T} of small jets; p_{T} (GeV)',50,0.0,350.0))
		self.jet_control_plots.append(ROOT.TH1F('lep_bjet_dR',
			'#Delta R(lepton) of loosely b-tagged small jets; #Delta R',50,0.,5.0))
		self.jet_control_plots.append(ROOT.TH1F('lep_bjet_comb_mass',
			'Best combined mass of (lepton, MET, leptonic bjet); M (GeV)',100,0.,500.))
		self.jet_control_plots.append(ROOT.TH1F('lep_bjet_CSV','CSV of small jets; CSV value',20,0.0,1.0))
		self.jet_control_plots.append(ROOT.TH1F('t1_top_pT','p_{T} of big jets; p_{T} (GeV)',100,0.0,500.0))
		self.jet_control_plots.append(ROOT.TH1F('t1_top_mass','Mass of big jets; M (GeV)',60,0.0,300.0))
		self.jet_control_plots.append(ROOT.TH1F('t1_top_tau32','#tau_{3}/#tau_{2} of big jets; #tau_{32}',20,0.0,1.0))
		self.jet_control_plots.append(ROOT.TH1F('t1_top_dR','#Delta R(lepton) of top candidates; #Delta R',50,0.0,5.0))
		self.jet_control_plots.append(ROOT.TH1F('t1_top_mult','number of hadronic top candidates',10,0.0,10.0))
		self.jet_control_plots.append(ROOT.TH1F('t2_top_W_pT','p_{T} of big jets; p_{T} (GeV)',100,0.0,500.0))
		self.jet_control_plots.append(ROOT.TH1F('t2_top_W_tau21','#tau_{2}/#tau_{1} of big jets; #tau_{21}',20,0.0,1.0))
		self.jet_control_plots.append(ROOT.TH1F('t2_top_W_mass','Mass of hadronic W jets; M (GeV)',60,0.0,300.0))
		self.jet_control_plots.append(ROOT.TH1F('t2_top_W_dR',
			'#Delta R(lepton) of hadronic W jets; #Delta R',50,0.0,5.0))
		self.jet_control_plots.append(ROOT.TH1F('t2_top_W_mult','number of hadronic W candidates',10,0.0,10.0))
		self.jet_control_plots.append(ROOT.TH1F('t2_top_comb_mass',
			'type2 top candidate combined mass; M (GeV)',60,0.0,300.0))
		self.jet_control_plots.append(ROOT.TH1F('t2_top_b_dR',
			'#Delta R(lepton) of hadronic side b candidates; #Delta R',50,0.0,5.0))
		self.jet_control_plots.append(ROOT.TH1F('t2_top_b_W_dR',
			'#Delta R(hadronic W) of hadronic side b candidates; #Delta R',50,0.0,5.0))
		self.jet_control_plots.append(ROOT.TH1F('t2_top_had_b_mult','number of hadronic side b candidates',10,0.0,10.0))
		self.all_control_plots = self.MET_control_plots + self.lepton_control_plots + self.jet_control_plots 

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

	################ function to fill fourvector branch values (type 2 tops) ###############
	def __fillFourVecsType2__(self,lepton,met,lepb,hadW,hadb) :
		self.lep_pt[0] 	= lepton.Pt(); 	self.lep_eta[0] = lepton.Eta()
		self.lep_phi[0] = lepton.Phi(); self.lep_M[0] 	= lepton.M()
		self.met_pt[0] 	= met.Pt(); 	self.met_eta[0] = met.Eta()
		self.met_phi[0] = met.Phi(); 	self.met_M[0] 	= met.M()
		lepW = lepton+met
		self.lepW_pt[0]  = lepW.Pt();  self.lepW_eta[0] = lepW.Eta()
		self.lepW_phi[0] = lepW.Phi(); self.lepW_M[0] 	= lepW.M()
		self.lepb_pt[0]  = lepb.Pt();  self.lepb_eta[0] = lepb.Eta()
		self.lepb_phi[0] = lepb.Phi(); self.lepb_M[0] 	= lepb.M()
		self.hadW_pt[0]  = hadW.Pt();  self.hadW_eta[0] = hadW.Eta()
		self.hadW_phi[0] = hadW.Phi(); self.hadW_M[0] 	= hadW.M()
		self.hadb_pt[0]  = hadb.Pt();  self.hadb_eta[0] = hadb.Eta()
		self.hadb_phi[0] = hadb.Phi(); self.hadb_M[0] 	= hadb.M()
		lept = lepW+lepb
		self.lept_pt[0]  = lept.Pt();  self.lept_eta[0] = lept.Eta()
		self.lept_phi[0] = lept.Phi(); self.lept_M[0] 	= lept.M()
		hadt = hadW+hadb
		self.hadt_pt[0]  = hadt.Pt();  self.hadt_eta[0] = hadt.Eta()
		self.hadt_phi[0] = hadt.Phi(); self.hadt_M[0] 	= hadt.M()

	################ function to fill fourvector branch values (type 1 tops) ###############
	def __fillFourVecsType1__(self,lepton,met,lepb,hadt) :
		self.lep_pt[0] 	= lepton.Pt(); 	self.lep_eta[0] = lepton.Eta()
		self.lep_phi[0] = lepton.Phi(); self.lep_M[0] 	= lepton.M()
		self.met_pt[0] 	= met.Pt(); 	self.met_eta[0] = met.Eta()
		self.met_phi[0] = met.Phi(); 	self.met_M[0] 	= met.M()
		lepW = lepton+met
		self.lepW_pt[0]  = lepW.Pt();  self.lepW_eta[0] = lepW.Eta()
		self.lepW_phi[0] = lepW.Phi(); self.lepW_M[0] 	= lepW.M()
		self.lepb_pt[0]  = lepb.Pt();  self.lepb_eta[0] = lepb.Eta()
		self.lepb_phi[0] = lepb.Phi(); self.lepb_M[0] 	= lepb.M()
		lept = lepW+lepb
		self.lept_pt[0]  = lept.Pt();  self.lept_eta[0] = lept.Eta()
		self.lept_phi[0] = lept.Phi(); self.lept_M[0] 	= lept.M()
		self.hadt_pt[0]  = hadt.Pt();  self.hadt_eta[0] = hadt.Eta()
		self.hadt_phi[0] = hadt.Phi(); self.hadt_M[0] 	= hadt.M()

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