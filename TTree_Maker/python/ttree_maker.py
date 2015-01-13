
#TTree Maker workhorse code
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
from lepHelper import leptonCuts, getLeptonFourVec
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
	#muons
	muHandles = []; muLabels = []
	muPtHandle  	= Handle('vector<float>'); muHandles.append(muPtHandle); 	  muPtLabel  	 = ('muons','muPt'); 				 muLabels.append(muPtLabel)
	muEtaHandle 	= Handle('vector<float>'); muHandles.append(muEtaHandle); 	  muEtaLabel 	 = ('muons','muEta'); 				 muLabels.append(muEtaLabel)
	muPhiHandle 	= Handle('vector<float>'); muHandles.append(muPhiHandle); 	  muPhiLabel 	 = ('muons','muPhi'); 				 muLabels.append(muPhiLabel)
	muMHandle   	= Handle('vector<float>'); muHandles.append(muMHandle); 	  muMLabel   	 = ('muons','muMass'); 				 muLabels.append(muMLabel)
	muChargeHandle  = Handle('vector<float>'); muHandles.append(muChargeHandle);  muChargeLabel  = ('muons','muCharge'); 			 muLabels.append(muChargeLabel)
	muSumCHPtHandle = Handle('vector<float>'); muHandles.append(muSumCHPtHandle); muSumCHPtLabel = ('muons','muSumChargedHadronPt'); muLabels.append(muSumCHPtLabel)
	muSumNHPtHandle = Handle('vector<float>'); muHandles.append(muSumNHPtHandle); muSumNHPtLabel = ('muons','muSumNeutralHadronPt'); muLabels.append(muSumNHPtLabel)
	muSumPhPtHandle = Handle('vector<float>'); muHandles.append(muSumPhPtHandle); muSumPhPtLabel = ('muons','muSumPhotonPt'); 		 muLabels.append(muSumPhPtLabel)
	muSumPUPtHandle = Handle('vector<float>'); muHandles.append(muSumPUPtHandle); muSumPUPtLabel = ('muons','muSumPUPt'); 			 muLabels.append(muSumPUPtLabel)
	muIsTightHandle = Handle('vector<float>'); muHandles.append(muIsTightHandle); muIsTightLabel = ('muons','muIsTightMuon'); 		 muLabels.append(muIsTightLabel)
	#electrons
	elHandles = []; elLabels = []
	elPtHandle  	= Handle('vector<float>'); elHandles.append(elPtHandle);     elPtLabel     = ('electrons','elPt');     elLabels.append(elPtLabel)
	elEtaHandle 	= Handle('vector<float>'); elHandles.append(elEtaHandle);    elEtaLabel    = ('electrons','elEta');    elLabels.append(elEtaLabel)
	elPhiHandle 	= Handle('vector<float>'); elHandles.append(elPhiHandle);    elPhiLabel    = ('electrons','elPhi');    elLabels.append(elPhiLabel)
	elMHandle   	= Handle('vector<float>'); elHandles.append(elMHandle);      elMLabel      = ('electrons','elMass');   elLabels.append(elMLabel)
	elChargeHandle  = Handle('vector<float>'); elHandles.append(elChargeHandle); elChargeLabel = ('electrons','elCharge'); elLabels.append(elChargeLabel)
	elIso03Handle   = Handle('vector<float>'); elHandles.append(elIso03Handle);  elIso03Label  = ('electrons','elIso03');  elLabels.append(elIso03Label)
	#MET
	metHandles = []; metLabels = []
	metPtHandle  = Handle('vector<float>'); metHandles.append(metPtHandle);  metPtLabel  = ('met','metPt');  metLabels.append(metPtLabel)
	metPhiHandle = Handle('vector<float>'); metHandles.append(metPhiHandle); metPhiLabel = ('met','metPhi'); metLabels.append(metPhiLabel)
	#AK4 Jets
	jetHandles_AK4 = [];	jetLabels_AK4 = []
	AK4jetPtHandle   = Handle('vector<float>'); jetHandles_AK4.append(jetPtHandle);  AK4jetPtLabel  = ('jets','jetPt');   jetLabels_AK4.append(jetPtLabel)
	AK4jetEtaHandle  = Handle('vector<float>'); jetHandles_AK4.append(jetEtaHandle); AK4jetEtaLabel = ('jets','jetEta');  jetLabels_AK4.append(jetEtaLabel)
	AK4jetPhiHandle  = Handle('vector<float>'); jetHandles_AK4.append(jetPhiHandle); AK4jetPhiLabel = ('jets','jetPhi');  jetLabels_AK4.append(jetPhiLabel)
	AK4jetMHandle    = Handle('vector<float>'); jetHandles_AK4.append(jetMHandle);   AK4jetMLabel   = ('jets','jetMass'); jetLabels_AK4.append(jetMLabel)
	AK4jetCSVHandle  = Handle('vector<float>'); jetHandles_AK4.append(jetCSVHandle); AK4jetCSVLabel = ('jets','jetCSV');  jetLabels_AK4.append(jetCSVLabel)
	#AK8 Jets
	jetHandles_AK8 = [];	jetLabels_AK8 = []
	AK8jetPtHandle	  = Handle('vector<float>'); jetHandles_AK8.append(jetPtHandle);   AK8jetPtLabel   = ('jetsAK8','jetAK8Pt');    jetLabels_AK8.append(jetPtLabel)
	AK8jetEtaHandle	  = Handle('vector<float>'); jetHandles_AK8.append(jetEtaHandle);  AK8jetEtaLabel  = ('jetsAK8','jetAK8Eta');   jetLabels_AK8.append(jetEtaLabel)
	AK8jetPhiHandle	  = Handle('vector<float>'); jetHandles_AK8.append(jetPhiHandle);  AK8jetPhiLabel  = ('jetsAK8','jetAK8Phi');   jetLabels_AK8.append(jetPhiLabel)
	AK8jetMHandle	  = Handle('vector<float>'); jetHandles_AK8.append(jetMHandle);	   AK8jetMLabel    = ('jetsAK8','jetAK8Mass');  jetLabels_AK8.append(jetMLabel)
	AK8jetCSVHandle	  = Handle('vector<float>'); jetHandles_AK8.append(jetCSVHandle);  AK8jetCSVLabel  = ('jetsAK8','jetAK8CSV');   jetLabels_AK8.append(jetCSVLabel)
	AK8jettau1Handle  = Handle('vector<float>'); jetHandles_AK8.append(jettau1Handle); AK8jettau1Label = ('jetsAK8','jetAK8tau1');  jetLabels_AK8.append(jettau1Label)
	AK8jettau2Handle  = Handle('vector<float>'); jetHandles_AK8.append(jettau2Handle); AK8jettau2Label = ('jetsAK8','jetAK8tau2');  jetLabels_AK8.append(jettau2Label)
	AK8jettau3Handle  = Handle('vector<float>'); jetHandles_AK8.append(jettau3Handle); AK8jettau3Label = ('jetsAK8','jetAK8tau3');  jetLabels_AK8.append(jettau3Label)

	##################################  ANALYZE FUNCTION  ##################################
	def analyze(self,event) :
		#keep track of whether event has been cut
		keepEvent = True
		#event type split
		if self.is_data == 0 :
			event.getByLabel(self.genLabel,self.genHandle)
			if not self.genHandle.isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			GenParticles = self.genHandle.product()
			if self.event_type != 4 :
				keepEvent,add_twice = eventTypeCheck(self.MC_generator,GenParticles,self.event_type) #function in eventTypeHelper.py
				if add_twice :
					self.addTwice[0] = 1
		if not keepEvent :
			return ERR_NONE
		#Mother particle (and MC truth top) assignment
		if self.is_data == 0 : #MC truth values only relevant for semileptonic qqbar->ttbar
			q_vec 	 = findInitialQuark(self.MC_generator,GenParticles) #function in eventTypeHelper.py
			qbar_vec = ROOT.TLorentzVector(q_vec.X(),q_vec.Y(),-1.0*q_vec.Z(),q_vec.E())
			MCt_vec, MCtbar_vec = findMCTops(GenParticles) #function in eventTypeHelper.py
		else : #if we don't have the MC truth information, we have to assign which is which later when we do the boost
			q_vec 	 = ROOT.TLorentzVector(0.000000001,0.0,sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY)
			qbar_vec = ROOT.TLorentzVector(0.000000001,0.0,-1.0*sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY)
			MCt_vec    = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
			MCtbar_vec = ROOT.TLorentzVector(-1.0,0.0,0.0,1.0)
		self.q_pt[0], 	   self.q_eta[0], 	   self.q_phi[0], 	   self.q_M[0] 		= q_vec.Pt(), 	   q_vec.Eta(), 	 q_vec.Phi(), 	   q_vec.M()
		self.qbar_pt[0],   self.qbar_eta[0],   self.qbar_phi[0],   self.qbar_M[0]   = qbar_vec.Pt(),   qbar_vec.Eta(),   qbar_vec.Phi(),   qbar_vec.M()
		self.MCt_pt[0],    self.MCt_eta[0],    self.MCt_phi[0],    self.MCt_M[0]    = MCt_vec.Pt(),    MCt_vec.Eta(),    MCt_vec.Phi(),    MCt_vec.M()
		self.MCtbar_pt[0], self.MCtbar_eta[0], self.MCtbar_phi[0], self.MCtbar_M[0] = MCtbar_vec.Pt(), MCtbar_vec.Eta(), MCtbar_vec.Phi(), MCtbar_vec.M()
		#get all the info from the event
		#leptons
		muVars = [];	elVars = []
		for i in range(len(self.muHandles)) :
			event.getByLabel(self.muLabels[i],self.muHandles[i])
			if not self.muHandles[i].isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			muVars.append(self.muHandles[i].product())
		for i in range(len(self.elHandles)) :
			event.getByLabel(self.elLabels[i],self.elHandles[i])
			if not self.elHandles[i].isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			elVars.append(self.elHandles[i].product())
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
		for i in range(len(self.jetHandles_AK4)) :
			event.getByLabel(self.jetLabels_AK4[i],self.jetHandles_AK4[i])
			if not self.jetHandles_AK4[i].isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			jetVars_AK4.append(self.jetHandles_AK4[i].product())
		#AK8 Jets
		jetVars_AK8 = []
		for i in range(len(self.jetHandles_AK8)) :
			event.getByLabel(self.jetLabels_AK8[i],self.jetHandles_AK8[i])
			if not self.jetHandles_AK8[i].isValid() :
				self.ERR_CODE = ERR_INVALID_HANDLE
				return self.ERR_CODE
			jetVars_AK8.append(self.jetHandles_AK8[i].product())
		#met cleaning
		met_cut = metCut(metVars) #function in metHelper.py
		if met_cut!=0 :
			return self.__closeout__(-1*met_cut)
		#lepton selection
		#get the index of the ONE valid lepton OR the cutflow failpoint
		lepIndex = leptonCuts(self.lep_type,self.side_band,muVars,elVars,jetVars) #function in lepHelper.py
		#set the fourvector of the lepton
		lep_vec, self.Q_l[0] = getLeptonFourVec(self.lep_type,muVars,elVars,lepIndex) #function in lepHelper.py
		self.lep_pt[0], self.lep_eta[0], self.lep_phi[0], self.lep_M[0] = lep_vec.Pt(), lep_vec.Eta(), lep_vec.Phi(), lep_vec.M()
		#neutrino handling and setup for fit
		met1_vec, met2_vec = setupMET(lep_vec,metVars) #function in metHelper.py
		self.met_pt[0], self.met_eta[0], self.met_phi[0], self.met_M[0] = met1_vec.Pt(), met1_vec.Eta(), met1_vec.Phi(), met1_vec.M()
		if met1_vec.Pz() == met2_vec.Pz() :
			self.nFits[0] = 2
		else :
			self.nFits[0] = 1
		#jet selection
		jetTuples = selectJets(self.top_type,lep_vec,met1_vec,met2_vec,jetVars_AK4,jetVars_AK8) #function in jetHelper.py
		#event reconstruction
		lep_vec, met_vec, jetTuples, self.chi2[0] = reconstruct(lep_vec,met1_vec,met2_vec,jetTuples) #function in ttbarReconstructor.py
		#populate the TTree with the fourvector variables, and angle and differential cross section variable reconstruction
		if self.top_type == 1 :
			self.__fillFourVecs__(lep_vec,met_vec,jetTuples[0][0],jetTuples[1][0]) 
		
			self.cstar[0], self.x_F[0], self.M[0] = getObservables(lep_vec+met_vec+jetTuples[0][0],jetTuples[1][0],self.Q_l[0]) #function in angleReconstructor.py
		elif self.top_type == 2 :
			self.__fillFourVecs__(lep_vec,met_vec,jetTuples[0][0],jetTuples[1][0],jetTuples[2][0])
			self.cstar[0], self.x_F[0], self.M[0] = getObservables(lep_vec+met_vec+jetTuples[0][0],jetTuples[1][0]+jetTuples[2][0],self.Q_l[0]) #function in angleReconstructor.py
		#MC Truth observable and reweighting calculation
		if self.is_data==0 :
			( self.cstar_MC[0],self.x_F_MC[0],self.M_MC[0],
				self.w_a[0],self.w_s_xi[0],self.w_a_xi[0],self.w_s_delta[0],self.w_a_delta[0],
				self.w_a_opp[0],self.w_s_xi_opp[0],self.w_a_xi_opp[0],self.w_s_delta_opp[0],self.w_a_delta_opp[0] ) = getMCObservables(q_vec,qbar_vec,MCt_vec,MCtbar_vec) 
		#scale factor and reweighting calculations
			#These are going to have to get added as we go with them, functions in eventWeightCalculator.py
		self.__closeout__(0) #yay! A successful event!

	##################################  #__init__ function  ##################################
	def __init__(self,fileName,isData,generator,eventType,sideband,lepType,topType,reweight) :
		self.ERR_CODE = ERR_NONE
		#handle input options
		self.__handleInput__(fileName,isData,generator,eventType,sideband,lepType,reweight)
		#book TTree
		self.__book__()
	
	##################################   #__handleInput__   ##################################   
	#############################   __init__ helper function   ###############################
	def __handleInput__(self,fileName,isData,generator,eventType,sideband,lepType,topType,reweight) :
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
		self.MC_generator = generator
		#event type?
		if eventType == 'qq_semilep' or eventType == 'qqbar_semilep' or eventType == 'qq' or eventType == 'qqbar' :
			print 'only SEMILEPTONIC qqbar EVENTS will be analyzed from this file'
			self.event_type = 0
		elif eventType == 'gg_semilep' or eventType == 'gg' :
			print 'only SEMILEPTONIC gg (qg,qiqbarj,etc.) EVENTS will be analyzed from this file'
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
		if topType == 1 :
			print 'File will be analyzed using TYPE 1 (FULLY MERGED) TOPS'
			self.top_type = 0
		elif topType == 2 :
			print 'File will be analyzed using TYPE 2 (PARTIALLY MERGED) TOPS'
			self.top_type = 1
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
	def __fillFourVecs__(self,lepton,met,lepb,hadW,hadb) :
		self.lep_pt[0]   = lepton.Pt();	self.lep_eta[0]   = lepton.Eta(); self.lep_phi[0]   = lepton.Phi(); self.lep_M[0]   = lepton.M()
		self.met_pt[0]   = met.Pt();	self.met_eta[0]   = met.Eta();	  self.met_phi[0]   = met.Phi();	self.met_M[0]   = met.M()
		lepW = lepton+met
		self.lepW_pt[0]  = lepW.Pt();	self.lepW_eta[0]  = lepW.Eta();	  self.lepW_phi[0]  = lepW.Phi();   self.lepW_M[0]  = lepW.M()
		self.lepb_pt[0]  = lepb.Pt();	self.lepb_eta[0]  = lepb.Eta();	  self.lepb_phi[0]  = lepb.Phi();   self.lepb_M[0]  = lepb.M()
		self.hadW_pt[0]  = hadW.Pt();	self.hadW_eta[0]  = hadW.Eta();   self.hadW_phi[0]  = hadW.Phi();	self.hadW_M[0]  = hadW.M()
		self.hadb_pt[0]  = hadb.Pt();	self.hadb_eta[0]  = hadb.Eta();   self.hadb_phi[0]  = hadb.Phi();	self.hadb_M[0]  = hadb.M()
		lept = lepW+lepb
		self.lept_pt[0]  = lept.Pt();	self.lept_eta[0]  = lept.Eta();   self.lept_phi[0]  = lept.Phi();	self.lept_M[0]  = lept.M()
		hadt = hadW+hadb
		self.hadt_pt[0]  = hadt.Pt();	self.hadt_eta[0]  = hadt.Eta();   self.hadt_phi[0]  = hadt.Phi();	self.hadt_M[0]  = hadt.M()

	################ function to fill fourvector branch values (type 2 tops) ###############
	def __fillFourVecs__(self,lepton,met,lepb,hadt) :
		self.lep_pt[0]   = lepton.Pt();	self.lep_eta[0]   = lepton.Eta(); self.lep_phi[0]   = lepton.Phi(); self.lep_M[0]   = lepton.M()
		self.met_pt[0]   = met.Pt();	self.met_eta[0]   = met.Eta();	  self.met_phi[0]   = met.Phi();	self.met_M[0]   = met.M()
		lepW = lepton+met
		self.lepW_pt[0]  = lepW.Pt();	self.lepW_eta[0]  = lepW.Eta();	  self.lepW_phi[0]  = lepW.Phi();   self.lepW_M[0]  = lepW.M()
		self.lepb_pt[0]  = lepb.Pt();	self.lepb_eta[0]  = lepb.Eta();	  self.lepb_phi[0]  = lepb.Phi();   self.lepb_M[0]  = lepb.M()
		lept = lepW+lepb
		self.lept_pt[0]  = lept.Pt();	self.lept_eta[0]  = lept.Eta();   self.lept_phi[0]  = lept.Phi();	self.lept_M[0]  = lept.M()
		self.hadt_pt[0]  = hadt.Pt();	self.hadt_eta[0]  = hadt.Eta();   self.hadt_phi[0]  = hadt.Phi();	self.hadt_M[0]  = hadt.M()

	########## function to close out the event, called before kicking back to runner #########
	def __closeout__(self,cut_flow) :
		#update cutflow
		self.cutflow[0] = cut_flow
		#update error code
		self.error_code[0] = self.ERR_CODE
		#fill ttree
		self.tree.Fill()
		#return error code
		return self.ERR_CODE

	################################## __del__ function  ###################################
	def __del__(self) :
		self.f.cd()
		self.f.Write()
		self.f.Close()