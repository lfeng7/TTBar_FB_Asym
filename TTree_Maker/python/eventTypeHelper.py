import ROOT
from math import *

#eventTypeHelper: holds helper functions for doing the event type split 
#and for setting the fourvector of the initial state quark if the file 
#is a Monte Carlo file

#global variables
#Beam energy
SQRT_S=13000.0
BEAM_ENERGY=SQRT_S/2.0
#PDG ID Numbering scheme http://pdg.lbl.gov/2002/montecarlorpp.pdf
PROTON_ID = 2212
TOP_ID    = 6
W_ID      = 24
ELECTRON_ID = 11
TAU_NEUTRINO_ID = 18

#################################  REFERENCED FUNCTIONS  ##################################

#eventTypeCheck
#takes in: generator GenParticle tree, event type indicator
#returns: boolean value of whether event type matches
def eventTypeCheck(generator,GenParticles,event_type) :
	#event types:
	#	0 = semileptonic qqbar
	#	1 = semileptonic gg, etc.
	#	2 = dileptonic
	#	3 = hadronic
	if generator == 'powheg' or generator == 'madgraph' or generator == 'mg5' or generator == 'mg' :
		return typeCheckPowheg(GenParticles,event_type)
	elif generator == 'mcatnlo' :
		return typeCheckMCAtNLO(GenParticles,event_type)
	elif generator == 'pythia8' :
		return typeCheckPythia8(GenParticles,event_type)
	else :
		print 'ERROR: GENERATOR '+generator+' NOT RECOGNIZED!!!'
		return False

#findInitialQuark
#takes in GenParticle tree
#returns TLorentzVector of initial state QUARK based on MC truth information
#NOTE: this just gives one direction or the other down the beam pipe
def findInitialQuark(generator,GenParticles) :
	if generator == 'powheg' or generator == 'madgraph' or generator == 'mg5' or generator == 'mg' :
		return findInitialQuarkPowheg(GenParticles)
	elif generator == 'mcatnlo' :
		return findInitialQuarkMCAtNLO(GenParticles)
	elif generator == 'pythia8' :
		return findInitialQuarkPythia8(GenParticles)
	else :
		print 'ERROR: GENERATOR '+generator+' NOT RECOGNIZED!!!'
		return ROOT.TLorentzVector(1.0,0.0,sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY) #DEBUG RETURN

#findMCTops
#takes in GenParticle tree
#returns a tuple of the (top, antitop) fourvectors from Monte Carlo
def findMCTops(GenParticles) :
	#same algorithm for every generator
	found_t    = False
	found_tbar = False
	tvec, tbarvec = (ROOT.TLorentzVector(1.0,0.0,0.0,1.0),ROOT.TLorentzVector(1.0,0.0,0.0,1.0))
	for p in GenParticles :
		if p.pt()<0 :
			continue
		if found_t and found_tbar :
			break
		if p.pdgId() == TOP_ID and not found_t :
			tvec.SetPtEtaPhiM(p.pt(),p.eta(),p.phi(),p.M())
		elif p.pdgId() == -1*TOP_ID and not found_tbar :
			tbarvec.SetPtEtaPhiM(p.pt(),p.eta(),p.phi(),p.M())
	return (tvec,tbarvec)
	#return (ROOT.TLorentzVector(1.0,0.0,0.0,1.0),ROOT.TLorentzVector(1.0,0.0,0.0,1.0)) #DEBUG RETURN

##################################  HELPER FUNCTIONS  ##################################

#Powheg eventTypeCheck function
def typeCheckPowheg(GenParticles,event_type) :
	#check the initial state partons if necessary
	if event_type<2 :
		initial_state_parton_ids = []
		#loop through and find the pdg IDs of the first daughters of the protons
		for p in GenParticles :
			if p.pt()<0 or p.status()!=3 :
				continue
			if p.pdgId()==PROTON_ID :
				initial_state_parton_ids.append(p.daughter(0).pdgId())
		#is it a qqbar event?
		is_qq = len(initial_state_parton_ids) == 2 and initial_state_parton_ids[0]+initial_state_parton_ids[1]==0
		#return false if the event type is incorrect
		if (event_type == 0 and not is_qq) or (event_type == 1 and is_qq) :
			return False
	return semilepCheck(GenParticles,event_type)
	#return True #DEBUG RETURN

#MC@NLO eventTypeCheck function
def typeCheckMCAtNLO(GenParticles,event_type) :
	#check the initial state partons if necessary
	if event_type<2 :
		initial_state_parton_ids = []
		#loop through and find the particles whose daughters include the ttbar pair
		for p in GenParticles :
			if p.pt()<0 or p.status()!=3 :
				continue
			ntopDaus = 0
			for i in range(p.numberOfDaughters()) :
				if fabs(p.daughter(i).pdgId()) == TOP_ID :
					ntopDaus+=1
			if ntopDaus == 2 :
				initial_state_parton_ids.append(p.pdgId())
		#is it a qqbar event?
		is_qq = len(initial_state_parton_ids) == 2 and initial_state_parton_ids[0]+initial_state_parton_ids[1]==0
		#return false if the event type is incorrect
		if (event_type == 0 and not is_qq) or (event_type == 1 and is_qq) :
			return False
	return semilepCheck(GenParticles,event_type)
	#return True #DEBUG RETURN

#pythia8 eventTypeCheck function
def typeCheckPythia8(GenParticles,event_type) :
	#check the initial state partons if necessary
	if event_type<2 :
		initial_state_parton_ids = []
		#INSERT ALGORITHM HERE
		#is it a qqbar event?
		is_qq = len(initial_state_parton_ids) == 2 and initial_state_parton_ids[0]+initial_state_parton_ids[1]==0
		#if (event_type == 0 and not is_qq) or (event_type == 1 and is_qq) :
		#	return False

	#check the decay type
	n_leps_from_Ws = 0
	for p in GenParticles :
		if p.pt()<0 or p.status()!=23 :
			continue
		#look for leptons
		if fabs(p.pdgId()) in range(ELECTRON_ID,TAU_NEUTRINO_ID+1) and fabs(p.mother(0).pdgId()) :
			n_leps_from_Ws+=1
	if (n_leps_from_Ws == 2 and event_type > 1) or (n_leps_from_Ws == 4 and event_type != 2) or (n_leps_from_Ws == 0 and event_type != 3) :
		return False
	
	return True #DEBUG RETURN

#semileptonic checker
def semilepCheck(GenParticles,event_type) :
	#check the decay type
	leptons_from_Ws_IDs = []
	for p in GenParticles :
		if p.pt()<0 or p.status()!=3 :
			continue
		#look for the Ws that are the decay products of the tops
		if fabs(p.pdgId())==TOP_ID :
			for i in range(p.numberOfDaughters()) :
				if fabs(p.daughter(i).pdgId()) == W_ID :
					w_p = p.daughter(i)
					#get the IDs of all its leptonic daughters and add them to the list
					for j in range(w_p.numberOfDaughters()) :
						if fabs(w_p.daughter(j).pdgId()) in range(ELECTRON_ID,TAU_NEUTRINO_ID+1) :
							leptons_from_Ws_IDs.append(w_p.daughter(j).pdgId())
	if (len(leptons_from_Ws_IDs) == 2 and event_type > 1) or (len(leptons_from_Ws_IDs) == 4 and event_type != 2) or (len(leptons_from_Ws_IDs) == 0 and event_type != 3) :
		return False
	return True
	#return True #DEBUG RETURN

#Powheg findInitialQuark function
def findInitialQuarkPowheg(GenParticles) :
	for p in GenParticles :
		if p.pt()<0 or p.status()!=3 :
			continue
		if p.pdgId() == PROTON_ID and p.daughter(0).pdgId() > 0:
			factor = 0.0
			if p.daughter(0).eta() > 0 :
				factor = 1.0
			elif p.daughter(0).eta() < 0 :
				factor = -1.0
			return ROOT.TLorentzVector(1.0,0.0,factor*sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY)
	#return ROOT.TLorentzVector(1.0,0.0,sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY) #DEBUG RETURN

#MC@NLO findInitialQuark function
def findInitialQuarkMCAtNLO(GenParticles) :
	for p in GenParticles :
		if p.pt()<0 or p.status()!=3 :
			continue
		ntopdaus = 0
		for i in range(p.numberOfDaughters()) :
			if fabs(p.daughter(i).pdgId()) == TOP_ID :
				ntopdaus+=1
		if ntopdaus == 2 and p.pdgId() > 0 :
			factor = 0.0
			if p.eta() > 0 :
				factor = 1.0
			elif p.eta() < 0 :
				factor = -1.0
			return ROOT.TLorentzVector(1.0,0.0,factor*sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY)
	#return ROOT.TLorentzVector(1.0,0.0,sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY) #DEBUG RETURN

#Pythia8 findInitialQuark function
def findInitialQuarkPythia8(GenParticles) :
	return ROOT.TLorentzVector(1.0,0.0,sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY) #DEBUG RETURN