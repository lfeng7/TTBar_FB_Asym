#eventTypeHelper: holds helper functions for doing the event type split 
#and for setting the fourvector of the initial state quark if the file 
#is a Monte Carlo file
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

import ROOT
from math import *

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
def eventTypeCheck(generator,GenParticles,genPartVars,event_type) :
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
		return typeCheckPythia8(GenParticles,genPartVars,event_type)
	else :
		print 'ERROR: GENERATOR '+generator+' NOT RECOGNIZED!!!'
		return (False,0)

#findInitialQuark
#takes in GenParticle tree
#returns TLorentzVector of initial state QUARK based on MC truth information
#NOTE: this just gives one direction or the other down the beam pipe
def findInitialQuark(generator,GenParticles,genPartVars) :
	if generator == 'powheg' or generator == 'madgraph' or generator == 'mg5' or generator == 'mg' :
		return findInitialQuarkPowheg(GenParticles)
	elif generator == 'mcatnlo' :
		return findInitialQuarkMCAtNLO(GenParticles)
	elif generator == 'pythia8' :
		return findInitialQuarkPythia8(GenParticles,genPartVars)
	else :
		print 'ERROR: GENERATOR '+generator+' NOT RECOGNIZED!!!'
		return ROOT.TLorentzVector(1.0,0.0,sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY) #DEBUG RETURN

#findMCTops
#takes in GenParticle tree
#returns a tuple of the (top, antitop) fourvectors from Monte Carlo
def findMCTops(generator,GenParticles) :
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
			if generator == 'pythia8' : #crank down the particle if it's pythia8
				while p.daughter(0).pdgId() == p.pdgId() :
					p = p.daughter(0)
			tvec.SetPtEtaPhiM(p.pt(),p.eta(),p.phi(),p.mass())
			found_t = True
		elif p.pdgId() == -1*TOP_ID and not found_tbar :
			if generator == 'pythia8' : #crank down the particle if it's pythia8
				while p.daughter(0).pdgId() == p.pdgId() :
					p = p.daughter(0)
			tbarvec.SetPtEtaPhiM(p.pt(),p.eta(),p.phi(),p.mass())
			found_tbar = True
	return (tvec,tbarvec)
	#return (ROOT.TLorentzVector(1.0,0.0,0.0,1.0),ROOT.TLorentzVector(1.0,0.0,0.0,1.0)) #DEBUG RETURN

##################################  HELPER FUNCTIONS  ##################################

#Powheg eventTypeCheck function
def typeCheckPowheg(GenParticles,event_type) :
	addTwice = False
	#check the initial state partons if necessary
	if event_type<2 :
		initial_parton_ids = []
		#loop through and find the pdg IDs of the first daughters of the protons
		for p in GenParticles :
			if p.pt()<0 or p.status()!=3 :
				continue
			if p.pdgId()==PROTON_ID :
				initial_parton_ids.append(p.daughter(0).pdgId())
		#is it a qqbar event?
		is_qq = len(initial_parton_ids) == 2 and initial_parton_ids[0]+initial_parton_ids[1]==0
		#regardless, did it have an initially symmetric state
		addTwice = event_type==0 or (len(initial_parton_ids)==2 and initial_parton_ids[0] == initial_parton_ids[1])
		#return false if the event type is incorrect
		if (event_type == 0 and not is_qq) or (event_type == 1 and is_qq) :
			return False,addTwice
	return semilepCheck(GenParticles,event_type),addTwice
	#return True #DEBUG RETURN

#MC@NLO eventTypeCheck function
def typeCheckMCAtNLO(GenParticles,event_type) :
	addTwice = False
	#check the initial state partons if necessary
	if event_type<2 :
		initial_parton_ids = []
		#loop through and find the particles whose daughters include the ttbar pair
		for p in GenParticles :
			if p.pt()<0 or p.status()!=3 :
				continue
			ntopDaus = 0
			for i in range(p.numberOfDaughters()) :
				if fabs(p.daughter(i).pdgId()) == TOP_ID :
					ntopDaus+=1
			if ntopDaus == 2 :
				initial_parton_ids.append(p.pdgId())
		#is it a qqbar event?
		is_qq = len(initial_parton_ids) == 2 and initial_parton_ids[0]+initial_parton_ids[1]==0
		#regardless, did it have an initially symmetric state
		addTwice = event_type==0 or (len(initial_parton_ids)==2 and initial_parton_ids[0] == initial_parton_ids[1])
		#return false if the event type is incorrect
		if (event_type == 0 and not is_qq) or (event_type == 1 and is_qq) :
			return False,addTwice
	return semilepCheck(GenParticles,event_type),addTwice
	#return True #DEBUG RETURN

#pythia8 eventTypeCheck function
def typeCheckPythia8(GenParticles,genPartVars,event_type) :
#	print 'new event---------------------------------' #DEBUGGING
	addTwice = False
	#check the initial state partons if necessary
	if event_type<2 :
		initial_parton_ids = []
		for i in range(len(genPartVars[0])) :
			if genPartVars[6][i] == 21 and genPartVars[5][i] == PROTON_ID :
				initial_parton_ids.append(genPartVars[4][i])
		#is it a qqbar event?
		is_qq = len(initial_parton_ids) == 2 and initial_parton_ids[0]+initial_parton_ids[1]==0
		#regardless, did it have an initially symmetric state
		addTwice = True #yes it did because it's raw pythia8 which only has qqbar or gg events in it
		if (event_type == 0 and not is_qq) or (event_type == 1 and is_qq) :
			return (False,addTwice)

	#check the decay type
	n_leps_from_Ws = 0
	for p in GenParticles :
		#look only for valid daughters of Ws
		if ( p.pt()<0 or (p.status()!=23 and p.status()!=1 and p.status()!=2) or 
			p.numberOfMothers()!=1 or fabs(p.mother(0).pdgId()) != W_ID ) :
			continue
		#look for leptons
		if fabs(p.pdgId()) in range(ELECTRON_ID,TAU_NEUTRINO_ID+1) :
			n_leps_from_Ws+=1
#	print '	n_leps_from_Ws = '+str(n_leps_from_Ws)+'' #DEBUGGING
	if ( (n_leps_from_Ws == 2 and event_type > 1) or (n_leps_from_Ws == 4 and event_type != 2) or 
		(n_leps_from_Ws == 0 and event_type != 3) ) :
		return (False,addTwice)
	
	return (True,addTwice) #DEBUG RETURN

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
	if ( (len(leptons_from_Ws_IDs) == 2 and event_type > 1) or (len(leptons_from_Ws_IDs) == 4 and event_type != 2) or 
		(len(leptons_from_Ws_IDs) == 0 and event_type != 3) ) :
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
def findInitialQuarkPythia8(GenParticles,genPartVars) :
	for i in range(len(genPartVars[0])) :
		if genPartVars[6][i] == 21 and genPartVars[5][i] == PROTON_ID and genPartVars[4][i] > 0 :
			factor = 0.0
			if genPartVars[1][i] > 0 :
				factor = 1.0
			elif genPartVars[1][i] < 0 :
				factor = -1.0
			return ROOT.TLorentzVector(1.0,0.0,factor*sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY)
	return ROOT.TLorentzVector(1.0,0.0,sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY) #DEBUG RETURN