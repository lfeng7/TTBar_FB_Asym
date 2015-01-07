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

#eventTypeCheck
#takes in: GenParticle tree, event type indicator
#returns: boolean value of whether event type matches
def eventTypeCheck(GenParticles,event_type) :
	#event types:
	#	0 = semileptonic qqbar
	#	1 = semileptonic gg, etc.
	#	2 = dileptonic
	#	3 = hadronic

	#check the initial state partons if necessary
	if event_type<2 :
		initial_state_parton_ids = []
		#loop through and find the pdg IDs of the first daughters of the protons
		for particle in GenParticles :
			if particle.pt()<0 :
				continue
			if particle.pdgId()==PROTON_ID :
				s = 'proton found! Daughters: '
				for i in range(particle.numberOfDaughters()) :
					s = s+str(particle.daughter(i).pdgId())+' '
				print s
				#initial_state_parton_ids.append(particle.daughter(0).pdgId())
		#is it a qqbar event?
		is_qq = len(initial_state_parton_ids) == 2 and initial_state_parton_ids[0]+initial_state_parton_ids[1]==0
		print 'initial_state_parton_ids='+str(initial_state_parton_ids)+', is_qq='+str(is_qq)
		#return false if the event type is incorrect
		if (event_type == 0 and not is_qq) or (event_type == 1 and is_qq) :
			return False

	#check the decay type
	leptons_from_Ws_IDs = []
	for particle in GenParticles :
		if particle.pt()<0 :
			continue
		#look for the Ws that are the decay products of the tops
		if fabs(particle.pdgId())==TOP_ID :
			for i in range(particle.numberOfDaughters()) :
				if fabs(particle.daughter(i).pdgId()) == W_ID :
					w_particle = particle.daughter(i)
					#get the IDs of all its leptonic daughters and add them to the list
					for j in range(w_particle.numberOfDaughters()) :
						if fabs(w_particle.daughter(j).pdgId()) in range(ELECTRON_ID,TAU_NEUTRINO_ID+1) :
							leptons_from_Ws_IDs.append(w_particle.daughter(j).pdgId())
	print 'leptons_from_Ws_IDs='+str(leptons_from_Ws_IDs)
	if (len(leptons_from_Ws_IDs) == 2 and event_type > 1) or (len(leptons_from_Ws_IDs) == 4 and event_type != 2) or (len(leptons_from_Ws_IDs) == 0 and event_type != 3) :
		return False
	return True
	#return True #DEBUG RETURN



#findInitialQuark
#takes in GenParticle tree
#returns TLorentzVector of initial state QUARK based on MC truth information
#NOTE: this just gives one direction or the other down the beam pipe
def findInitialQuark(GenParticles) :
	return ROOT.TLorentzVector(1.0,0.0,sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY) #DEBUG RETURN

#findMCTops
#takes in GenParticle tree
#returns a tuple of the (top, antitop) fourvectors from Monte Carlo
def findMCTops(GenParticles) :
	return (ROOT.TLorentzVector(1.0,0.0,0.0,1.0),ROOT.TLorentzVector(1.0,0.0,0.0,1.0)) #DEBUG RETURN