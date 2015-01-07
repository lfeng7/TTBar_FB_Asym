import ROOT

#eventTypeHelper: holds helper functions for doing the event type split 
#and for setting the fourvector of the initial state quark if the file 
#is a Monte Carlo file

#global variables
#Beam energy
SQRT_S=13000.0
BEAM_ENERGY=SQRT_S/2.0

#eventTypeCheck
#takes in: GenParticle tree, event type indicator
#returns: boolean value of whether event type matches
def eventTypeCheck(GenParticles,event_type) :
	#event types:
	#	0 = semileptonic qqbar
	#	1 = semileptonic gg, etc.
	#	2 = dileptonic
	#	3 = hadronic

	return True #DEBUG RETURN



#findInitialQuark
#takes in GenParticle tree
#returns TLorentzVector of initial state QUARK based on MC truth information
#NOTE: this just gives one direction or the other down the beam pipe
def findInitialQuark(GenParticles) :
	return ROOT.TLorentzVector(0.0,0.0,sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY) #DEBUG RETURN

#findMCTops
#takes in GenParticle tree
#returns a tuple of the (top, antitop) fourvectors from Monte Carlo
def findMCTops(GenParticles) :
	return (ROOT.TLorentzVector(1.0,0.0,0.0,1.0),ROOT.TLorentzVector(1.0,0.0,0.0,1.0)) #DEBUG RETURN