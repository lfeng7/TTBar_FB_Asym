import ROOT

#metHelper: contains helper functions for dealing with the MET/neutrino

#setupMET
#takes in the fourvector of the selected lepton and the list of 
#variables describing the MET
#returns the fourvector of the neutrino assuming just the WMass constraint
def setupMET(lep_vec,metVars) :
	return ROOT.TLorentzVector(1.0,0.0,0.0,1.0) #DEBUG RETURN