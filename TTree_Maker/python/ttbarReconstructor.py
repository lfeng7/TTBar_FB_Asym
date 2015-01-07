import ROOT
from math import *

#ttbar Reconstructor contains helper functions to reconstruct the ttbar system

#reconstruct takes in the lepton and met fourvectors and the list of jet tuples
#	with fourvectors and CSV values
#returns a 3-tuple of:
#	1) the corrected lepton fourvector
#	2) the corrected met fourvector
#	3) a list of tuples of (jet fourvectors, jet CSVs) in the order:
#	   leptonic b, hadronic W subjet 1, hadronic W subjet 2, hadronic b, extra jets
#	4) the final Chi2 value from the kinematic fit
def reconstruct(lepton,met,jettuples) :
	return ( (ROOT.TLorentzVector(1.0,0.0,0.0,1.0),ROOT.TLorentzVector(1.0,0.0,0.0,1.0),
				[(ROOT.TLorentzVector(1.0,0.0,0.0,1.0),0.5),(ROOT.TLorentzVector(1.0,0.0,0.0,1.0),0.5),
				 (ROOT.TLorentzVector(1.0,0.0,0.0,1.0),0.5),(ROOT.TLorentzVector(1.0,0.0,0.0,1.0),0.5)],55.) ) #DEBUG RETURN