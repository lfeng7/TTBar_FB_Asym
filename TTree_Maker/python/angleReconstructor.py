import ROOT

#angleReconstructor calculates differential cross section observables given event fourvectors

#getObservables takes in the (naive) initial parton and reconstructed top quark vectors 
#returns a 3-tuple of:
#	1) costheta
#	2) feynman x
#	3) ttbar mass
def getObservables(q_vec,qbar_vec,lept_vec,hadt_vec,lepton_charge) :
	return (0.0,0.2,450.) #DEBUG RETURN

#getMCObservables takes in the MC TRUTH initial parton, t, and tbar fourvectors
#returns a 13-tuple of:
#	1-3) MC truth costheta, feynman x, and ttbar mass
#	4-8) antisymmetric, symmetric/antisymmetric xi, and symmetric/antisymmetic delta reweighting factors
#	9-13) the same reweighting factors calculated with the opposite sign angle
def getMCObservables(q_vec,qbar_vec,t_vec,tbar_vec) :
	return (0.0,0.2,450.,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0) #DEBUG RETURN