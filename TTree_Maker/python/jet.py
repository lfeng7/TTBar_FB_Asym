#Jet Class

#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JUNE 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

import ROOT
from jetHelper import matchedUnprunedVec, getSF

class jet :

	CSVL_WORKING_POINT = 0.244

	#init function
	def __init__(self,fourvec,unpruned_fourvecs,tau1s,tau2s,tau3s,csv_value,flav) :
		self.vec = ROOT.TLorentzVector(fourvec)
		matchedUnprunedJetIndex = matchUnprunedVec(self.vec,unpruned_fourvecs)
		self.matchedUnprunedVec = unpruned_fourvecs[matchedUnprunedJetIndex]
		self.tau32 = tau3s[matchedUnprunedJetIndex]/tau2s[matchedUnprunedJetIndex]
		self.tau21 = tau3s[matchedUnprunedJetIndex]/tau2s[matchedUnprunedJetIndex]
		self.csv = csv_value
		self.flavor = flav
		self.btagSF, self.btagSFlow, self.btagSFhigh = getSF(self.vec, self.csv, self.flavor, CSVL_WORKING_POINT)