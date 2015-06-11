#Jet Class

#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JUNE 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

import ROOT
from jetHelper import matchUnprunedVec, getSF

class jet :

	CSVL_WORKING_POINT = 0.244

	#init function
	def __init__(self,fourvec,unpruned_fourvecs,tau1s,tau2s,tau3s,csv_value,flav) :
		self.vec = ROOT.TLorentzVector(fourvec.X(),fourvec.Y(),fourvec.Z(),fourvec.T())
		self.csv = csv_value
		self.flavor = flav
		self.btagSF, self.btagSFlow, self.btagSFhigh = getSF(self.vec, self.flavor, self.CSVL_WORKING_POINT)
		matchedUnprunedJetIndex = matchUnprunedVec(self.vec,unpruned_fourvecs)
		if matchedUnprunedJetIndex == -1 :
			self.matchedUnprunedVec = self.vec
			self.tau32 = -1.
			self.tau21 = -1.
		else :
			thisVec = unpruned_fourvecs[matchedUnprunedJetIndex]
			self.matchedUnprunedVec = ROOT.TLorentzVector(thisVec.X(),thisVec.Y(),thisVec.Z(),thisVec.T())
			tau1 = tau1s[matchedUnprunedJetIndex]
			tau2 = tau2s[matchedUnprunedJetIndex]
			tau3 = tau3s[matchedUnprunedJetIndex]
			if tau2!=0. :
				self.tau32 = tau3/tau2
			else :
				self.tau32 = -1.
			if tau1!=0. :
				self.tau21 = tau2/tau1
			else :
				self.tau21 = -1.