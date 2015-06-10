#Electron class

#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JUNE 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

import ROOT
from jet import jet
from lepHelper import findNearestJet

class electron :

	DR_MIN = 0.5
	REL_PT_MIN = 25. #GeV
	TRIANGLE_CUT_RATIO = (1.5/75.)
	ELECTRON_PT_CUT = 25. #GeV
	ELECTRON_ETA_CUT = 2.4

	#init function
	def __init__(self,fourvec,charge,tight,loose,metvec,jets_list) :
		self.vec = ROOT.TLorentzVector(fourvec)
		self.charge = charge
		self.isTight = tight
		self.isLoose = loose
		if fourvec.M()<=0 or charge==0 or (tight!=1 and tight!=0) or (loose!=1 and loose!=0) :
			print '		WARNING, VERY STRANGE ELECTRON ADDED!'
			print '		electron = (%.2f,%.2f,%.2f,%.2f), charge = %d, loose = %d, tight = %d'%(fourvec.Pt(),
				fourvec.Eta(),fourvec.Phi(),fourvec.M(),charge,loose,tight)
		nearestJetVec = findNearestJet(self.vec,jets_list)
		self.relPt = nearestJetVec.Pt(self.vec.Vect())
		self.dR = nearestJetVec.DeltaR(self.vec)
		self.triangle_el_val  = abs(self.vec.DeltaPhi(metvec)-1.5)
		self.triangle_jet_val = abs(jets_list[0].vec.DeltaPhi(metvec)-1.5)
		self.triangle_cut_val = TRIANGLE_CUT_RATIO*metvec.E()