#Muon class

#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JUNE 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

import ROOT
from jet import jet
from lepHelper import findNearestJet

class muon :

	DR_MIN = 0.5
	REL_PT_MIN = 25. #GeV
	MUON_PT_CUT = 40. #GeV
	MUON_ETA_CUT = 2.1

	#init function
	def __init__(self,fourvec,charge,tight,loose,jets_list) :
		self.vec = ROOT.TLorentzVector(fourvec)
		self.charge = charge
		self.isTight = tight
		self.isLoose = loose
		if fourvec.M()<=0 or charge==0 or (tight!=1 and tight!=0) or (loose!=1 and loose!=0) :
			print '		WARNING, VERY STRANGE MUON ADDED!'
			print '		muon = (%.2f,%.2f,%.2f,%.2f), charge = %d, loose = %d, tight = %d'%(fourvec.Pt(),fourvec.Eta(),
				fourvec.Phi(),fourvec.M(),charge,loose,tight)
		nearestJetVec = findNearestJet(self.vec,jets_list)
		self.relPt = nearestJetVec.Pt(self.vec.Vect())
		self.dR = nearestJetVec.DeltaR(self.vec)


