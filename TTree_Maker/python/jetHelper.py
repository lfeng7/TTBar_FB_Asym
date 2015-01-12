import ROOT
from math import *

#jetHelper contains helper functions for jet selection

#Global variables
#cutflow
CUTFLOW_EXACTLY_ONE_LEPTONIC_BCAND = 4
CUTFLOW_EXACTLY_ONE_HADRONIC_WCAND = 5
CUTFLOW_EXACTLY_ONE_HADRONIC_BCAND = 6
CUTFLOW_EXACTLY_ONE_HADRONIC_TOPCAND = 7
#cut values
MIN_LEPB_PT = 25 #GeV
CSVL_WORKING_POINT = 0.244
MIN_LEP_TOP_MASS_WINDOW = 140 #GeV
MAX_LEP_TOP_MASS_WINDOW = 250 #GeV
MIN_HAD_TOP_PT = 400 #GeV
MIN_HAD_TOP_MASS = 130 #GeV
MAX_HAD_TOP_TAU32 = 0.7
MIN_HAD_W_PT = 200 #GeV
MIN_HAD_W_MASS_WINDOW = 70 #GeV
MAX_HAD_W_MASS_WINDOW = 100 #GeV
MAX_HAD_W_TAU21 = 0.75
MIN_HAD_B_PT = 25 #GeV
MIN_HAD_TOP_MASS_WINDOW = 130 #GeV
MAX_HAD_TOP_MASS_WINDOW = 270 #GeV
MIN_HAD_W_B_DELTAR = 0.6

#selectJets takes in the top type (1 or 2), the lepton vector, the met vectors, and lists of AK4 and AK8 jet variables
#returns:
#	1) a 2 or 3 element list of tuples in the ordering
#		-leptonic b (fourvector, CSV value), hadronic top (fourvector, CSV value, tau1, tau2, tau3) for type 1 or
#		-leptonic b (fourvector, CSV value), hadronic b (fourvector, CSV value), hadronic W (fourvector, CSV value, tau1, tau2, tau3) for type 2 OR
#	2) a 1-element list of negative value of a cutflow failure point
def selectJets(toptype,lepvec,met1vec,met2vec,jetvars_AK4,jetvars_AK8) :
	return_list = []
	#first find the leptonic b by looking for loosely b-tagged jets in the lepton hemisphere
	lepbcands = []
	for i in range(len(jetvars_AK4[0])) :
		#min pT cut and btag requirement
		if jetvars_AK4[0][i] < MIN_LEPB_PT or jetvars_AK4[4][i] < CSVL_WORKING_POINT :
			continue
		thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		thisJet.SetPtEtaPhiM(jetvars_AK4[0][i],jetvars_AK4[1][i],jetvars_AK4[2][i],jetvars_AK4[3][i])
		#check that it's in the leptonic hemisphere
		if lepvec.DeltaR(thisJet) > pi/2. :
			continue
		#check that it combines with the lepton and met vectors to give a mass within the top mass window
		if (lepvec+met1vec+thisJet).M() in range(MIN_LEP_TOP_MASS_WINDOW,MAX_LEP_TOP_MASS_WINDOW) or (lepvec+met2vec+thisJet).M() in range(MIN_LEP_TOP_MASS_WINDOW,MAX_LEP_TOP_MASS_WINDOW) :
			lepbcands.append((thisJet,jetvars_AK4[4][i]))
	#check that only one jet satisfies the selection criteria
	if len(lepbcands) != 1 :
		return [-1*CUTFLOW_EXACTLY_ONE_LEPTONIC_BCAND]
	#add the b candidate to the return list
	return_list.append(lepbcands[0])

	#now get the one or two jets making up the hadronic top 
	if toptype == 1 :
		return_list.append(selectJetsType1Tops(lepvec,jetvars_AK8))
	elif toptype == 2 :
		twojetstuple = selectJetsType2Tops(lepvec,jetvars_AK4,jetvars_AK8)
		return_list.append(twojetstuple[0]); return_list.append(twojetstuple[1])

	#check if the hadronic top reconstruction failed and return the cutflow value
	if len(return_list[1]) == 1 :
		return return_list[1][0]

	return return_list

#selectJetsType1Tops does the selection for type 1 (fully merged) tops
def selectJetsType1Tops(lepvec,jetvars_AK8) :
	#look through all the AK8 jets
	hadtopcands = []
	for i in range(len(jetvars_AK8[0])) :
		#hard cuts on pT, mass, and tau3/tau2
		if jetvars_AK8[0][i] < MIN_HAD_TOP_PT or jetvars_AK8[3][i] < MIN_HAD_TOP_MASS or jetvars_AK8[7][i]/jetvars_AK8[6][i] > MAX_HAD_TOP_TAU32 :
			continue
		#check that it's on the hadronic hemisphere
		thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		thisJet.SetPtEtaPhiM(jetvars_AK8[0][i],jetvars_AK8[1][i],jetvars_AK8[2][i],jetvars_AK8[3][i])
		if lepvec.DeltaR(thisJet) < pi/2. :
			continue
		hadtopcands.append((thisJet,jetvars_AK8[4][i],jetvars_AK8[5][i],jetvars_AK8[6][i],jetvars_AK8[7][i]))
	#require exactly one hadronic top candidate
	if len(hadtopcands) != 1 :
		return (-1*CUTFLOW_EXACTLY_ONE_HADRONIC_TOPCAND)
	#return the hadronic top candidate
	return hadtopcands[0]
	#return (ROOT.TLorentzVector(1.0,0.0,0.0,1.0),0.5,0.5,0.5,0.5) #DEBUG RETURN

#selectJetsType2Tops does the selection for type 2 (partially merged) tops
def  selectJetsType2Tops(lepvec,jetvars_AK4,jetvars_AK8) :
	#look through the AK8 jets for W candidates
	hadWcands = []
	for i in range(len(jetvars_AK8[0])) :
		#hard cuts on pT and tau2/tau1
		if jetvars_AK8[0][i] < MIN_HAD_W_PT or jetvars_AK8[6][i]/jetvars_AK8[5][i] > MAX_HAD_W_TAU21 :
			continue
		#check that it's in the W mass window
		thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		thisJet.SetPtEtaPhiM(jetvars_AK8[0][i],jetvars_AK8[1][i],jetvars_AK8[2][i],jetvars_AK8[3][i])
		if thisJet.M() not in range(MIN_HAD_W_MASS_WINDOW,MAX_HAD_W_MASS_WINDOW) :
			continue
		#check that it's in the hadronic hemisphere
		if lepvec.DeltaR(thisJet) < pi/2. :
			continue
		#append to the list of W candidates
		hadWcands.append((thisJet,jetvars_AK8[4][i],jetvars_AK8[5][i],jetvars_AK8[6][i],jetvars_AK8[7][i]))
	#require exactly one hadronic W candidate
	if len(hadWcands) != 1 :
		return (-1*CUTFLOW_EXACTLY_ONE_HADRONIC_WCAND)

	#now look through the AK4 jets for hadronic b candidates
	hadbcands = []
	for i in range(len(jetvars_AK4[0])) :
		#hard cut on pT 
		if jetvars_AK4[0][i] < MIN_HAD_B_PT :
			continue
		#check that it combines with the W to give a mass in the top mass window
		thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		thisJet.SetPtEtaPhiM(jetvars_AK4[0][i],jetvars_AK4[1][i],jetvars_AK4[2][i],jetvars_AK4[3][i])
		if (thisJet+hadWcands[0][0]).M() not in range(MIN_HAD_TOP_MASS_WINDOW,MAX_HAD_TOP_MASS_WINDOW) :
			continue
		#check that it's in the hadronic hemisphere
		if lepvec.DeltaR(thisJet) < pi/2. :
			continue
		#check that it's sufficiently separated from the W jet
		if thisJet.DeltaR(hadWcands[0][0]) < MIN_HAD_W_B_DELTAR :
			continue
		#append to the list of b candidates
		hadbcands.append((thisJet,jetvars_AK4[4][i]))
	#require exactly one hadronic b candidate
	if len(hadbcands) != 1 :
		return (-1*CUTFLOW_EXACTLY_ONE_HADRONIC_BCAND)

	return (hadWcands[0],hadbcands[0])
	#return [(ROOT.TLorentzVector(1.0,0.0,0.0,1.0),0.5),(ROOT.TLorentzVector(1.0,0.0,0.0,1.0),0.5,0.5,0.5,0.5)] #DEBUG RETURN