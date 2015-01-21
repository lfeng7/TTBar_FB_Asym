#jetHelper contains helper functions for jet selection
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

import ROOT
from math import *

#Global variables
#constants
TOP_MASS = 172.5
#cutflow
CUTFLOW_EXACTLY_ONE_LEPTONIC_BCAND = 6
CUTFLOW_EXACTLY_ONE_HADRONIC_WCAND = 7
CUTFLOW_EXACTLY_ONE_HADRONIC_BCAND = 8
CUTFLOW_EXACTLY_ONE_HADRONIC_TOPCAND = 9
#cut values
MIN_LEPB_PT = 25. #GeV
CSVL_WORKING_POINT = 0.244
MIN_LEP_TOP_MASS_WINDOW = 140. #GeV
MAX_LEP_TOP_MASS_WINDOW = 250. #GeV
MIN_HAD_TOP_PT = 400. #GeV
MIN_HAD_TOP_MASS = 130. #GeV
MAX_HAD_TOP_TAU32 = 0.7
MIN_HAD_W_PT = 200. #GeV
MIN_HAD_W_MASS_WINDOW = 70. #GeV
MAX_HAD_W_MASS_WINDOW = 100. #GeV
MAX_HAD_W_TAU21 = 0.75
MIN_HAD_B_PT = 25. #GeV
MIN_HAD_TOP_MASS_WINDOW = 130. #GeV
MAX_HAD_TOP_MASS_WINDOW = 270. #GeV
MIN_HAD_W_B_DELTAR = 0.6

#selectJets takes in the top type (1 or 2), the lepton vector, the met vectors, and lists of AK4 and AK8 jet variables
#returns:
#	1) a 2 or 3 element list of tuples in the ordering
#		-leptonic b (fourvector, CSV value), hadronic top (fourvector, CSV value, tau1, tau2, tau3) for type 1 or
#		-leptonic b (fourvector, CSV value), hadronic b (fourvector, CSV value), 
#			hadronic W (fourvector, CSV value, tau1, tau2, tau3) for type 2 OR
#	2) a 1-element list of negative value of a cutflow failure point
def selectJets(toptype,lepvec,met1vec,met2vec,jetvars_AK4,jetvars_AK8,jet_control_plots) :
	return_list = []
	#first find the leptonic b by looking for loosely b-tagged jets in the lepton hemisphere
	lepbcands = []
	for i in range(len(jetvars_AK4[0])) :
		#min pT cut and btag requirement
		jet_control_plots[0].Fill(jetvars_AK4[4][i]); jet_control_plots[1].Fill(jetvars_AK4[0][i])
		if jetvars_AK4[0][i] < MIN_LEPB_PT or jetvars_AK4[4][i] < CSVL_WORKING_POINT :
			continue
		thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		thisJet.SetPtEtaPhiM(jetvars_AK4[0][i],jetvars_AK4[1][i],jetvars_AK4[2][i],jetvars_AK4[3][i])
		#check that it's in the leptonic hemisphere
		lep_dR = lepvec.DeltaR(thisJet)
		jet_control_plots[2].Fill(lep_dR)
		if lep_dR > pi/2. :
			continue
		#check that it combines with the lepton and met vectors to give a mass within the top mass window
		comb_mass_1 = (lepvec+met1vec+thisJet).M(); comb_mass_2 = (lepvec+met2vec+thisJet).M()
		if fabs(comb_mass_1-TOP_MASS) < fabs(comb_mass_2-TOP_MASS) :
			jet_control_plots[3].Fill(comb_mass_1)
		else :
			jet_control_plots[3].Fill(comb_mass_2)
		if ( (comb_mass_1 > MIN_LEP_TOP_MASS_WINDOW and comb_mass_1 < MAX_LEP_TOP_MASS_WINDOW) or 
			(comb_mass_2 > MIN_LEP_TOP_MASS_WINDOW and comb_mass_2 < MAX_LEP_TOP_MASS_WINDOW) ) :
			lepbcands.append((thisJet,jetvars_AK4[4][i]))
	#check that only one jet satisfies the selection criteria
	jet_control_plots[4].Fill(len(lepbcands))
	if len(lepbcands) != 1 :
		return [-1*CUTFLOW_EXACTLY_ONE_LEPTONIC_BCAND]
	#add the b candidate to the return list
	return_list.append(lepbcands[0])

	#now get the one or two jets making up the hadronic top 
#	print 'toptype = '+str(toptype)+'' #DEBUGGING
	if toptype == 1 :
		return_list.append(selectJetsType1Tops(lepvec,jetvars_AK8,jet_control_plots))
	elif toptype == 2 :
		twojetstuple = selectJetsType2Tops(lepvec,jetvars_AK4,jetvars_AK8,jet_control_plots)
		for i in range(len(twojetstuple)) :
			return_list.append(twojetstuple[i])

	#check if the hadronic top reconstruction failed and return the cutflow value
#	print 'return_list = '+str(return_list)+'' #DEBUGGING
	if len(return_list[1]) == 1 :
		return return_list[1]

	return return_list

#selectJetsType1Tops does the selection for type 1 (fully merged) tops
def selectJetsType1Tops(lepvec,jetvars_AK8,jet_control_plots) :
	#look through all the AK8 jets
	hadtopcands = []
	for i in range(len(jetvars_AK8[0])) :
		#hard cuts on pT, mass, and tau3/tau2
		jet_control_plots[5].Fill(jetvars_AK8[0][i]) 
		jet_control_plots[6].Fill(jetvars_AK8[3][i]) 
		jet_control_plots[7].Fill(jetvars_AK8[7][i]/jetvars_AK8[6][i])
		if ( jetvars_AK8[0][i] < MIN_HAD_TOP_PT or jetvars_AK8[3][i] < MIN_HAD_TOP_MASS or 
			jetvars_AK8[7][i]/jetvars_AK8[6][i] > MAX_HAD_TOP_TAU32 ) :
			continue
		#check that it's on the hadronic hemisphere
		thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		thisJet.SetPtEtaPhiM(jetvars_AK8[0][i],jetvars_AK8[1][i],jetvars_AK8[2][i],jetvars_AK8[3][i])
		jet_control_plots[8].Fill(lepvec.DeltaR(thisJet))
		if lepvec.DeltaR(thisJet) < pi/2. :
			continue
		hadtopcands.append((thisJet,jetvars_AK8[4][i],jetvars_AK8[5][i],jetvars_AK8[6][i],jetvars_AK8[7][i]))
#	print 'len(hadtopcands) = '+str(len(hadtopcands))+'' #DEBUGGING
	#require exactly one hadronic top candidate
	jet_control_plots[9].Fill(len(hadtopcands))
	if len(hadtopcands) != 1 :
		return [-1*CUTFLOW_EXACTLY_ONE_HADRONIC_TOPCAND]
	#return the hadronic top candidate
	return hadtopcands[0]
	#return (ROOT.TLorentzVector(1.0,0.0,0.0,1.0),0.5,0.5,0.5,0.5) #DEBUG RETURN

#selectJetsType2Tops does the selection for type 2 (partially merged) tops
def  selectJetsType2Tops(lepvec,jetvars_AK4,jetvars_AK8,jet_control_plots) :
	#look through the AK8 jets for W candidates
	hadWcands = []
	for i in range(len(jetvars_AK8[0])) :
		#hard cuts on pT and tau2/tau1
		jet_control_plots[10].Fill(jetvars_AK8[0][i]); jet_control_plots[11].Fill(jetvars_AK8[6][i]/jetvars_AK8[5][i])
		if jetvars_AK8[0][i] < MIN_HAD_W_PT or jetvars_AK8[6][i]/jetvars_AK8[5][i] > MAX_HAD_W_TAU21 :
			continue
		#check that it's in the W mass window
		thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		thisJet.SetPtEtaPhiM(jetvars_AK8[0][i],jetvars_AK8[1][i],jetvars_AK8[2][i],jetvars_AK8[3][i])
		jetMass = thisJet.M()
		jet_control_plots[12].Fill(jetMass)
		if jetMass < MIN_HAD_W_MASS_WINDOW or jetMass > MAX_HAD_W_MASS_WINDOW :
			continue
		#check that it's in the hadronic hemisphere
		jet_control_plots[13].Fill(lepvec.DeltaR(thisJet))
		if lepvec.DeltaR(thisJet) < pi/2. :
			continue
		#append to the list of W candidates
		hadWcands.append((thisJet,jetvars_AK8[4][i],jetvars_AK8[5][i],jetvars_AK8[6][i],jetvars_AK8[7][i]))
	#require exactly one hadronic W candidate
	jet_control_plots[14].Fill(len(hadWcands))
	if len(hadWcands) != 1 :
		return ([-1*CUTFLOW_EXACTLY_ONE_HADRONIC_WCAND],[-1*CUTFLOW_EXACTLY_ONE_HADRONIC_WCAND])

	#now look through the AK4 jets for hadronic b candidates
	hadbcands = []
	for i in range(len(jetvars_AK4[0])) :
		#hard cut on pT 
		if jetvars_AK4[0][i] < MIN_HAD_B_PT :
			continue
		#check that it combines with the W to give a mass in the top mass window
		thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		thisJet.SetPtEtaPhiM(jetvars_AK4[0][i],jetvars_AK4[1][i],jetvars_AK4[2][i],jetvars_AK4[3][i])
		combMass = (thisJet+hadWcands[0][0]).M()
		jet_control_plots[15].Fill(combMass)
		if combMass < MIN_HAD_TOP_MASS_WINDOW or combMass > MAX_HAD_TOP_MASS_WINDOW :
			continue
		#check that it's in the hadronic hemisphere
		jet_control_plots[16].Fill(lepvec.DeltaR(thisJet))
		if lepvec.DeltaR(thisJet) < pi/2. :
			continue
		#check that it's sufficiently separated from the W jet
		jet_control_plots[17].Fill(thisJet.DeltaR(hadWcands[0][0]))
		if thisJet.DeltaR(hadWcands[0][0]) < MIN_HAD_W_B_DELTAR :
			continue
		#append to the list of b candidates
		hadbcands.append((thisJet,jetvars_AK4[4][i]))
	#require exactly one hadronic b candidate
	jet_control_plots[18].Fill(len(hadbcands))
	if len(hadbcands) != 1 :
		return ([-1*CUTFLOW_EXACTLY_ONE_HADRONIC_BCAND],[-1*CUTFLOW_EXACTLY_ONE_HADRONIC_BCAND])

	return (hadWcands[0],hadbcands[0])
	#return [(ROOT.TLorentzVector(1.0,0.0,0.0,1.0),0.5),(ROOT.TLorentzVector(1.0,0.0,0.0,1.0),0.5,0.5,0.5,0.5)] #DEBUG 