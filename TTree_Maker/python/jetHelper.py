#jetHelper contains helper functions for jet selection
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

import ROOT
from math import *

#Global variables
#constants
TOP_MASS = 172.5
#cutflow
CUTFLOW_EXACTLY_ONE_LEPTONIC_BCAND = 7
CUTFLOW_EXACTLY_ONE_HADRONIC_WCAND = 8
CUTFLOW_EXACTLY_ONE_HADRONIC_BCAND = 9
CUTFLOW_EXACTLY_ONE_HADRONIC_TOPCAND = 10
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
#	1) a 4 or 5 element list of the following
#		-leptonic b (fourvector, CSV value, flavour), hadronic top (fourvector, CSV value, tau1, tau2, tau3,flavour) 
#			for type 1 or
#		-leptonic b (fourvector, CSV value, flavour), hadronic b (fourvector, CSV value,flavour), 
#			hadronic W (fourvector, CSV value, tau1, tau2, tau3,flavour) for type 2 THEN
#		-the btagging efficiency scale factor, low error, and hi error for MC events OR
#	2) a 1-element list of negative value of a cutflow failure point
def selectJets(isdata,toptype,lepvec,met1vec,met2vec,jetvars_AK4,jetvars_AK8,jet_control_plots) :
	return_list = []
	#first find the leptonic b by looking for loosely b-tagged jets in the lepton hemisphere
	lepbcand = ROOT.TLorentzVector(1.0,0.0,0.0,1.0); lepbcand_CSV = 0.5; lepbcand_flavour = -100
	best_comb_mass_offset = 100000.
#	print 'NEW EVENT LEPTONIC B SELECTION: ' #DEBUGGING
	for i in range(len(jetvars_AK4[0])) :
		#min pT cut
		jet_control_plots[0].Fill(jetvars_AK4[0][i])
		if jetvars_AK4[0][i] < MIN_LEPB_PT :
			continue
		thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		thisJet.SetPtEtaPhiM(jetvars_AK4[0][i],jetvars_AK4[1][i],jetvars_AK4[2][i],jetvars_AK4[3][i])
		#check that it's in the leptonic hemisphere
		lep_dR = lepvec.DeltaR(thisJet)
		jet_control_plots[1].Fill(lep_dR)
		if lep_dR > pi/2. :
			continue
		#check that it combines with the lepton and met vectors to give a mass within the top mass window
		comb_mass_1 = (lepvec+met1vec+thisJet).M(); comb_mass_2 = (lepvec+met2vec+thisJet).M()
		if ( (comb_mass_1 > MIN_LEP_TOP_MASS_WINDOW and comb_mass_1 < MAX_LEP_TOP_MASS_WINDOW) or 
			  (comb_mass_2 > MIN_LEP_TOP_MASS_WINDOW and comb_mass_2 < MAX_LEP_TOP_MASS_WINDOW) ) :
			comb_mass_offset = min(abs(comb_mass_1-TOP_MASS),abs(comb_mass_2-TOP_MASS))
			if comb_mass_offset == abs(comb_mass_1-TOP_MASS) :
				jet_control_plots[2].Fill(comb_mass_1)
			else :
				jet_control_plots[2].Fill(comb_mass_2) 
			if comb_mass_offset<best_comb_mass_offset :
				best_comb_mass_offset = comb_mass_offset
#				print '	best offset = '+str(best_comb_mass_offset)+'' #DEBUGGING
				lepbcand.SetPtEtaPhiM(thisJet.Pt(),thisJet.Eta(),thisJet.Phi(),thisJet.M())
				lepbcand_CSV = jetvars_AK4[4][i]
				lepbcand_flavour = jetvars_AK4[5][i]
	if best_comb_mass_offset == 100000. : #didn't find any jets in the mass window
		return ([-1*CUTFLOW_EXACTLY_ONE_LEPTONIC_BCAND],1.0,1.0,1.0)
	jet_control_plots[3].Fill(lepbcand_CSV)
	#cut that the candidate jet is loosely b-tagged and calculate the scalefactor and error
	if lepbcand_CSV < CSVL_WORKING_POINT :
		return ([-1*CUTFLOW_EXACTLY_ONE_LEPTONIC_BCAND],1.0,1.0,1.0)
	btag_weight = 0.0; btag_weight_low = 0.0; btag_weight_hi = 0.0
	if isdata == 0 :
		btag_weight, btag_weight_low, btag_weight_hi = getSF((lepbcand,lepbcand_CSV,lepbcand_flavour),CSVL_WORKING_POINT)
	if btag_weight == 0.0 :
		btag_weight = 1.0; btag_weight_low = 1.0; btag_weight_hi = 1.0
	#add the b candidate to the return list
	return_list.append((lepbcand,lepbcand_CSV,lepbcand_flavour))
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
		return (return_list[1],1.0,1.0,1.0)

	return (return_list, btag_weight, btag_weight_low, btag_weight_hi)

#selectJetsType1Tops does the selection for type 1 (fully merged) tops
def selectJetsType1Tops(lepvec,jetvars_AK8,jet_control_plots) :
	#look through all the AK8 jets
	hadtopcands = []
	for i in range(len(jetvars_AK8[0])) :
		#hard cuts on pT, mass, and tau3/tau2
		jet_control_plots[4].Fill(jetvars_AK8[0][i]) 
		jet_control_plots[5].Fill(jetvars_AK8[3][i]) 
		if jetvars_AK8[6][i] == 0 :
			continue
		jet_control_plots[6].Fill(jetvars_AK8[7][i]/jetvars_AK8[6][i])
		if ( jetvars_AK8[0][i] < MIN_HAD_TOP_PT or jetvars_AK8[3][i] < MIN_HAD_TOP_MASS or 
			jetvars_AK8[7][i]/jetvars_AK8[6][i] > MAX_HAD_TOP_TAU32 ) :
			continue
		#check that it's on the hadronic hemisphere
		thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		thisJet.SetPtEtaPhiM(jetvars_AK8[0][i],jetvars_AK8[1][i],jetvars_AK8[2][i],jetvars_AK8[3][i])
		jet_control_plots[7].Fill(lepvec.DeltaR(thisJet))
		if lepvec.DeltaR(thisJet) < pi/2. :
			continue
		hadtopcands.append(
			(thisJet,jetvars_AK8[4][i],jetvars_AK8[5][i],jetvars_AK8[6][i],jetvars_AK8[7][i],jetvars_AK8[8][i]) )
#	print 'len(hadtopcands) = '+str(len(hadtopcands))+'' #DEBUGGING
	#require exactly one hadronic top candidate
	jet_control_plots[8].Fill(len(hadtopcands))
	if len(hadtopcands) != 1 :
		return [-1*CUTFLOW_EXACTLY_ONE_HADRONIC_TOPCAND]
	#return the hadronic top candidate
	return hadtopcands[0]
	#return (ROOT.TLorentzVector(1.0,0.0,0.0,1.0),0.5,0.5,0.5,0.5,0) #DEBUG RETURN

#selectJetsType2Tops does the selection for type 2 (partially merged) tops
def  selectJetsType2Tops(lepvec,jetvars_AK4,jetvars_AK8,jet_control_plots) :
	#look through the AK8 jets for W candidates
	hadWcands = []
	for i in range(len(jetvars_AK8[0])) :
		#hard cuts on pT and tau2/tau1
		jet_control_plots[9].Fill(jetvars_AK8[0][i]) 
		if jetvars_AK8[5][i] == 0 :
			continue
		jet_control_plots[10].Fill(jetvars_AK8[6][i]/jetvars_AK8[5][i])
		if jetvars_AK8[0][i] < MIN_HAD_W_PT or jetvars_AK8[6][i]/jetvars_AK8[5][i] > MAX_HAD_W_TAU21 :
			continue
		#check that it's in the W mass window
		thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		thisJet.SetPtEtaPhiM(jetvars_AK8[0][i],jetvars_AK8[1][i],jetvars_AK8[2][i],jetvars_AK8[3][i])
		jetMass = thisJet.M()
		jet_control_plots[11].Fill(jetMass)
		if jetMass < MIN_HAD_W_MASS_WINDOW or jetMass > MAX_HAD_W_MASS_WINDOW :
			continue
		#check that it's in the hadronic hemisphere
		jet_control_plots[12].Fill(lepvec.DeltaR(thisJet))
		if lepvec.DeltaR(thisJet) < pi/2. :
			continue
		#append to the list of W candidates
		hadWcands.append(
			(thisJet,jetvars_AK8[4][i],jetvars_AK8[5][i],jetvars_AK8[6][i],jetvars_AK8[7][i],jetvars_AK8[8][i]) )
	#require exactly one hadronic W candidate
	jet_control_plots[13].Fill(len(hadWcands))
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
		jet_control_plots[14].Fill(combMass)
		if combMass < MIN_HAD_TOP_MASS_WINDOW or combMass > MAX_HAD_TOP_MASS_WINDOW :
			continue
		#check that it's in the hadronic hemisphere
		jet_control_plots[15].Fill(lepvec.DeltaR(thisJet))
		if lepvec.DeltaR(thisJet) < pi/2. :
			continue
		#check that it's sufficiently separated from the W jet
		jet_control_plots[16].Fill(thisJet.DeltaR(hadWcands[0][0]))
		if thisJet.DeltaR(hadWcands[0][0]) < MIN_HAD_W_B_DELTAR :
			continue
		#append to the list of b candidates
		hadbcands.append((thisJet,jetvars_AK4[4][i],jetvars_AK4[5][i]))
	#require exactly one hadronic b candidate
	jet_control_plots[17].Fill(len(hadbcands))
	if len(hadbcands) != 1 :
		return ([-1*CUTFLOW_EXACTLY_ONE_HADRONIC_BCAND],[-1*CUTFLOW_EXACTLY_ONE_HADRONIC_BCAND])

	return (hadWcands[0],hadbcands[0])
	#return [(ROOT.TLorentzVector(1.0,0.0,0.0,1.0),0.5,0),
	#(ROOT.TLorentzVector(1.0,0.0,0.0,1.0),0.5,0.5,0.5,0.5,0)] #DEBUG RETURN

##################################################################################################
##########						b-tagging efficiency stuff 								##########
##################################################################################################

#btagging efficiency stuff
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation53XReReco

#getSF takes a list of jets in the format (4vector,CSV value,flavour) and a CSV working point 
#and returns a tuple of (SF, SF low, SF hi)
def getSF(jet,CSV_wp) :
	if CSV_wp == CSVL_WORKING_POINT :
		sf_func_b = get_CSVL_SFb
		sf_func_c = get_CSVL_SFc
		sf_func_light = get_CSVL_SFlight
	else :
		print 'WARNING: CSV WORKING POINT NOT RECOGNIZED FOR BTAG EFFICIENCY SCALEFACTOR CALCULATION'
		return (1.0,1.0,1.0)
	#get SF, low, hi from functions for jets depending on flavor
	if abs(jet[2]) == 5 : #b jets
		new_SF_tuple = sf_func_b(jet[0].Pt())
	elif abs(jet[2]) == 4 : #c jets
		new_SF_tuple = sf_func_c(jet[0].Pt())
	elif abs(jet[2]) == 0 or abs(jet[2]) == 1 or abs(jet[2]) == 2 or abs(jet[2]) == 3 or abs(jet[2]) == 21 : #light
		new_SF_tuple = sf_func_light(jet[0].Pt(),jet[0].Eta())
	elif abs(jet[2]) == 100 : #skip jets without flavor info
		return (1.0,1.0,1.0)
	else :
		print 'WARNING: UNRECOGNIZED JET FLAVOUR IN BTAG EFFICIENCY SF CALCULATION: flavour = '+str(jet[2])+''
	return_SF = new_SF_tuple[0]
	return_SF_low = new_SF_tuple[1]
	return_SF_hi  = new_SF_tuple[2]
	return (return_SF,return_SF_low,return_SF_hi)
	#return (1.0,1.0,1.0) #DEBUG RETURN

#payload here: https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation53XReReco/SFb-pt_WITHttbar_payload_EPS13.txt
btag_eff_pt_bins = [20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800]
def get_CSVL_SFb(pT) : 
	SF = 1.00572*((1.+(0.013676*pT))/(1.+(0.0143279*pT)))
	CSVL_SFb_error = ( [0.033408,0.015446,0.0146992,0.0183964,0.0185363,0.0145547,0.0176743,0.0203609,0.0143342,
					0.0148771,0.0157936,0.0176496,0.0209156,0.0278529,0.0346877,0.0350101] )
	SF_err = 1.0
	if pT < btag_eff_pt_bins[0] :
		SF_err = 2.0*CSVL_SFb_error[0]
	elif pT > btag_eff_pt_bins[len(btag_eff_pt_bins)-1] :
		SF_err = 2.0*CSVL_SFb_error[len(CSVL_SFb_error)-1]
	else:
		for i in range(len(btag_eff_pt_bins)-1) :
			if pT > btag_eff_pt_bins[i] and pT < btag_eff_pt_bins[i+1] :
				SF_err = CSVL_SFb_error[i] 
	return SF,SF-SF_err,SF+SF_err

def get_CSVL_SFc(pT) : 
	SFb_tuple = get_CSVL_SFb(pT)
	return SFb_tuple[0],SFb_tuple[0]-(2.0*(SFb_tuple[0]-SFb_tuple[1])),SFb_tuple[0]+(2.0*(SFb_tuple[2]-SFb_tuple[0]))

#payload here: https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation53XReReco/SFlightFuncs_EPS2013.C
def get_CSVL_SFlight(pT,eta) : 
	eta = abs(eta)
	SF = 1.0; SF_low = 1.0; SF_hi = 1.0
	if eta > 0.0 and eta < 0.5 :
		SF = ((1.01177+(0.0023066*pT))+(-4.56052e-06*(pT*pT)))+(2.57917e-09*(pT*(pT*pT)))
		SF_low = ((0.977761+(0.00170704*pT))+(-3.2197e-06*(pT*pT)))+(1.78139e-09*(pT*(pT*pT)))
		SF_hi  = ((1.04582+(0.00290226*pT))+(-5.89124e-06*(pT*pT)))+(3.37128e-09*(pT*(pT*pT)))
	elif eta > 0.5 and eta < 1.0 :
		SF = ((0.975966+(0.00196354*pT))+(-3.83768e-06*(pT*pT)))+(2.17466e-09*(pT*(pT*pT)))
		SF_low = ((0.945135+(0.00146006*pT))+(-2.70048e-06*(pT*pT)))+(1.4883e-09*(pT*(pT*pT)))
		SF_hi  = ((1.00683+(0.00246404*pT))+(-4.96729e-06*(pT*pT)))+(2.85697e-09*(pT*(pT*pT)))
	elif eta > 1.0 and eta < 1.5 :
		SF = ((0.93821+(0.00180935*pT))+(-3.86937e-06*(pT*pT)))+(2.43222e-09*(pT*(pT*pT)))
		SF_low = ((0.911657+(0.00142008*pT))+(-2.87569e-06*(pT*pT)))+(1.76619e-09*(pT*(pT*pT)))
		SF_hi  = ((0.964787+(0.00219574*pT))+(-4.85552e-06*(pT*pT)))+(3.09457e-09*(pT*(pT*pT)))
	elif eta > 1.5 and eta < 2.4 :
		SF = ((1.00022+(0.0010998*pT))+(-3.10672e-06*(pT*pT)))+(2.35006e-09*(pT*(pT*pT)))
		SF_low = ((0.970045+(0.000862284*pT))+(-2.31714e-06*(pT*pT)))+(1.68866e-09*(pT*(pT*pT)))
		SF_hi  = ((1.03039+(0.0013358*pT))+(-3.89284e-06*(pT*pT)))+(3.01155e-09*(pT*(pT*pT)))
	else :
		print 'WARNING: EVENT WITH JET OUTSIDE OF ETA RANGE! CANNOT GET BTAG EFFICIENCY! |eta| = '+str(eta)+''
	return SF,SF_low,SF_hi
