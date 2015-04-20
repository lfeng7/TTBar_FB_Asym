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
MAX_LEP_TOP_MASS_WINDOW = 210. #GeV
MIN_HAD_TOP_PT = 350. #GeV
MIN_HAD_TOP_MASS = 130. #GeV
MAX_HAD_TOP_TAU32 = 0.7
MIN_HAD_W_PT = 200. #GeV
MIN_HAD_W_MASS_WINDOW = 65. #GeV
MAX_HAD_W_MASS_WINDOW = 105. #GeV
MAX_HAD_W_TAU21 = 0.75
MIN_HAD_B_PT = 25. #GeV
MIN_HAD_TOP_MASS_WINDOW = 140. #GeV
MAX_HAD_TOP_MASS_WINDOW = 210. #GeV
MIN_HAD_W_B_DELTAR = 0.8

#selectJets takes in the top type (1 or 2), the lepton vector, the met vectors, and lists of AK4 and AK8 jet variables
#returns:
#	1) a 4 or 5 element list of the following
#		-leptonic b (fourvector, CSV value, flavour), hadronic top (fourvector, 1.0), 
#			3 hadronic subjets (fourvector, 1.0) for type 1 or
#		-leptonic b (fourvector, CSV value, flavour), hadronic W (fourvector,tau1, tau2, tau3),
#			hadronic b (fourvector,1.0) (needed to maintain a parallel structure here), for type 2 THEN
#		-the btagging efficiency scale factor, low error, and hi error for MC events OR
#	2) a 1-element list of negative value of a cutflow failure point
def selectJets(isdata,lepvec,met1vec,met2vec,jetvars_small,jetvars_large,jetvars_tt,jet_control_plots) :
	toptype = 0
	return_list = []
	#first find the leptonic b by looking for loosely b-tagged jets in the lepton hemisphere
	lepbcand = ROOT.TLorentzVector(1.0,0.0,0.0,1.0); lepbcand_CSV = 0.5; lepbcand_flavour = -100
	best_comb_mass_offset = 100000.
#	print 'NEW EVENT LEPTONIC B SELECTION: ' #DEBUGGING
	for i in range(len(jetvars_small[0])) :
		#min pT cut
		jet_control_plots[0].Fill(jetvars_small[0][i])
		if jetvars_small[0][i] < MIN_LEPB_PT :
			continue
		thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		thisJet.SetPtEtaPhiM(jetvars_small[0][i],jetvars_small[1][i],jetvars_small[2][i],jetvars_small[3][i])
		#Fill a control plot with the jet deltaphi from the lepton
		lep_dPhi = abs(lepvec.DeltaPhi(thisJet))
		jet_control_plots[1].Fill(lep_dPhi)
		#see if this jet is the best one yet in terms of combined mass
		comb_mass_1 = (lepvec+met1vec+thisJet).M(); comb_mass_2 = (lepvec+met2vec+thisJet).M()
		comb_mass_offset = min(abs(comb_mass_1-TOP_MASS),abs(comb_mass_2-TOP_MASS))
		if comb_mass_offset<best_comb_mass_offset :
			best_comb_mass_offset = comb_mass_offset
#			print '	best offset = '+str(best_comb_mass_offset)+'' #DEBUGGING
			lepbcand.SetPtEtaPhiM(thisJet.Pt(),thisJet.Eta(),thisJet.Phi(),thisJet.M())
			lepbcand_CSV = jetvars_small[4][i]
			if isdata == 0 :
				lepbcand_flavour = jetvars_small[5][i]
			else :
				lepbcand_flavour = 0
	if best_comb_mass_offset == 100000. : #didn't find any jets with the minimum pT
		return ([-1*CUTFLOW_EXACTLY_ONE_LEPTONIC_BCAND],1.0,1.0,1.0,toptype)
	#Check that this best combined mass is within the mass window
	leptopvec = ROOT.TLorentzVector()
	comb_mass_1 = (lepvec+met1vec+lepbcand).M(); comb_mass_2 = (lepvec+met2vec+lepbcand).M()
	comb_mass_offset = min(abs(comb_mass_1-TOP_MASS),abs(comb_mass_2-TOP_MASS))
	if comb_mass_offset == abs(comb_mass_1-TOP_MASS) :
		jet_control_plots[2].Fill(comb_mass_1)
		if comb_mass_1 < MIN_LEP_TOP_MASS_WINDOW or comb_mass_1 > MAX_LEP_TOP_MASS_WINDOW :
			return ([-1*CUTFLOW_EXACTLY_ONE_LEPTONIC_BCAND],1.0,1.0,1.0,toptype)
		leptopvec = lepvec+met1vec+lepbcand
	else :
		jet_control_plots[2].Fill(comb_mass_2) 
		if comb_mass_2 < MIN_LEP_TOP_MASS_WINDOW or comb_mass_2 > MAX_LEP_TOP_MASS_WINDOW :
			return ([-1*CUTFLOW_EXACTLY_ONE_LEPTONIC_BCAND],1.0,1.0,1.0,toptype)
		leptopvec = lepvec+met2vec+lepbcand
	#cut that the candidate jet is loosely b-tagged and calculate the scalefactor and error
	jet_control_plots[3].Fill(lepbcand_CSV)
	if lepbcand_CSV < CSVL_WORKING_POINT :
		return ([-1*CUTFLOW_EXACTLY_ONE_LEPTONIC_BCAND],1.0,1.0,1.0,toptype)
	btag_weight = 0.0; btag_weight_low = 0.0; btag_weight_hi = 0.0
	if isdata == 0 :
		btag_weight, btag_weight_low, btag_weight_hi = getSF((lepbcand,lepbcand_CSV,lepbcand_flavour),CSVL_WORKING_POINT)
	if btag_weight == 0.0 :
		btag_weight = 1.0; btag_weight_low = 1.0; btag_weight_hi = 1.0
	#add the b candidate to the return list
	return_list.append((lepbcand,lepbcand_CSV,lepbcand_flavour))
	
	#now get the one or two jets making up the hadronic top 
	#see if the even has any fully merged top jet candidates
	hascand = False
	for i in range(len(jetvars_tt[0])) :
		jetpt = jetvars_tt[0][i]
		nsub = jetvars_tt[27][i]
		if jetpt > MIN_HAD_TOP_PT and nsub == 3 :
			hascand = True
			break
	if not hascand :
		toptype = 2 #if there are no fully merged jet candidates, then we have to reconstruct the partially merged top
#	print 'toptype = '+str(toptype)+'' #DEBUGGING
	if toptype == 2 :
		twojetstuple = selectJetsType2Tops(leptopvec,jetvars_large,jet_control_plots)
		for jettuple in twojetstuple :
			return_list.append(jettuple)
	else :
		toptype = 1
		fourjetstuple = selectJetsType1Tops(leptopvec,jetvars_tt,jet_control_plots)
		for jettuple in fourjetstuple :
			return_list.append(jettuple)

	#check if the hadronic top reconstruction failed and return the cutflow value
#	print 'return_list = '+str(return_list)+'' #DEBUGGING
	if len(return_list[1]) == 1 :
		return (return_list[1],1.0,1.0,1.0,toptype)

	#otherwise return the tuple of the jet list and the btag weights and the top type
	return (return_list, btag_weight, btag_weight_low, btag_weight_hi, toptype)

#selectJetsType1Tops does the selection for type 1 (fully merged) tops
def selectJetsType1Tops(leptopvec,jetvars_tt,jet_control_plots) :
	#look through all the top tagger jets
	hadtopcands = []
	for i in range(len(jetvars_tt[0])) :
		#hard cut on jet pT
		jet_control_plots[4].Fill(jetvars_tt[0][i]) 
		if jetvars_tt[0][i] < MIN_HAD_TOP_PT :
			continue
		#hard cut on nsubjets
		jet_control_plots[5].Fill(jetvars_tt[27][i]) 
		if jetvars_tt[27][i] != 3 :
			continue
		#mass window cut near the W mass for mmin
		mmin = jetvars_tt[25][i]
		jet_control_plots[6].Fill(mmin)
#		if mmin < MIN_HAD_W_MASS_WINDOW or mmin > MAX_HAD_W_MASS_WINDOW :
#			continue
		#mass window cut near the top mass for mjet
		mjet = jetvars_tt[26][i]
		jet_control_plots[7].Fill(mjet)
#		if mjet < MIN_HAD_TOP_MASS_WINDOW or mjet > MAX_HAD_TOP_MASS_WINDOW :
#			continue
		#add this jet index to the list of candidates
		hadtopcands.append(i)
#	print 'len(hadtopcands) = '+str(len(hadtopcands))+'' #DEBUGGING
	#require exactly one hadronic top candidate
	if len(hadtopcands) != 1 :
		return ([-1*CUTFLOW_EXACTLY_ONE_HADRONIC_TOPCAND],[-1*CUTFLOW_EXACTLY_ONE_HADRONIC_TOPCAND])
	#return the hadronic top candidate with its subjets
	index = hadtopcands[0]
	topjet = ROOT.TLorentzVector(); sub0 = ROOT.TLorentzVector(); sub1 = ROOT.TLorentzVector(); sub2 = ROOT.TLorentzVector()
	topjet.SetPtEtaPhiM(jetvars_tt[0][index],jetvars_tt[1][index],jetvars_tt[2][index],jetvars_tt[3][index])
	sub0.SetPtEtaPhiM(jetvars_tt[4][index],jetvars_tt[5][index],jetvars_tt[6][index],jetvars_tt[7][index])
	sub1.SetPtEtaPhiM(jetvars_tt[8][index],jetvars_tt[9][index],jetvars_tt[10][index],jetvars_tt[11][index])
	sub2.SetPtEtaPhiM(jetvars_tt[12][index],jetvars_tt[13][index],jetvars_tt[14][index],jetvars_tt[15][index])
	return ((topjet,1.0),(sub0,1.0),(sub1,1.0),(sub2,1.0))
	#return (ROOT.TLorentzVector(1.0,0.0,0.0,1.0),0.5,0.5,0.5,0.5,0) #DEBUG RETURN

#selectJetsType2Tops does the selection for type 2 (partially merged) tops
def  selectJetsType2Tops(leptopvec,jetvars_large,jet_control_plots) :
	#look through the large jets for W candidates
	hadWcands = []
	for i in range(len(jetvars_large[0])) :
		#hard cut on pT
		jet_control_plots[8].Fill(jetvars_large[0][i]) 
		if jetvars_large[0][i] < MIN_HAD_W_PT :
			continue
		#check that it's in the W mass window
		jetMass = jetvars_large[3][i]
		jet_control_plots[9].Fill(jetMass)
		if jetMass < MIN_HAD_W_MASS_WINDOW or jetMass > MAX_HAD_W_MASS_WINDOW :
			continue
		#check that it's on the hadronic side
		thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		thisJet.SetPtEtaPhiM(jetvars_large[0][i],jetvars_large[1][i],jetvars_large[2][i],jetvars_large[3][i])
		jet_control_plots[10].Fill(abs(leptopvec.DeltaPhi(thisJet)))
		if abs(leptopvec.DeltaPhi(thisJet)) < pi/2. :
			continue
		#get the matched unpruned jet's tau variables and cut on tau2/tau1
		matchedJet_taus = get_taus(thisJet,jetvars_large)
		tau21 = matchedJet_taus[1]/matchedJet_taus[0]
		jet_control_plots[11].Fill(tau21)
		if tau21 > MAX_HAD_W_TAU21 :
			continue
		#append to the list of W candidates
		hadWcands.append((thisJet,matchedJet_taus[0],matchedJet_taus[1],matchedJet_taus[2]))
	#require exactly one hadronic W candidate
	jet_control_plots[12].Fill(len(hadWcands))
	if len(hadWcands) != 1 :
		return ([-1*CUTFLOW_EXACTLY_ONE_HADRONIC_WCAND],[-1*CUTFLOW_EXACTLY_ONE_HADRONIC_WCAND])

	#now look through the large jets for hadronic b candidates
	hadbcand = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
	best_comb_mass_offset = 100000.
	for i in range(len(jetvars_large[0])) :
		#hard cut on pT 
		if jetvars_large[0][i] < MIN_HAD_B_PT :
			continue
		thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		thisJet.SetPtEtaPhiM(jetvars_large[0][i],jetvars_large[1][i],jetvars_large[2][i],jetvars_large[3][i])
		#check that it's in the hadronic hemisphere
		jet_control_plots[13].Fill(abs(leptopvec.DeltaPhi(thisJet)))
		if abs(leptopvec.DeltaPhi(thisJet)) < pi/2. :
			continue
		#check that it combines with the W to give a mass in the top mass window
		combMass = (thisJet+hadWcands[0][0]).M()
		jet_control_plots[14].Fill(combMass)
		if combMass < MIN_HAD_TOP_MASS_WINDOW or combMass > MAX_HAD_TOP_MASS_WINDOW :
			continue
		#check that it's sufficiently separated from the W jet
		jet_control_plots[15].Fill(thisJet.DeltaR(hadWcands[0][0]))
		if thisJet.DeltaR(hadWcands[0][0]) < MIN_HAD_W_B_DELTAR :
			continue
		#if it's giving the best combined mass yet, it's the new hadronic b jet
		comb_mass_offset = abs(combMass-TOP_MASS)
		if comb_mass_offset < best_comb_mass_offset :
			best_comb_mass_offset = comb_mass_offset
			hadbcand.SetPtEtaPhiM(thisJet.Pt(),thisJet.Eta(),thisJet.Phi(),thisJet.M())
	#require exactly one hadronic b candidate
	if best_comb_mass_offset == 100000. :
		return ([-1*CUTFLOW_EXACTLY_ONE_HADRONIC_BCAND],[-1*CUTFLOW_EXACTLY_ONE_HADRONIC_BCAND])
	jet_control_plots[16].Fill((hadbcand+hadWcands[0][0]).M()) 		

	return (hadWcands[0],(hadbcand,1.0))
	#return [(ROOT.TLorentzVector(1.0,0.0,0.0,1.0),0.5,0),
	#(ROOT.TLorentzVector(1.0,0.0,0.0,1.0),0.5,0.5,0.5,0.5,0)] #DEBUG RETURN

#get_taus takes a pruned jet vector and the list of all the other jet variables. It returns a 3-tuple 
#of the "tau" substructure variables for the unpruned jet that is the best match to the given pruned jet
def get_taus(jetVec,jetVars) :
	#look through all the Unpruned CA8 Jets and try to match to the given jet
	best_dR = 10000000.
	bestJetIndex = -1
	for i in range(len(jetVars[4])) :
		newJet = ROOT.TLorentzVector()
		newJet.SetPtEtaPhiM(jetVars[4][i],jetVars[5][i],jetVars[6][i],jetVars[7][i])
		if jetVec.DeltaR(newJet) < best_dR :
			bestJetIndex = i
			best_dR = jetVec.DeltaR(newJet)
	return (jetVars[8][bestJetIndex],jetVars[9][bestJetIndex],jetVars[10][bestJetIndex])

##################################################################################################
##########						b-tagging efficiency stuff 								##########
##################################################################################################

#btagging efficiency stuff
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation53XReReco
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods

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
	if eta < 0.5 :
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
