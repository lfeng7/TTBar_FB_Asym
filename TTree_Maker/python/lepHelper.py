#lepHelper: contains helper functions for lepton selection
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

import ROOT
from math import *

#Global variables
#cutflow
CUTFLOW_EXACTLY_ONE_LEPTON = 2
CUTFLOW_LEPTON_2D = 3
#cut values
MU_PT_MIN = 45 #GeV
MU_ETA_MAX = 2.1
EL_PT_MIN = 35 #GeV
EL_ETA_MAX = 2.5
OTHER_MU_PT_MIN = 35 #GeV
OTHER_MU_ETA_MAX = 2.1
OTHER_EL_PT_MIN = 35 #GeV
OTHER_EL_ETA_MAX = 2.5
JET_PT_MIN = 25 #GeV
DR_MIN = 0.5
REL_PT_MIN = 25 #GeV

#leptonCuts
#takes in lepton type, sideband switch, lists of muon and electron variables
#returns: 
#	1) Index of single valid lepton OR
#	2) Negative value of cutflow failure point
def leptonCuts(lep_type,sideband,muVars,elVars,jetVars) :
	#make list of lepton candidate and "other" lepton indices
	lep_cands = []
	other_leps = []
	if lep_type == 0 : #muons
		for i in range(len(muVars[0])) :
			if muVars[0][i] > MU_PT_MIN and fabs(muVars[1][i]) < MU_ETA_MAX and muVars[9][i] == 1 :
				lep_cands.append(i)
		for i in range(len(elVars[0])) :
			if elVars[0][i] > OTHER_EL_PT_MIN and fabs(elVars[1][i]) < OTHER_EL_ETA_MAX :
				other_leps.append(i)
	elif lep_type == 1 : #electrons
		for i in range(len(elVars[0])) :
			if elVars[0][i] > EL_PT_MIN and fabs(elVars[1][i]) < EL_ETA_MAX :
				lep_cands.append(i)
		for i in range(len(muVars[0])) :
			if muVars[0][i] > OTHER_MU_PT_MIN and fabs(muVars[1][i]) < OTHER_MU_ETA_MAX :
				other_leps.append(i)
	#check for exactly one lepton and no "other" leptons
	if len(lep_cands) != 1 or len(other_leps) != 0 :
		return -1*CUTFLOW_EXACTLY_ONE_LEPTON
	lepton_index = lep_cands[0]

	#find the nearest jet with minimum pT
	nearest_jet_dR = 10000.
	nearest_jet_vec = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
	lepvec = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
	if lep_type == 0 :
		lepvec.SetPtEtaPhiM(muVars[0][lepton_index],muVars[1][lepton_index],muVars[2][lepton_index],muVars[3][lepton_index])
	elif lep_type == 1 :
		lepvec.SetPtEtaPhiM(elVars[0][lepton_index],elVars[1][lepton_index],elVars[2][lepton_index],elVars[3][lepton_index])
	for i in range(len(jetVars[0])) :
		if jetVars[0][i] < JET_PT_MIN :
			continue
		thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		thisJet.SetPtEtaPhiM(jetVars[0][i],jetVars[1][i],jetVars[2][i],jetVars[3][i])
		thisJetDeltaR = lepvec.DeltaR(thisJet)
		if thisJetDeltaR < nearest_jet_dR :
			nearest_jet_vec.SetPtEtaPhiM(thisJet.Pt(),thisJet.Eta(),thisJet.Phi(),thisJet.M())
			nearest_jet_dR = thisJetDeltaR
	#make lepton 2D cut
	if nearest_jet_dR < DR_MIN and lepvec.Pt(nearest_jet_vec.Vect()) < REL_PT_MIN :
		return -1*CUTFLOW_LEPTON_2D

	#return the index of the single valid lepton
	return lepton_index
	#return 0 #DEBUG RETURN

#getLeptonFourVec
#takes in lepton type, lists of muon/electron variables, index of single valid lepton
#returns a tuple of (TLorentzVector of single valid lepton, lepton charge)
def getLeptonFourVec(lep_type,muVars,elVars,lep_index) :
	lepvec = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
	lepcharge = 0
	if lep_type == 0 : #muons
		lepvec.SetPtEtaPhiM(muVars[0][lep_index],muVars[1][lep_index],muVars[2][lep_index],muVars[3][lep_index])
		lepcharge = muVars[4][lep_index]
	elif lep_type == 1 : #electrons
		lepvec.SetPtEtaPhiM(elVars[0][lep_index],elVars[1][lep_index],elVars[2][lep_index],elVars[3][lep_index])
		lepcharge = elVars[4][lep_index]
	return (lepvec,lepcharge)
	#return (ROOT.TLorentzVector(1.0,0.0,0.0,1.0), 1) #DEBUG RETURN