#lepHelper: contains helper functions for lepton selection
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

import ROOT
from math import *

#Global variables
#cutflow
CUTFLOW_EXACTLY_ONE_LEPTON = 2
CUTFLOW_OTHER_LEPTON_VETO = 3
CUTFLOW_JET_PRESELECTION = 4
CUTFLOW_LEPTON_2D = 5
CUTFLOW_LEPTON_TRIANGLE = 6
#cut values
MU_PT_MIN = 45 #GeV
MU_ETA_MAX = 2.1
EL_PT_MIN = 35 #GeV
EL_ETA_MAX = 2.5
OTHER_MU_PT_MIN = 35 #GeV
OTHER_MU_ETA_MAX = 2.1
OTHER_EL_PT_MIN = 35 #GeV
OTHER_EL_ETA_MAX = 2.5
ONE_LARGE_JET_PT_MIN = 150. #GeV
ONE_LARGE_ONE_SMALL_JET_PT_MIN = 50. #GeV
JETS_ETA_MAX = 2.4
JET_PT_MIN = 25. #GeV
DR_MIN = 0.5
REL_PT_MIN = 25. #GeV

#leptonCuts
#takes in lepton type, sideband switch, lists of muon and electron variables
#returns: 
#	1) Index of single valid lepton OR
#	2) Negative value of cutflow failure point
def leptonCuts(lep_type,sideband,muVars,elVars,metVars,jetVars,jetVars_large,control_plots) :
	#make list of lepton candidate and "other" lepton indices
	lep_cands = []
	other_leps = []
#	print 'new event ----------------------------------------' #DEBUGGING
	if lep_type == 0 : #muons
#		print '# of gen muons = '+str(len(muVars[11]))+'' #DEBUGGING
		for i in range(len(muVars[0])) :
#			s = '	istight cut '   #DEBUGGING
#			if muVars[9][i] == 1 : #DEBUGGING
#				s+='passed, ' #DEBUGGING
#			else :  #DEBUGGING
#				s+= 'FAILED ('+str(muVars[9][i])+'), ' #DEBUGGING
#			s+= ' pt cut '   #DEBUGGING
#			if muVars[0][i] > MU_PT_MIN : #DEBUGGING
#				s+='passed, ' #DEBUGGING
#			else : #DEBUGGING
#				s+= 'FAILED, ('+str(muVars[0][i])+'), ' #DEBUGGING
#			s+= ' eta cut '  #DEBUGGING
#			if fabs(muVars[1][i]) < MU_ETA_MAX : #DEBUGGING
#				s+='passed' #DEBUGGING
#			else : #DEBUGGING
#				s+= 'FAILED ('+str(muVars[1][i])+')' #DEBUGGING
#			print s  #DEBUGGING
			if i<2 and muVars[9][i] == 1 :
				control_plots[2*i].Fill(muVars[0][i]); control_plots[2*i+1].Fill(muVars[1][i])
			if muVars[0][i] > MU_PT_MIN and fabs(muVars[1][i]) < MU_ETA_MAX and muVars[9][i] == 1 :
				lep_cands.append(i)
		for i in range(len(elVars[0])) :
			if i<2 and len(lep_cands)==1 and elVars[7][i] == 1 :
				control_plots[4+2*i].Fill(elVars[0][i]); control_plots[5+2*i].Fill(elVars[1][i])
			if elVars[0][i] > OTHER_EL_PT_MIN and fabs(elVars[1][i]) < OTHER_EL_ETA_MAX and elVars[7][i] ==1 :
				other_leps.append(i)
	elif lep_type == 1 : #electrons
		for i in range(len(elVars[0])) :
			if i<2 and elVars[6][i] == 1 :
				control_plots[2*i].Fill(elVars[0][i]); control_plots[2*i+1].Fill(elVars[1][i])
			if elVars[0][i] > EL_PT_MIN and fabs(elVars[1][i]) < EL_ETA_MAX and elVars[6][i] ==1 :
				lep_cands.append(i)
		for i in range(len(muVars[0])) :
			if i<2 and len(lep_cands)==1 and muVars[10][i] == 1 :
				control_plots[4+2*i].Fill(muVars[0][i]); control_plots[5+2*i].Fill(muVars[1][i])
			if muVars[0][i] > OTHER_MU_PT_MIN and fabs(muVars[1][i]) < OTHER_MU_ETA_MAX and muVars[10][i] == 1 :
				other_leps.append(i)
	#check for exactly one lepton and no "other" leptons
	if len(lep_cands) != 1 :
#		print '# of lepton candidates = '+str(len(lep_cands))+'' #DEBUGGING
		return -1*CUTFLOW_EXACTLY_ONE_LEPTON
	if len(other_leps) != 0 :
#		print ' # of other leptons = '+str(len(other_leps))+'' #DEBUGGING'
		return -1*CUTFLOW_OTHER_LEPTON_VETO
	lepton_index = lep_cands[0]

	#check that there are at least one large and one small jet
	if (len(jetVars[0])<1 or len(jetVars_large[0])<1) :
		return -1*CUTFLOW_JET_PRESELECTION
	control_plots[8].Fill(jetVars_large[0][0],jetVars[0][0])
	#cut that there's at least one very hard large jet, or at least one hard large jet and one hard small jet
	#in the eta range
	isGood = False
	if jetVars_large[0][0] > ONE_LARGE_JET_PT_MIN and fabs(jetVars_large[1][0]) < JETS_ETA_MAX :
		isGood = True
	elif ( jetVars_large[0][0] > ONE_LARGE_ONE_SMALL_JET_PT_MIN and jetVars[0][0] > ONE_LARGE_ONE_SMALL_JET_PT_MIN
		and fabs(jetVars_large[1][0]) < JETS_ETA_MAX and fabs(jetVars[1][0]) < JETS_ETA_MAX ) :
		isGood = True
	if not isGood :
		return -1*CUTFLOW_JET_PRESELECTION

	#find the nearest jet with pT > some minimum 
	nearest_jet_dR = 1000000.
	nearest_jet_vec = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
	lepvec = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
	if lep_type == 0 :
		lepvec.SetPtEtaPhiM(muVars[0][lepton_index],muVars[1][lepton_index],
			muVars[2][lepton_index],muVars[3][lepton_index])
	elif lep_type == 1 :
		lepvec.SetPtEtaPhiM(elVars[0][lepton_index],elVars[1][lepton_index],
			elVars[2][lepton_index],elVars[3][lepton_index])
#	print 'lepvec = ('+str(lepvec.Px())+','+str(lepvec.Py())+','+str(lepvec.Pz())+','+str(lepvec.E())+')' #DEBUGGING
	for i in range(len(jetVars[0])) :
		if jetVars[0][i] < JET_PT_MIN :
			continue
		thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		thisJet.SetPtEtaPhiM(jetVars[0][i],jetVars[1][i],jetVars[2][i],jetVars[3][i])
		thisJetdR = lepvec.DeltaR(thisJet)
		if thisJetdR < nearest_jet_dR :#and thisJetdR > 0.1 :
			nearest_jet_vec.SetPtEtaPhiM(thisJet.Pt(),thisJet.Eta(),thisJet.Phi(),thisJet.M())
			nearest_jet_dR = thisJetdR
#	print 'nearest_jet_dR = '+str(nearest_jet_dR)+'' #DEBUGGING
	#make lepton 2D cut
	control_plots[9].Fill(nearest_jet_dR,lepvec.Pt(nearest_jet_vec.Vect()))
	if nearest_jet_dR < DR_MIN and lepvec.Pt(nearest_jet_vec.Vect()) < REL_PT_MIN : #yes, this SHOULD say "and".
		return -1*CUTFLOW_LEPTON_2D

	#make electron triangle cut
	if lep_type == 1 :
		metVec = ROOT.TLorentzVector() 
		metVec.SetPtEtaPhiM(metVars[0][0],0.0,metVars[1][0],0.0)
		metE = metVec.E()
		firstJet = ROOT.TLorentzVector() 
		firstJet.SetPtEtaPhiM(jetVars[0][0],jetVars[1][0],jetVars[2][0],jetVars[3][0])
		el_val = abs(lepvec.DeltaPhi(metVec)-1.5); jet_val = abs(firstJet.DeltaPhi(metVec)-1.5)
		control_plots[10].Fill(metE,el_val); control_plots[11].Fill(metE,jet_val)
		cut_val = (1.5/75.)*metE
		if not (el_val < cut_val and jet_val < cut_val) :
			return -1*CUTFLOW_LEPTON_TRIANGLE

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
	return (lepvec,int(lepcharge))
	#return (ROOT.TLorentzVector(1.0,0.0,0.0,1.0), 1) #DEBUG RETURN

#leptonCleaningKludge implements not completely correct lepton cleaning for AK4 and AK8 jets
#takes in the lepton and jet fourvectors
#subtracts off the fourvectors of loose, non-isolated leptons from jets within 0.5
#returns the subtracted jet variable lists
def leptonCleaningKludge(muVars,elVars,jetVars_small,jetVars_large) :
	newJetVars_AK4 = []; newJetVars_AK8 = []
	for i in range(len(jetVars_small)) :
		newJetVars_AK4.append([])
		for j in range(len(jetVars_small[i])) :
			newJetVars_AK4[i].append(jetVars_small[i][j])
	for i in range(len(jetVars_large)) :
		newJetVars_AK8.append([])
		for j in range(len(jetVars_large[i])) :
			newJetVars_AK8[i].append(jetVars_large[i][j])
#	print 'old newJetVars_AK4 = '+str(newJetVars_AK4)+'' #DEBUGGING
	#look at all the leptons and find the non-isolated ones
	for i in range(len(muVars[0])) :
		#only look at the loose leptons
		if muVars[7][i] != 1 :
			continue
		#calculate the muon isolation
		muIso = (muVars[5][i] + max(0.,muVars[6][i]+muVars[7][i]-0.5*muVars[8][i]))/muVars[0][i]
		#for non-isolated muons
		if muIso < 0.12 :
			continue
		thisMuon = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		thisMuon.SetPtEtaPhiM(muVars[0][i],muVars[1][i],muVars[2][i],muVars[3][i])
		#find AK4 jets within 0.5 of the lepton and subtract off the lepton pT
		for j in range(len(newJetVars_AK4[0])) :
			thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0) 
			thisJet.SetPtEtaPhiM(jetVars_small[0][j],jetVars_small[1][j],jetVars_small[2][j],jetVars_small[3][j])
			if thisJet.DeltaR(thisMuon) < 0.4 and thisMuon.Pt()<thisJet.Pt() :
#				print 'removing a muon from an AK4 jet' #DEBUGGING
				thisNewJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0) 
				thisNewJet.SetPtEtaPhiM(newJetVars_AK4[0][j],newJetVars_AK4[1][j],
					newJetVars_AK4[2][j],newJetVars_AK4[3][j])
				newJet = thisNewJet - thisMuon
				newJetVars_AK4[0][j], newJetVars_AK4[1][j] = newJet.Pt(), newJet.Eta() 
				newJetVars_AK4[2][j], newJetVars_AK4[3][j] = newJet.Phi(), newJet.M()
		#find AK8 jets within 0.5 of the lepton and subtract off the lepton pT
		for j in range(len(newJetVars_AK8[0])) :
			thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0) 
			thisJet.SetPtEtaPhiM(jetVars_large[0][j],jetVars_large[1][j],jetVars_large[2][j],jetVars_large[3][j])
			if thisJet.DeltaR(thisMuon) < 0.8 and thisMuon.Pt()<thisJet.Pt() :
				thisNewJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0) 
				thisNewJet.SetPtEtaPhiM(newJetVars_AK8[0][j],newJetVars_AK8[1][j],
					newJetVars_AK8[2][j],newJetVars_AK8[3][j])
				newJet = thisNewJet - thisMuon
				newJetVars_AK8[0][j], newJetVars_AK8[1][j] = newJet.Pt(), newJet.Eta() 
				newJetVars_AK8[2][j], newJetVars_AK8[3][j] = newJet.Phi(), newJet.M()
	for i in range(len(elVars[0])) :
		#only look at the loose leptons
		#if elVars[7][i] != 1 and elVars[6][i] != 1:
		#	continue
		#for non-isolated electrons
		#if elVars[5][i] < 0.10 :
		#	continue
		thisElectron = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		thisElectron.SetPtEtaPhiM(elVars[0][i],elVars[1][i],elVars[2][i],elVars[3][i])
		#find AK4 jets within 0.5 of the lepton and subtract off the lepton pT
		for j in range(len(newJetVars_AK4[0])) :
			thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0) 
			thisJet.SetPtEtaPhiM(jetVars_small[0][j],jetVars_small[1][j],jetVars_small[2][j],jetVars_small[3][j])
			if thisJet.DeltaR(thisElectron) < 0.4 and thisElectron.Pt()<thisJet.Pt() :
#					print 'removed an electron from an AK4 jet' #DEBUGGING
				thisNewJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0) 
				thisNewJet.SetPtEtaPhiM(newJetVars_AK4[0][j],newJetVars_AK4[1][j],
					newJetVars_AK4[2][j],newJetVars_AK4[3][j])
				newJet = thisNewJet - thisElectron
				newJetVars_AK4[0][j], newJetVars_AK4[1][j] = newJet.Pt(), newJet.Eta() 
				newJetVars_AK4[2][j], newJetVars_AK4[3][j] = newJet.Phi(), newJet.M()
		#find AK8 jets within 0.5 of the lepton and subtract off the lepton pT
		for j in range(len(newJetVars_AK8[0])) :
			thisJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0) 
			thisJet.SetPtEtaPhiM(jetVars_large[0][j],jetVars_large[1][j],jetVars_large[2][j],jetVars_large[3][j])
			if thisJet.DeltaR(thisElectron) < 0.8 and thisElectron.Pt()<thisJet.Pt() :
				thisNewJet = ROOT.TLorentzVector(1.0,0.0,0.0,1.0) 
				thisNewJet.SetPtEtaPhiM(newJetVars_AK8[0][j],newJetVars_AK8[1][j],
					newJetVars_AK8[2][j],newJetVars_AK8[3][j])
				newJet = thisNewJet - thisElectron
				newJetVars_AK8[0][j], newJetVars_AK8[1][j] = newJet.Pt(), newJet.Eta() 
				newJetVars_AK8[2][j], newJetVars_AK8[3][j] = newJet.Phi(), newJet.M()
#	print 'new newJetVars_AK4 = '+str(newJetVars_AK4)+'' #DEBUGGING
	return (newJetVars_AK4,newJetVars_AK8)
