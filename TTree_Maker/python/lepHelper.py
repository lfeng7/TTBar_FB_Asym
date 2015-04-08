#lepHelper: contains helper functions for lepton selection
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

import ROOT
from math import *

#Global variables
#cutflow
CUTFLOW_JET_PRESELECTION = 2
CUTFLOW_LEPTON_2D = 3
CUTFLOW_LEPTON_TRIANGLE = 4
CUTFLOW_EXACTLY_ONE_LEPTON = 5
CUTFLOW_OTHER_LEPTON_VETO = 6
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
	#make list of lepton candidate and "other" lepton indices after quick kinematic cuts
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

	#do jet preselection cut
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
#		print '	failed jet preselection' #DEBUGGING
		return -1*CUTFLOW_JET_PRESELECTION
	
	#do the lepton-jet cleaning and make "isolation" cuts
	iso_leps = []
	iso_other_leps = []
	#for each lepton, find the nearest jet with pT > some minimum 
	metvec = ROOT.TLorentzVector()
	metvec.SetPtEtaPhiM(metVars[0][0],0.0,metVars[1][0],0.0)
	for i in range(len(lep_cands)) :
		nearest_jet_dR = 1000000.
		nearest_jet_vec = ROOT.TLorentzVector()
		lepvec = ROOT.TLorentzVector()
		if lep_type == 0 :
			lepvec.SetPtEtaPhiM(muVars[0][i],muVars[1][i],muVars[2][i],muVars[3][i])
		elif lep_type == 1 :
			lepvec.SetPtEtaPhiM(elVars[0][i],elVars[1][i],elVars[2][i],elVars[3][i])
		#check the isolation
		iso_func_return = lepIsIso(lep_type,sideband,lepvec,metvec,jetVars,control_plots)
#		print 'iso_func_return = '+str(iso_func_return)+'' #DEBUGGING
		if iso_func_return == 0 :
			iso_leps.append(lep_cands[i])
	for i in range(len(other_leps)) :
		nearest_jet_dR = 1000000.
		nearest_jet_vec = ROOT.TLorentzVector()
		lepvec = ROOT.TLorentzVector()
		if lep_type == 1 :
			lepvec.SetPtEtaPhiM(muVars[0][i],muVars[1][i],muVars[2][i],muVars[3][i])
			iso_func_return = lepIsIso(0,sideband,lepvec,metvec,jetVars,control_plots)
		elif lep_type == 0 :
			lepvec.SetPtEtaPhiM(elVars[0][i],elVars[1][i],elVars[2][i],elVars[3][i])
			iso_func_return = lepIsIso(1,sideband,lepvec,metvec,jetVars,control_plots)
		if iso_func_return == 0 :
			iso_other_leps.append(other_leps[i])

	print '      # leptons = %d, # iso leptons = %d\n      # other leptons = %d, # iso other leptons = %d'%(len(lep_cands),len(iso_leps),len(other_leps),len(iso_other_leps)) #DEBUGGING

	#check for exactly one lepton that passed the isolation cuts
	if len(iso_leps) != 1 :
#		print '# of lepton candidates = '+str(len(iso_leps))+'' #DEBUGGING
		return -1*CUTFLOW_EXACTLY_ONE_LEPTON
	#check for no "other" leptons passing the isolation cut
	if len(iso_other_leps) != 0 :
#		print ' # of other leptons = '+str(len(iso_other_leps))+'' #DEBUGGING'
		return -1*CUTFLOW_OTHER_LEPTON_VETO

	#return the index of the single valid lepton
	return iso_leps[0]
	#return 0 #DEBUG RETURN

#getLeptonFourVec
#takes in lepton type, lists of muon/electron variables, index of single valid lepton
#returns a tuple of (TLorentzVector of single valid lepton, lepton charge)
def getLeptonFourVec(lep_type,muVars,elVars,lep_index) :
	lepvec = ROOT.TLorentzVector()
	lepcharge = 0
	if lep_type == 0 : #muons
		lepvec.SetPtEtaPhiM(muVars[0][lep_index],muVars[1][lep_index],muVars[2][lep_index],muVars[3][lep_index])
		lepcharge = muVars[4][lep_index]
	elif lep_type == 1 : #electrons
		lepvec.SetPtEtaPhiM(elVars[0][lep_index],elVars[1][lep_index],elVars[2][lep_index],elVars[3][lep_index])
		lepcharge = elVars[4][lep_index]
	return (lepvec,int(lepcharge))
	#return (ROOT.TLorentzVector(), 1) #DEBUG RETURN

#lepIsIso takes in lepton type, whether or not it's a sideband file, a lepton fourvector, the MET fourvector,
#the list of jet variables, and the list of control plots
#it makes the 2D cut (and the triangle cut in the case of electrons) and returns
#0 if the lepton passes
#a cutflow if it fails
def lepIsIso(lep_type,sideband,lepvec,metvec,jetVars,control_plots) :
	#find the nearest jet with pT > some minimum 
	nearest_jet_dR = 1000000.
	nearest_jet_vec = ROOT.TLorentzVector()
#	print 'lepvec = ('+str(lepvec.Px())+','+str(lepvec.Py())+','+str(lepvec.Pz())+','+str(lepvec.E())+')' #DEBUGGING
	for i in range(len(jetVars[0])) :
		thisJet = ROOT.TLorentzVector()
		thisJet.SetPtEtaPhiM(jetVars[0][i],jetVars[1][i],jetVars[2][i],jetVars[3][i])
		thisJetdR = lepvec.DeltaR(thisJet)
		#if the lepton is inside the jet, subtract the lepton fourvector from the jet fourvector
		if thisJetdR < 0.25 :
#			old_pT = thisJet.Pt() #DEBUGGING
#			old_dR = thisJetdR #DEBUGGING
			thisJet = thisJet-lepvec
			thisJetdR = lepvec.DeltaR(thisJet)
#			print '	lepton removed. lepton pT = %.4f, old pT = %.4f, new pT = %.4f, old dR = %.4f, new dR = %.4f'%(lepvec.Pt(),old_pT,thisJet.Pt(),old_dR,thisJetdR) #DEBUGGING
		#check that the jet, subtracted if necessary, has at least the minimum pT
		if thisJet.Pt() < JET_PT_MIN :
			continue
		#set the nearest jet to this jet if necessary
		if thisJetdR < nearest_jet_dR :
			nearest_jet_vec.SetPtEtaPhiM(thisJet.Pt(),thisJet.Eta(),thisJet.Phi(),thisJet.M())
			nearest_jet_dR = thisJetdR
#		print '		nearest_jet_dR = '+str(nearest_jet_dR)+'' #DEBUGGING
	#make lepton 2D cut
	control_plots[9].Fill(nearest_jet_dR,lepvec.Pt(nearest_jet_vec.Vect()))
	if sideband == 0 and (nearest_jet_dR < DR_MIN and lepvec.Pt(nearest_jet_vec.Vect()) < REL_PT_MIN) :
		return -1*CUTFLOW_LEPTON_2D
	if sideband == 1 and not (nearest_jet_dR < DR_MIN and lepvec.Pt(nearest_jet_vec.Vect()) < REL_PT_MIN) :
		return -1*CUTFLOW_LEPTON_2D
	#make electron triangle cut
	if lep_type == 1 and sideband == 0 :
		metE = metvec.E()
		firstJet = ROOT.TLorentzVector() 
		firstJet.SetPtEtaPhiM(jetVars[0][0],jetVars[1][0],jetVars[2][0],jetVars[3][0])
		el_val = abs(lepvec.DeltaPhi(metvec)-1.5); jet_val = abs(firstJet.DeltaPhi(metvec)-1.5)
		control_plots[10].Fill(metE,el_val); control_plots[11].Fill(metE,jet_val)
		cut_val = (1.5/75.)*metE
		if not (el_val < cut_val and jet_val < cut_val) :
			return -1*CUTFLOW_LEPTON_TRIANGLE
	#otherwise return success
	return 0


