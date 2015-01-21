#metHelper: contains helper functions for dealing with the MET/neutrino
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

import ROOT
from math import *

#Global variables
#cutflow
CUTFLOW_MINIMUM_MET = 1
#cut values
MET_PT_MIN = 10. #GeV
#constants
MW = 80.4

#metCut
#takes in the MET variables
#returns 0 if cut is passed or negative value of cutflow failure point
def metCut(metVars,control_plots) :
	metPt = metVars[0][0]
	control_plots[0].Fill(metPt)
	if metPt > MET_PT_MIN :
		return 0
	else :
		return -1*CUTFLOW_MINIMUM_MET


#setupMET
#takes in the fourvector of the selected lepton and the list of 
#variables describing the MET
#returns the a tuple of the two possible fourvectors of the neutrino 
#assuming just the WMass constraint
def setupMET(lep_vec,metVars) :
	met1 = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
	met2 = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
	met1.SetPtEtaPhiM(metVars[0][0],0.0,metVars[1][0],0.0)
	met2.SetPtEtaPhiM(metVars[0][0],0.0,metVars[1][0],0.0)
	pTv    = metVars[0][0]
	phivec = [cos(metVars[1][0]),sin(metVars[1][0])]
	Elep   = lep_vec.E()
	plep   = lep_vec.Vect().Mag()
	pZlep  = lep_vec.Pz()
	pPhi   = lep_vec.Px()*phivec[0]+lep_vec.Py()*phivec[1]
	arg0   = MW*MW+plep*plep-Elep*Elep+2.*pTv*pPhi
	arg    = Elep*Elep*(4.*pTv*pTv*(pZlep*pZlep-Elep*Elep)+arg0*arg0) #discriminant in the quadratic equation solution
#	print ' arg = %.4f = (%.4f)^2*(4*(%.4f)^2*((%.4f)^2-(%.4f)^2)+(%.4f)^2'%(arg,Elep,pTv,pZlep,Elep,arg0) #DEBUGGING
	if not arg > 0 : #If discriminant is imaginary
		pzv1 = pZlep*arg0/(2.*(Elep*Elep-pZlep*pZlep))
		met1.SetPz(pzv1)
		met1.SetE(sqrt(met1.Px()*met1.Px()+met1.Py()*met1.Py()+met1.Pz()*met1.Pz()))
		met2.SetPtEtaPhiM(met1.Pt(),met1.Eta(),met1.Phi(),met1.M())
	else : #have two choices for the neutrino Pz from the quadratic equation
		pzv1 = (pZlep*arg0+sqrt(arg))/(2.*(Elep*Elep-pZlep*pZlep))
		met1.SetPz(pzv1)
		met1.SetE(sqrt(met1.Px()*met1.Px()+met1.Py()*met1.Py()+met1.Pz()*met1.Pz()))
		pzv2 = (pZlep*arg0-sqrt(arg))/(2.*(Elep*Elep-pZlep*pZlep))
		met2.SetPz(pzv2)
		met2.SetE(sqrt(met2.Px()*met2.Px()+met2.Py()*met2.Py()+met2.Pz()*met2.Pz()))
	return (met1,met2)
	#return (ROOT.TLorentzVector(1.0,0.0,0.0,1.0),ROOT.TLorentzVector(1.0,0.0,0.0,1.0)) #DEBUG RETURN