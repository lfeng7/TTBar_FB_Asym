#angleReconstructor calculates differential cross section observables given event fourvectors
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

import ROOT
from math import *

#Global variables
#Alpha value for adjustment due to longitudinal gluon polarization
ALPHA = -0.129 #This is the value for the 8TeV Powheg sample as far as we can tell. . .
#Default Lorentz Rotation
S = ROOT.TLorentzRotation()
#Beam energy
SQRT_S=8000.0
BEAM_ENERGY=SQRT_S/2.0

#getObservables takes in the reconstructed top quark vectors and lepton charge
#returns a 3-tuple of:
#	1) costheta
#	2) feynman x
#	3) ttbar mass
def getObservables(lept_vec,hadt_vec,lepton_charge) :
	#find which of the leptonic/hadronic top vectors is the t/tbar
	if lepton_charge == 1 :
		Top = lept_vec
		ATop = hadt_vec
	elif lepton_charge == -1 :
		ATop = lept_vec
		Top = hadt_vec
	#initialize the rotation to the identity
	R = ROOT.TLorentzRotation()
	#Make the 4-vector of the ttbar pair, get its mass, calculate x_F
	Q = Top+ATop
	ttbar_mass=Q.Mag()
	x_f = 2*Q.Pz()/SQRT_S
	#defining the Px, Py,and Pz, and energies to boost into the ttbar rest frame
	Bx = -1*Q.Px()/Q.E(); By = -1*Q.Py()/Q.E(); Bz = -1*Q.Pz()/Q.E()
	#calculating beta for the boost
	M2_1 = Top.Mag2(); M2_2 = ATop.Mag2()
	num  = ( 1. - 2.*(M2_1+M2_2)/(ttbar_mass*ttbar_mass) + 
		(M2_1-M2_2)*(M2_1-M2_2)/(ttbar_mass*ttbar_mass*ttbar_mass*ttbar_mass) )
	denom_1 = (1. + (M2_1-M2_2)/(ttbar_mass*ttbar_mass))*(1. + (M2_1-M2_2)/(ttbar_mass*ttbar_mass))
	denom_2 = (1. + (M2_2-M2_1)/(ttbar_mass*ttbar_mass))*(1. + (M2_2-M2_1)/(ttbar_mass*ttbar_mass))
	beta = sqrt(sqrt((num*num)/(denom_1*denom_1) * (num*num)/(denom_2*denom_2)))
	#Doing the boost
	R = R.Boost(Bx,By,Bz)
	Top = R*Top; ATop = R*ATop
	Proton1 = R*ROOT.TLorentzVector(0.0,0.0,sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY)
	Proton2 = R*ROOT.TLorentzVector(0.0,0.0,-1.0*sqrt(BEAM_ENERGY*BEAM_ENERGY -1*1),BEAM_ENERGY)
	#Reset the boost
	R=S
	#Define three-vectors for the top and protons in the ttbar rest frame
	top = Top.Vect(); proton1 = Proton1.Vect(); proton2 = Proton2.Vect()
	#Flip the larger one between proton1 and proton2
	if proton1.Mag()>proton2.Mag() :
		proton1=-1.0*proton1
	else :
		proton2=-1.0*proton2
	#Normalize vectors
	top = top*(1.0/top.Mag()); proton1 = proton1*(1.0/proton1.Mag()); proton2 = proton2*(1.0/proton2.Mag())
	#find the unit bisectors
	bisector = (proton1+proton2)*(1.0/(proton1+proton2).Mag())
#	print 'bisector = ('+str(bisector.X())+','+str(bisector.Y())+','+str(bisector.Z())+')' #DEBUGGING
	#find the CS angle
	cos_theta_cs=cos(top.Angle(bisector))
	return (cos_theta_cs,x_f,ttbar_mass)
	#return (0.0,0.2,450.) #DEBUG RETURN

#getMCObservables takes in the MC TRUTH initial parton, t, and tbar fourvectors
#returns a 13-tuple of:
#	1-3) MC truth costheta, feynman x, and ttbar mass
#	4-8) antisymmetric, symmetric/antisymmetric xi, and symmetric/antisymmetic delta reweighting factors
#	9-13) the same reweighting factors calculated with the opposite sign angle
def getMCObservables(q_vec,qbar_vec,t_vec,tbar_vec) :
	#initialize the rotation to the identity
	R = ROOT.TLorentzRotation()
	#Make the 4-vector of the ttbar pair, get its mass, calculate x_F
	Q = t_vec+tbar_vec
	ttbar_mass=Q.Mag()
	x_f = 2*Q.Pz()/SQRT_S
	#defining the Px, Py,and Pz, and energies to boost into the ttbar rest frame
	Bx = -1*Q.Px()/Q.E(); By = -1*Q.Py()/Q.E(); Bz = -1*Q.Pz()/Q.E()
	#calculating beta for the boost
	M2_1 = t_vec.Mag2(); M2_2 = tbar_vec.Mag2()
	num  = ( 1. - 2.*(M2_1+M2_2)/(ttbar_mass*ttbar_mass) + 
		(M2_1-M2_2)*(M2_1-M2_2)/(ttbar_mass*ttbar_mass*ttbar_mass*ttbar_mass) )
	denom_1 = (1. + (M2_1-M2_2)/(ttbar_mass*ttbar_mass))*(1. + (M2_1-M2_2)/(ttbar_mass*ttbar_mass))
	denom_2 = (1. + (M2_2-M2_1)/(ttbar_mass*ttbar_mass))*(1. + (M2_2-M2_1)/(ttbar_mass*ttbar_mass))
	beta = sqrt(sqrt((num*num)/(denom_1*denom_1) * (num*num)/(denom_2*denom_2)))
	#Doing the boost
	R = R.Boost(Bx,By,Bz)
	t_vec = R*t_vec; tbar_vec = R*tbar_vec
	q_vec = R*q_vec; qbar_vec = R*qbar_vec
	#Reset the boost
	R=S
	#Define normalized three-vectors for the top and protons in the ttbar rest frame
	top = t_vec.Vect(); q = q_vec.Vect(); qbar = qbar_vec.Vect()
	#Normalize vectors (and flip the qbar direction)
	top = top*(1.0/top.Mag()); q = q*(1.0/q.Mag()); qbar = -1.0*qbar*(1.0/qbar.Mag())
	#find the unit bisectors
	bisector = (q+qbar)*(1.0/(q+qbar).Mag())
	#find the CS angle
	cos_theta_cs=cos(top.Angle(bisector))
	#calculate the reweighting factors
	one_m_b2 = 1.0-beta*beta;
	b2c2 = beta*beta*cos_theta_cs*cos_theta_cs;
	otb2 = (1.0/3.0)*beta*beta;
	denom = 1.0+b2c2+one_m_b2+ALPHA*(1.0-b2c2);
	w_a = 2.0 * ((1.0+otb2+one_m_b2+ALPHA*(1.0-otb2))/denom) * cos_theta_cs; 
	w_s_xi = one_m_b2/denom;
	w_a_xi = 2.0*(one_m_b2/denom)*cos_theta_cs;
	w_s_delta = (1.0-b2c2)/denom;
	w_a_delta = 2.0*((1.0-otb2)/denom)*cos_theta_cs;
	w_a_opp = 2.0 * ((1.0+otb2+one_m_b2+ALPHA*(1.0-otb2))/denom) * (-1.0*cos_theta_cs); 
	w_s_xi_opp = one_m_b2/denom;
	w_a_xi_opp = 2.0*(one_m_b2/denom)*(-1.0*cos_theta_cs);
	w_s_delta_opp = (1.0-b2c2)/denom;
	w_a_delta_opp = 2.0*((1.0-otb2)/denom)*(-1.0*cos_theta_cs);
	return (cos_theta_cs,x_f,ttbar_mass,w_a,w_s_xi,w_a_xi,w_s_delta,w_a_delta,
		w_a_opp,w_s_xi_opp,w_a_xi_opp,w_s_delta_opp,w_a_delta_opp)
	#return (0.0,0.2,450.,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0) #DEBUG RETURN