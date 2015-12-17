#angleReconstructor calculates differential cross section observables given event fourvectors
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

import ROOT
from math import *

#Global variables
#Alpha value for adjustment due to longitudinal gluon polarization
ALPHA = -0.129 #This is the value for the 8TeV Powheg sample as far as we can tell. . .
#epsilon value for gg cross section correction
EPSILON = 0.820 #Value from fit
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
	beta = sqrt(1. - 2.*(M2_1+M2_2)/(ttbar_mass*ttbar_mass) + (M2_1-M2_2)*(M2_1-M2_2)/(ttbar_mass*ttbar_mass*ttbar_mass*ttbar_mass))
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
def getMCObservables(q_vec,qbar_vec,t_vec,tbar_vec,eventtype) :
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
	beta = sqrt(1. - 2.*(M2_1+M2_2)/(ttbar_mass*ttbar_mass) + (M2_1-M2_2)*(M2_1-M2_2)/(ttbar_mass*ttbar_mass*ttbar_mass*ttbar_mass))
	#Doing the boost
	R = R.Boost(Bx,By,Bz)
	t_vec = R*t_vec; tbar_vec = R*tbar_vec
	q_vec = R*q_vec; qbar_vec = R*qbar_vec
	#Reset the boost
	R=S
	#Define normalized three-vectors for the top and protons in the ttbar rest frame
	top = t_vec.Vect(); q = q_vec.Vect(); qbar = qbar_vec.Vect()
	#Normalize vectors (and flip the qbar direction)
	if eventtype != 0 and q.Mag()>qbar.Mag() :
		top = top*(1.0/top.Mag()); q = -1.0*q*(1.0/q.Mag()); qbar = qbar*(1.0/qbar.Mag())
	else :
		top = top*(1.0/top.Mag()); q = q*(1.0/q.Mag()); qbar = -1.0*qbar*(1.0/qbar.Mag())
	#find the unit bisectors
	bisector = (q+qbar)*(1.0/(q+qbar).Mag())
	#find the CS angle
	cos_theta_cs=cos(top.Angle(bisector))
	#calculate the reweighting factors
	b2 = beta*beta;
	b2c2 = beta*beta*cos_theta_cs*cos_theta_cs;
	denom = 1.+b2c2+(1.-b2)+ALPHA*(1.-b2c2)
	Afunc = (7.+9.*b2c2)/(1.-b2c2)
	Bfunc = (1.-b2c2*b2c2+2.*b2*(1.-b2)*(1.-cos_theta_cs*cos_theta_cs))/(2.*(1.-b2c2))
	wg1 = 1./(Bfunc*(1.+EPSILON*b2c2))
	wg2 = 56./((1.-b2)*Afunc*Bfunc*(1.+EPSILON*b2c2))
	wg3 = 4./((1.-b2c2)*Afunc*Bfunc*(1.+EPSILON*b2c2))
	wg4 = (8./(Afunc*Bfunc*(1.+EPSILON*b2c2)))*(1./(1.-b2c2) + 1./(1.-b2) + ((4.*(1.-b2c2))/((1.-b2)*(1.-b2))))
	wqs1 = (4./denom)
	wqs2 = (4./denom)*(1.-b2c2)/(1.-b2)
	wqa0 = (2.*(2.+ALPHA)*(1.-b2/3.)*cos_theta_cs)/denom
	wqa1 = (8.*cos_theta_cs)/denom
	wqa2 = ((8.*cos_theta_cs)/denom)*((1.-b2/3.)/(1.-b2))
	wg1_opp = 1./(Bfunc*(1.+EPSILON*b2c2))
	wg2_opp = 56./((1.-b2)*Afunc*Bfunc*(1.+EPSILON*b2c2))
	wg3_opp = 4./((1.-b2c2)*Afunc*Bfunc*(1.+EPSILON*b2c2))
	wg4_opp = (8./(Afunc*Bfunc*(1.+EPSILON*b2c2)))*(1./(1.-b2c2) + 1./(1.-b2) + ((4.*(1.-b2c2))/((1.-b2)*(1.-b2))))
	wqs1_opp = (4./denom)
	wqs2_opp = (4./denom)*(1.-b2c2)/(1.-b2)
	wqa0_opp = (2.*(2.+ALPHA)*(1.-b2/3.)*(-1.*cos_theta_cs))/denom
	wqa1_opp = (-1.*8.*cos_theta_cs)/denom
	wqa2_opp = ((-1.*8.*cos_theta_cs)/denom)*((1.-b2/3.)/(1.-b2))
	return (cos_theta_cs,x_f,ttbar_mass,wg1,wg2,wg3,wg4,wqs1,wqs2,wqa0,wqa1,wqa2,
		wg1_opp,wg2_opp,wg3_opp,wg4_opp,wqs1_opp,wqs2_opp,wqa0_opp,wqa1_opp,wqa2_opp)
	#return (0.0,0.2,450.,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0) #DEBUG RETURN