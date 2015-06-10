#ttbar Reconstructor contains helper functions to reconstruct the ttbar system
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

import ROOT
from math import *
from array import array

#Global variables
MW = 80.4 #W mass
MT = 172.5 #top mass
QW = (MW*MW)/(MT*MT)
ZW = (2.0*2.0)/(MT*MT)
ZT = (1.4*1.4)/(MT*MT)
CDFT  = (1.0+atan(1./sqrt(ZT)))/sqrt(ZT)
CDFW  = 0.5+2.*QW+(1.5*QW*QW-0.5*ZW*QW-1.5)*log(((1.-QW)*(1.-QW)+ZW*QW)/(QW*QW+QW*ZW))
CDFW += ((QW*QW*QW-3.*ZW*QW*QW-3.*QW+2.)/sqrt(ZW*QW))*(atan((1.-QW)/sqrt(ZW*QW))+atan(QW/sqrt(ZW*QW)))
SIGMAJ  = 0.10 #jet momentum resolution
SIGMAL  = 0.03 #lepton momentum resolution
#global fourvectors for the fit
lep_global_vec = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
met_global_vec = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
blep_global_vec = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
thad_global_vec = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)

#reconstruct takes in the lepton and met fourvectors and the list of jet tuples
#	with fourvectors and CSV and tau_n values
#returns a 5-tuple of:
#	1) the corrected lepton fourvector
#	2) the corrected (and selected) met fourvector
#	3) the corrected leptonic b jet fourvector
#	4) the corrected hadronic top jet fourvector
#	5) the final Chi2 value from the kinematic fit
def reconstruct(lepton,met1,met2,jetlist) :
	#lists of final parameters and chi2 values
	bestParValues = [[met1.Pz(),1.0,1.0,1.0],[met2.Pz(),1.0,1.0,1.0]]
	parNames = ['pZv','scalelep','scaleblep','scaletop']
	parerrs = [0.0,0.0,0.0,0.0]
	finalChi2s = [1000000000.,1000000000.]
	#see whether we have to fit twice based on how many solutions for the neutrino Pz we have
	nFits = 2
	if met1.Pz() == met2.Pz() :
		nFits = 1
	#fit setup stuff common to both iterations
	minuit = ROOT.TMinuit(4)
	minuit.SetFCN(fcn)
	ierflag = ROOT.Long(1)
	arglist = array( 'd', [-1.0] )
	minuit.mnexcm('SET PRINT', arglist, 1,ierflag)
	minuit.mnexcm('SET NOWARNINGS',arglist,1,ierflag)
	arglist[0] = 100000.
	#perform the fit once for each unique neutrino solution
	for iFit in range(nFits) :
		#Set which met solution we're looking at
		if iFit == 0 :
			met = met1
		elif iFit == 1 :
			met = met2
		#set the parameters in minuit
		for i in range(len(bestParValues[0])) :
			minuit.mnparm(i,parNames[i],bestParValues[iFit][i],1.0,0,0,ierflag)
		#set the global fourvector variables
		lep_global_vec.SetPtEtaPhiM(lepton.Pt(),lepton.Eta(),lepton.Phi(),lepton.M())
		met_global_vec.SetPtEtaPhiM(met.Pt(),met.Eta(),met.Phi(),met.M())
		blep_global_vec.SetPtEtaPhiM(jetlist[1].vec.Pt(),jetlist[1].vec.Eta(),
			jetlist[1].vec.Phi(),jetlist[1].vec.M())
		thad_global_vec.SetPtEtaPhiM(jetlist[0].vec.Pt(),jetlist[0].vec.Eta(),
			jetlist[0].vec.Phi(),jetlist[0].vec.M())
		#minimize
		minuit.mnexcm('MIGRAD', arglist, 1,ierflag)
		if ierflag != 0 :
			print 'PROBLEM IN FIT: ierflag = '+str(ierflag)+''
			continue
		#Get the best parameters back from minuit
		for i in range(len(bestParValues[0])) :
			tmp = ROOT.Double(1.0)
			minuit.GetParameter(i,tmp,ROOT.Double(parerrs[i]))
			bestParValues[iFit][i] = tmp
		#Set fit Chi2 for this pZ solution
		tmp1 = ROOT.Double(1.0); tmp2 = ROOT.Double(1.0); tmp3 = ROOT.Double(1.0)
		minuit.mnstat(tmp1,tmp2,tmp3,ROOT.Long(1),ROOT.Long(1),ROOT.Long(1))
		finalChi2s[iFit] = tmp1
	#find which pZ solution gave better results and record best parameter values
#	print 'finalChi2s = '+str(finalChi2s)+'' #DEBUGGING
	final_par_vals = []
	final_met = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
	if finalChi2s[0] < finalChi2s[1] :
		for i in range(len(bestParValues[0])) :
			final_par_vals.append(bestParValues[0][i])
		final_met.SetPtEtaPhiM(met1.Pt(),met1.Eta(),met1.Phi(),met1.M())
	else :
		for i in range(len(bestParValues[1])) :
			final_par_vals.append(bestParValues[1][i])
		final_met.SetPtEtaPhiM(met2.Pt(),met2.Eta(),met2.Phi(),met2.M())
		finalChi2s[0] = finalChi2s[1]
	#rescale the lepton and jet four vectors based on the final parameters
	lep_return = rescale(lepton,final_par_vals[1])
	lepb_return = rescale(jetlist[1],final_par_vals[2])
	hadt_return = rescale(jetlist[0],final_par_vals[3])
	#rebuild the neutrino post-rescaling
	newmetx = final_met.Px()+ (1.0-final_par_vals[1])*lep_return.Px()
	newmety = final_met.Py()+ (1.0-final_par_vals[1])*lep_return.Py()
	newmetx += (1.0-final_par_vals[2])*lepb_return.Px()
	newmety += (1.0-final_par_vals[2])*lepb_return.Py()
	newmetx += (1.0-final_par_vals[3])*hadt_return.Px()
	newmety += (1.0-final_par_vals[3])*hadt_return.Py()
	final_met.SetPx(newmetx); final_met.SetPy(newmety); final_met.SetPz(final_par_vals[0])
	final_met.SetE(final_met.Vect().Mag())
	#return everything
	return (lep_return,final_met,lepb_return,hadt_return,finalChi2s[0])

#top fitting function
def fcn(npar, deriv, f, par, flag) :
	#Build rescaled versions of the vectors involved
	l = rescale(lep_global_vec,par[1])
	bl = rescale(blep_global_vec,par[2])
	th = rescale(thad_global_vec,par[3])
	#rebuild the neutrino from the met post-rescaling
	newmetx = ( met_global_vec.Px()+(1.0-par[1])*lep_global_vec.Px()+(1.0-par[2])*blep_global_vec.Px()
					+(1.0-par[3])*thad_global_vec.Px() )
	newmety = ( met_global_vec.Py()+(1.0-par[1])*lep_global_vec.Py()+(1.0-par[2])*blep_global_vec.Py()
					+(1.0-par[3])*thad_global_vec.Py() )
	v = rescale(met_global_vec,1.0)
	v.SetPx(newmetx); v.SetPy(newmety); v.SetPz(par[0])
	v.SetE(v.Vect().Mag())
	wl = v + l; tl = wl + bl
	mwl2 = wl.M2(); mtl2 = tl.M2(); mth2 = th.M2();
	ql = mwl2/(MT*MT); xl = mtl2/(MT*MT); xh = mth2/(MT*MT)
	pdftl = 1./((xl - 1.)*(xl - 1.) + ZT)
	pdfth = 1./((xh - 1.)*(xh - 1.) + ZT)
	pdfwl = (1. - ql)*(1. - ql)*(2. + ql)/((ql-QW)*(ql-QW)+ZW*QW)
	pdf = pdftl*pdfth*pdfwl/(CDFT*CDFT*CDFW)
	lnL = 0.
	if pdf > 0.0 :
		lnL += log(pdf)    #need positive f
	else :
		print('WARNING -- pdf is negative!!!')
		pdf = 1.e-50
		lnL += log(pdf)
	f[0] = ( -2.0*lnL+(par[1]-1.)*(par[1]-1.)/(SIGMAL*SIGMAL)+(par[2]-1.)*(par[2]-1.)/(SIGMAJ*SIGMAJ)
				+(par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ) )
	#and we don't need to return anything because minuit

#fourvector rescaling function
def rescale(vec,fac) :
	p2 = fac*fac*vec.Vect().Mag2()
	m2 = vec.M()*vec.M()
	newE = sqrt(p2+m2)
	#note that the function returns a NEW TLorentzVector, meaning the original is unaltered
	return ROOT.TLorentzVector(fac*vec.Px(),fac*vec.Py(),fac*vec.Pz(),newE)