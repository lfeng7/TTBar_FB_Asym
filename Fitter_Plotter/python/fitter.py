#fitter.py is the fitter and plotter workhorse code that fits the MC templates to the data
#and makes and saves a bunch of plots too I guess
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

from ROOT import *
from array import array
import os
import sys
from math import *

##############################         Fitter Class         ##############################

class fitter :
	
	#docstring
	"""Fitter class, fits MC templates to data TTrees, makes plots"""

	#__init__ function
	def __init__(self,runName,templates_filename,output_name) :
		self.run_name = runName
		self.template_file = TFile(templates_filename)
		self.data_tree = self.template_file.Get('data_tree')
		self.template_list = []
		self.final_templates = []
		self.histo_list = []
		self.canv_list = []

	#run theta to fit MC to data
	def fit(self,analysis_filename) :
		print 'running theta from /uscms_data/d3/eminizer/ttbar_run2/CMSSW_7_2_0/src/Analysis/theta/utils2/theta-auto.py'
		print 'analysis script is '+analysis_filename
		os.system('/uscms_data/d3/eminizer/ttbar_run2/CMSSW_7_2_0/src/Analysis/theta/utils2/theta-auto.py '+analysis_filename)

	#makes the file with all of the theta-formatted templates
	def makeTemplateFile(self,filename) :
		#Set some variables
		if self.template_file.GetListOfKeys().Contains('mu__fg0') :
			self.chargeSummed = True; self.lepprefix = 'mu'
		elif self.template_file.GetListOfKeys().Contains('mu__fg0') :
			self.chargeSummed = True; self.lepprefix = 'ele'
		elif self.template_file.GetListOfKeys().Contains('muplus__fg0') and self.template_file.GetListOfKeys().Contains('muminus__fg0') :
			self.chargeSummed = False; self.lepprefix = 'mu'
		elif self.template_file.GetListOfKeys().Contains('eleplus__fg0') and self.template_file.GetListOfKeys().Contains('eleminus__fg0') :
			self.chargeSummed = False; self.lepprefix = 'ele'
		#make all of the MC templates by permuting the parameters and stuff
		self.__organize_templates__()
		for i in range(len(self.final_dists)) :
			new3DTemplate = self.final_dists[i].Clone(self.names[i])
			new3DTemplate.Scale(self.factors[i])
			new3DTemplate.SetDirectory(0)
			self.final_templates.append(new3DTemplate)
		#make the data template(s) by looping over the events in the tree
		self.__loadDataTreeBranches__()
		print 'adding to data template(s)'
		nEntries = self.data_tree.GetEntriesFast()
		print '	# of data events: '+str(nEntries)
		for entry in range(nEntries) :
			percent_done = 100.*entry/nEntries
			if percent_done%10 < 100./nEntries :
				print '	'+str((int)(percent_done))+'%'
			self.data_tree.GetEntry(entry)
			if self.chargeSummed or self.Q_l[0]>0 :
				self.final_data_templates[0].Fill(self.cstar[0],self.x_F[0],self.M[0])
			elif self.Q_l[0]<0 :
				self.final_data_templates[1].Fill(self.cstar[0],self.x_F[0],self.M[0])
		#write out templates to file
		self.outfile = TFile(filename,'recreate')
		for template in self.final_data_templates :
			convertTo1D(template).Write()
		for template in self.final_templates :
			convertTo1D(template).Write()
		self.outfile.Close()

	#makeComparisonPlots function makes the four projection plots and saves them
	def makeComparisonPlots(self,save_pdfs) :
		x_stack = THStack('x_stack','Fit Comparison, X Projection')
		y_stack = THStack('y_stack','Fit Comparison, Y Projection')
		z_stack = THStack('z_stack','Fit Comparison, Z Projection')
		q_stack = THStack('q_stack','Fit Comparison, Lepton Charge')
		#oh my goodness so much annoying code needs to go here

	def __organize_templates__(self) :
		nEntries = self.data_tree.GetEntriesFast()
		#fit parameters and errors
		pars = []
#		pars.append(['Rbck',0.064297,0.064297-0.2,0.064297+0.2])
#		pars.append(['Rntmj',0.224518,0.224518-0.2,0.224518+0.2])
#		pars.append(['Rqqbar',0.057765,0.057765-0.2,0.057765+0.2])
#		pars.append(['Afb',0.000000,-0.2,0.2])
		pars.append(['Rbck',1.0,0.8,1.2])
		pars.append(['Rntmj',1.0,0.8,1.2])
		pars.append(['Rqqbar',1.0,0.8,1.2])
		pars.append(['Afb',0.000000,-0.2,0.2])
		dists, funcs = self.__getDistributions__()
		for dist in dists :
			dist.SetDirectory(0)
		#build the templates
		self.final_dists = []
		#lists of new template prefactors and names
		self.factors = []; self.names = []
		for i in range(len(dists)) : #for each distribution
			func = funcs[i]; args = func.split('#'); argispar = []
			for j in range(len(args)) :
				argispar.append(False)
			for j in range(len(args)) : #replace the parameters within the arguments
				for par in pars :
					if par[0] == args[j] :
						args[j] = par
						argispar[j] = True
			#nominal distribution
			self.names.append(dists[i].GetName())
			factorstring = ''
			for j in range(len(args)) :
				if not argispar[j] :
					factorstring+=args[j]
				else :
					factorstring+=str(args[j][1])
			self.factors.append(eval(factorstring))
			self.final_dists.append(dists[i])
			print 'adding template with name %s and factor %s=%.4f'%(self.names[len(self.names)-1],factorstring,self.factors[len(self.factors)-1])
			#all the up and down distributions
			for k in range(argispar.count(True)) :
				factorstring_down = ''; factorstring_up = ''
				countpars = 0
				for j in range(len(args)) :
					if not argispar[j] :
						factorstring_down+=args[j]; factorstring_up+=args[j]
					elif argispar[j] and countpars!=k :
						factorstring_down+=str(args[j][1]); factorstring_up+=str(args[j][1])
						countpars+=1
					elif argispar[j] and countpars==k :
						factorstring_down+=str(args[j][2]); factorstring_up+=str(args[j][3])
						self.names.append(dists[i].GetName()+'__'+args[j][0]+'__'+'down')
						self.names.append(dists[i].GetName()+'__'+args[j][0]+'__'+'up')
						countpars+=1
				self.factors.append(eval(factorstring_down)); self.factors.append(eval(factorstring_up))
				self.final_dists.append(dists[i]); self.final_dists.append(dists[i])
				print 'adding template with name %s and factor %s=%.4f'%(self.names[len(self.names)-2],factorstring_down,self.factors[len(self.factors)-2])
				print 'adding template with name %s and factor %s=%.4f'%(self.names[len(self.names)-1],factorstring_up,  self.factors[len(self.factors)-1])

	def __getDistributions__(self) :
		dists = []; funcs = []
		if self.chargeSummed :
			fg0 	   = self.template_file.Get(self.lepprefix+'__fg0'); dists.append(fg0); funcs.append('(3.-#Rbck#-#Rntmj#)*(2.-#Rqqbar#)')
#			fg1 	   = self.template_file.Get(self.lepprefix+'__fg1')
#			fg2 	   = self.template_file.Get(self.lepprefix+'__fg2')
#			fg3 	   = self.template_file.Get(self.lepprefix+'__fg3')
#			fg4 	   = self.template_file.Get(self.lepprefix+'__fg4')
			fqs0 	   = self.template_file.Get(self.lepprefix+'__fqs0'); dists.append(fqs0); funcs.append('(3.-#Rbck#-#Rntmj#)*#Rqqbar#')
#			fqs1 	   = self.template_file.Get(self.lepprefix+'__fqs1')
#			fqs2 	   = self.template_file.Get(self.lepprefix+'__fqs2')
			fqa0 	   = self.template_file.Get(self.lepprefix+'__fqa0'); dists.append(fqa0); funcs.append('(3.-#Rbck#-#Rntmj#)*#Rqqbar#*#Afb#')
#			fqa1 	   = self.template_file.Get(self.lepprefix+'__fqa1')
#			fqa2 	   = self.template_file.Get(self.lepprefix+'__fqa2')
			fbck 	   = self.template_file.Get(self.lepprefix+'__fbck'); dists.append(fbck); funcs.append('#Rbck#')
			fntmj 	   = self.template_file.Get(self.lepprefix+'__fntmj'); dists.append(fntmj); funcs.append('#Rntmj#')
		else :
			fg0_plus 	   = self.template_file.Get(self.lepprefix+'plus__fg0'); dists.append(fg0_plus); funcs.append('(3.-#Rbck#-#Rntmj#)*(2.-#Rqqbar#)')
#			fg1_plus 	   = self.template_file.Get(self.lepprefix+'plus__fg1')
#			fg2_plus 	   = self.template_file.Get(self.lepprefix+'plus__fg2')
#			fg3_plus 	   = self.template_file.Get(self.lepprefix+'plus__fg3')
#			fg4_plus 	   = self.template_file.Get(self.lepprefix+'plus__fg4')

			fqs0_plus 	   = self.template_file.Get(self.lepprefix+'plus__fqs0'); dists.append(fqs0_plus); funcs.append('(3.-#Rbck#-#Rntmj#)*#Rqqbar#')
#			fqs1_plus 	   = self.template_file.Get(self.lepprefix+'plus__fqs1')
#			fqs2_plus 	   = self.template_file.Get(self.lepprefix+'plus__fqs2')

			fqa0_plus 	   = self.template_file.Get(self.lepprefix+'plus__fqa0'); dists.append(fqa0_plus); funcs.append('(3.-#Rbck#-#Rntmj#)*#Rqqbar#*#Afb#')
#			fqa1_plus 	   = self.template_file.Get(self.lepprefix+'plus__fqa1')
#			fqa2_plus 	   = self.template_file.Get(self.lepprefix+'plus__fqa2')

			fbck_plus 	   = self.template_file.Get(self.lepprefix+'plus__fbck'); dists.append(fbck_plus); funcs.append('#Rbck#')
			fntmj_plus    = self.template_file.Get(self.lepprefix+'plus__fntmj'); dists.append(fntmj_plus); funcs.append('#Rntmj#')

			fg0_minus 	   = self.template_file.Get(self.lepprefix+'minus__fg0'); dists.append(fg0_minus); funcs.append('(3.-#Rbck#-#Rntmj#)*(2.-#Rqqbar#)')
#			fg1_minus 	   = self.template_file.Get(self.lepprefix+'minus__fg1')
#			fg2_minus 	   = self.template_file.Get(self.lepprefix+'minus__fg2')
#			fg3_minus 	   = self.template_file.Get(self.lepprefix+'minus__fg3')
#			fg4_minus 	   = self.template_file.Get(self.lepprefix+'minus__fg4')

			fqs0_minus    = self.template_file.Get(self.lepprefix+'minus__fqs0'); dists.append(fqs0_minus); funcs.append('(3.-#Rbck#-#Rntmj#)*#Rqqbar#')
#			fqs1_minus    = self.template_file.Get(self.lepprefix+'minus__fqs1')
#			fqs2_minus    = self.template_file.Get(self.lepprefix+'minus__fqs2')

			fqa0_minus    = self.template_file.Get(self.lepprefix+'minus__fqa0'); dists.append(fqa0_minus); funcs.append('(3.-#Rbck#-#Rntmj#)*#Rqqbar#*#Afb#')
#			fqa1_minus    = self.template_file.Get(self.lepprefix+'minus__fqa1')
#			fqa2_minus    = self.template_file.Get(self.lepprefix+'minus__fqa2')

			fbck_minus    = self.template_file.Get(self.lepprefix+'minus__fbck'); dists.append(fbck_minus); funcs.append('#Rbck#')
			fntmj_minus   = self.template_file.Get(self.lepprefix+'minus__fntmj'); dists.append(fntmj_minus); funcs.append('#Rntmj#')

		return dists, funcs

	def __loadDataTreeBranches__(self) :
		#lepton charge
		self.Q_l = array('i',[0]); self.data_tree.SetBranchAddress('Q_l',self.Q_l)
		#kinematic fit chi2
		self.chi2 = array('d',[0.0]); self.data_tree.SetBranchAddress('chi2',self.chi2)
		#cosine(theta)
		self.cstar = array('d',[100.0]); self.data_tree.SetBranchAddress('cstar',self.cstar)
		#Feynman x
		self.x_F = array('d',[100.0]); self.data_tree.SetBranchAddress('x_F',self.x_F)
		#ttbar invariant mass
		self.M = array('d',[-1.0]); self.data_tree.SetBranchAddress('M',self.M)
		#and make the final data templates, why not?
		self.final_data_templates = []
		if self.chargeSummed :
			self.final_data_templates.append(self.final_templates[0].Clone(self.lepprefix+'__DATA'))
		else :
			self.final_data_templates.append(self.final_templates[0].Clone(self.lepprefix+'plus__DATA'))
			self.final_data_templates.append(self.final_templates[0].Clone(self.lepprefix+'minus__DATA'))
		for template in self.final_data_templates :
			template.SetDirectory(0)
			template.Reset()

	#__del__ function
	def __del__(self) :
		self.outfile.cd()
		for histo in self.histo_list :
			histo.Write()
		for canv in self.canv_list :
			canv.Write()
		#self.outfile.Write()
		self.outfile.Close()

def convertTo1D(original) :
	nBins = original.GetNbinsX()*original.GetNbinsY()*original.GetNbinsZ()
	newHisto = TH1F(original.GetName(),original.GetTitle(),nBins,0.,nBins-1.)
	newHisto.SetDirectory(0)
	for k in range(nBins) :
		newHisto.SetBinContent(k,original.GetBinContent(k))
	return newHisto