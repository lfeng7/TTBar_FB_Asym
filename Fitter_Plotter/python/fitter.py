#fitter.py is the fitter and plotter workhorse code that fits the MC templates to the data
#and makes and saves a bunch of plots too I guess
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

import ROOT
from array import array
import os
import sys
from math import *

#global variables
INITIAL_RQQBAR = 0.15
INITIAL_RBCK   = 0.30
INITIAL_AFB    = 0.10
INITIAL_XI 	   = 0.00
INITIAL_DELTA  = 0.00
TOTAL_DATA_FILENAME = 'total_data_tree.root'

#fcn is the fitting function for MC to data
def fcn(npar, deriv, f, par, flag) :
	lnL = 0.
	#loop over all the events in the tree
	#parameters, unambiguously
	ev_Rbck = par[0]; ev_Rqqbar = par[1]; ev_Afb = par[2]; ev_xi = par[3]; ev_delta = par[4];
	for entry in tree.GetEntriesFast() :
		tree.GetEntry(entry)
		#skip the event if the cutflow isn't zero
		if cutflow[0] != 0 :
			continue
		#get the event likelihood depending on whether the analysis is charge summed or separated
		ev_fqqs = 0.; ev_fqqa = 0.; ev_fqqs_delta = 0.; ev_fqqs_xi = 0.; ev_fqqa_delta = 0.; ev_fqqa_xi = 0.; ev_fgg = 0.; ev_fbck = 0.
		if chargeSummed :
			bin = fqqs.FindFixBin(cstar[0],x_F[0],M[0])
			ev_fqqs = fqqs.GetBinContent(bin)
			ev_fqqa = fqqa.GetBinContent(bin)
			ev_fqqs_delta = fqqs_delta.GetBinContent(bin)
			ev_fqqs_xi = fqqs_xi.GetBinContent(bin)
			ev_fqqa_delta = fqqa_delta.GetBinContent(bin)
			ev_fqqa_xi = fqqa_xi.GetBinContent(bin)
			ev_fgg = fgg.GetBinContent(bin)
			ev_fbck = fbck.GetBinContent(bin)
		else :
			bin = fqqs_plus.FindFixBin(cstar[0],x_F[0],M[0])
			if Q_l[0] > 0 :
				ev_fqqs = fqqs_plus.GetBinContent(bin)
				ev_fqqa = fqqa_plus.GetBinContent(bin)
				ev_fqqs_delta = fqqs_delta_plus.GetBinContent(bin)
				ev_fqqs_xi = fqqs_xi_plus.GetBinContent(bin)
				ev_fqqa_delta = fqqa_delta_plus.GetBinContent(bin)
				ev_fqqa_xi = fqqa_xi_plus.GetBinContent(bin)
				ev_fgg = fgg_plus.GetBinContent(bin)
				ev_fbck = fbck_plus.GetBinContent(bin)
			elif Q_l[0] < 0 :
				ev_fqqs = fqqs_minus.GetBinContent(bin)
				ev_fqqa = fqqa_minus.GetBinContent(bin)
				ev_fqqs_delta = fqqs_delta_minus.GetBinContent(bin)
				ev_fqqs_xi = fqqs_xi_minus.GetBinContent(bin)
				ev_fqqa_delta = fqqa_delta_minus.GetBinContent(bin)
				ev_fqqa_xi = fqqa_xi_minus.GetBinContent(bin)
				ev_fgg = fgg_minus.GetBinContent(bin)
				ev_fbck = fbck_minus.GetBinContent(bin)
			else :
				print 'WARNING: EVENT SKIPPED BECAUSE LEPTON CHARGE INFORMATION IS UNAVAILABLE'
		#skip events in places where any pos-def template is unpopulated
		if ev_fbck == 0. or ev_fgg == 0. or ev_fqqs == 0. :
			continue
		ev_L = ev_Rbck*ev_fbck																					#background
		ev_L+= (1.0-ev_Rbck)*(1.0-ev_Rqqbar)*ev_fgg 															#glu-glu
		xi_delta_factor = (1.0/(1.0+ev_delta*F_delta+ev_xi*F_xi))
		ev_L+= (1.0-ev_Rbck)*ev_Rqqbar*xi_delta_factor*(ev_fqqs+ev_delta*ev_fqqs_delta+ev_xi*ev_fqqs_xi) 		#symmetric qqbar
		ev_L+= (1.0-ev_Rbck)*ev_Rqqbar*xi_delta_factor*ev_Afb*(ev_fqqa+ev_delta*ev_fqqa_delta+ev_xi*ev_fqqa_xi) #antisymmetric qqbar
		#make sure event likelihood is positive, then add log to lnL to return
		if not ev_L > 0. :
			print 'PARAMETERS ARE IN A BAD SPOT: NEGATIVE EVENT LIKELIHOOD!'
			ev_L = sys.float_info.epsilon #smallest possible positive value
		lnL+=log(ev_L)
	#set the return value
	f[0] = lnL


##############################         Fitter Class         ##############################

class fitter :
	
	#docstring
	"""Fitter class, fits MC templates to data TTrees, makes plots"""

	#__init__ function
	def __init__(self,runName,onGrid,MC_templates_filename,data_TTree_filename,output_name) :
		self.run_name = runName
		self.template_file = TFile(MC_templates_filename)
		self.__mergeAllDataTTrees__(data_TTree_filename)
		self.f = TFile(output_name,'recreate')
		self.on_grid = onGrid.lower()
		self.histo_list = []
		self.canv_list = []

	#fit function fits MC templates to data
	def fit(self) :
		#Get the MC templates
		#charge summed case
		self.template_list = []
		if self.f.GetListOfKeys().Contains('fqqs') :
			chargeSummed = True
			fqqs 	   = self.f.Get('fqqs'); 	   self.template_list.append(fqqs)
			fqqs_xi    = self.f.Get('fqqs_xi');	   self.template_list.append(fqqs_xi)
			fqqs_delta = self.f.Get('fqqs_delta'); self.template_list.append(fqqs_delta)
			fqqa 	   = self.f.Get('fqqa'); 	   self.template_list.append(fqqa)
			fqqa_xi    = self.f.Get('fqqa_xi');	   self.template_list.append(fqqa_xi)
			fqqa_delta = self.f.Get('fqqa_delta'); self.template_list.append(fqqa_delta)
			fgg 	   = self.f.Get('fgg'); 	   self.template_list.append(fgg)
			fbck 	   = self.f.Get('fbck'); 	   self.template_list.append(fbck)
			F_xi 	= fqqs_xi.Integral()
			F_delta = fqqs_delta.Integral()
		#charge separated case
		if self.f.GetListOfKeys().Contains('fqqs_plus') and self.f.GetListOfKeys().Contains('fqqs_minus') :
			chargeSummed = False
			fqqs_plus 		 = self.f.Get('fqqs_plus');		   self.template_list.append(fqqs_plus)
			fqqs_xi_plus 	 = self.f.Get('fqqs_xi_plus');	   self.template_list.append(fqqs_xi_plus)
			fqqs_delta_plus  = self.f.Get('fqqs_delta_plus');  self.template_list.append(fqqs_delta_plus)
			fqqa_plus 		 = self.f.Get('fqqa_plus');		   self.template_list.append(fqqa_plus)
			fqqa_xi_plus 	 = self.f.Get('fqqa_xi_plus');	   self.template_list.append(fqqa_xi_plus)
			fqqa_delta_plus  = self.f.Get('fqqa_delta_plus');  self.template_list.append(fqqa_delta_plus)
			fgg_plus 		 = self.f.Get('fgg_plus');		   self.template_list.append(fgg_plus)
			fbck_plus 		 = self.f.Get('fbck_plus');		   self.template_list.append(fbck_plus)
			fqqs_minus 		 = self.f.Get('fqqs_minus');	   self.template_list.append(fqqs_minus)
			fqqs_xi_minus 	 = self.f.Get('fqqs_xi_minus');	   self.template_list.append(fqqs_xi_minus)
			fqqs_delta_minus = self.f.Get('fqqs_delta_minus'); self.template_list.append(fqqs_delta_minus)
			fqqa_minus 		 = self.f.Get('fqqa_minus');	   self.template_list.append(fqqa_minus)
			fqqa_xi_minus 	 = self.f.Get('fqqa_xi_minus');	   self.template_list.append(fqqa_xi_minus)
			fqqa_delta_minus = self.f.Get('fqqa_delta_minus'); self.template_list.append(fqqa_delta_minus)
			fgg_minus 		 = self.f.Get('fgg_minus');		   self.template_list.append(fgg_minus)
			fbck_minus 		 = self.f.Get('fbck_minus');	   self.template_list.append(fbck_minus)
			F_xi 	= fqqs_xi_plus.Integral()+fqqs_xi_minus.Integral()
			F_delta = fqqs_delta_plus.Integral()+fqqs_delta_minus.Integral()
		#open the total data ttree to use outside of the class
		tree = self.data_tree
		self.__loadDataTreeBranches__()
		#fit parameters and errors
		self.par_names = []; self.par_names_tex = []; self.pars = []; self.par_errs = []
		par_ini_vals = []; par_mins = []; par_maxs = []; par_fixs = []
		newParam('Rbck',  '$R_{\mathrm{bk}}$',INITIAL_RBCK,  0.5,INITIAL_RBCK,  0.0, 1.0,False)
		newParam('Rqqbar','$R_q\\bar q}$',	  INITIAL_RQQBAR,0.5,INITIAL_RQQBAR,0.0, 1.0,False)
		newParam('Afb',   '$A^{(1)}_{FB}$',	  INITIAL_AFB, 	 0.5,INITIAL_AFB, 	-1.0,1.0,False)
		newParam('xi', 	  '$\\xi$',			  INITIAL_XI, 	 0.5,INITIAL_XI, 	-1.0,1.0,True)
		newParam('delta', '$\\delta$',		  INITIAL_DELTA, 0.5,INITIAL_DELTA, -1.0,1.0,True)
		#set up Minuit for fitting
		minuit = ROOT.TMinuit(len(self.pars))
		minuit.SetFCN(fcn)
		ierflag = ROOT.Long(1)
		arglist = array('d', [-1.0])
		arglist[0] = 100000.
		for i in range(len(self.pars)) :
			minuit.mnparm(i,self.par_names[i],self.pars[i],self.par_errs[i],par_ini_vals[i],par_mins[i],par_maxs[i],ierflag)
			if par_fixs[i] :
				#minuit.mnfixp(i,ierflag)
				minuit.FixParameter(i)
		#minimize
		minuit.mnexcm('MIGRAD', arglist, 1,ierflag)
		if ierflag != 0 :
			print 'PROBLEM IN FIT: error flag = '+str(ierflag)+''
		#get final parameters values and errors
		for i in range(len(self.pars)) :
			tmp = ROOT.Double(1.0);		tmp_err = ROOT.Double(1.0)
			minuit.GetParameter(0,tmp,tmp_err)
			self.pars[i] = tmp;		self.par_errs[i] = tmp_err
		#print out final values
		print 'MIGRAD COMPLETE! FINAL PARAMETER VALUES: '
		for i in range(len(self.pars)) :
			print '	'+self.par_names[i]+' = %.4f +/- %.4f'%(self.pars[i],self.par_errs[i])
		#save final values to normal and LaTeX-formatted summary table files in the directory immediately above
		writeToSummaryFile()

	#makeComparisonPlots function makes the four projection plots and saves them
	def makeComparisonPlots(self,save_pdfs) :
		x_stack = ROOT.THStack('x_stack','Fit Comparison, X Projection')
		y_stack = ROOT.THStack('y_stack','Fit Comparison, Y Projection')
		z_stack = ROOT.THStack('z_stack','Fit Comparison, Z Projection')
		q_stack = ROOT.THStack('q_stack','Fit Comparison, Lepton Charge')
		#oh my goodness so much annoying code needs to go here

	def __loadDataTreeBranches__(self) :
		#cutflow
		cutflow = array('I',[0]); tree.SetBranchAddress('cutflow',cutflow)
		#lepton charge
		Q_l = array('i',[0]); tree.SetBranchAddress('Q_l',Q_l)
		#kinematic fit chi2
		chi2 = array('d',[0.0]); tree.SetBranchAddress('chi2',chi2)
		#cosine(theta)
		cstar = array('d',[100.0]); tree.SetBranchAddress('cstar',cstar)
		#Feynman x
		x_F = array('d',[100.0]); tree.SetBranchAddress('x_F',x_F)
		#ttbar invariant mass
		M = array('d',[-1.0]); tree.SetBranchAddress('M',M)

	#newParam function adds a parameter to a whole bunch of lists and stuff
	def newParam(self,name,tex_name,val,err,ini_val,min_val,max_val,fix) :
		self.par_names.append(name)
		self.par_names_tex.append(tex_name)
		self.pars.append(val)
		self.par_errs.append(err)
		par_ini_vals.append(ini_val)
		par_mins.append(min_val)
		par_maxs.append(max_val)
		par_fixs.append(fix)

	#writeToSummaryFile function writes final fit parameter values to .tex and .txt table files in the above directory
	def writeToSummaryFile(self) :
		refdir = ''
		if self.on_grid == 'yes' :
			refdir+='tardir/'
		else :
			refdir+='../'
		if 'summary_table.tex' in os.listdir(refdir) :
			summary_table_tex_file = open(refdir+'summary_table.tex','r')
			new_file_lines = []
			for i in range(5) : #first part of the headers
				new_file_lines.append(summary_table_tex_file.readline())
			new_line = summary_table_tex_file.readline() #throw this old one out
			n_cols = new_line.count('c|')+1
			new_line = '\\begin{tabular}{|'
			for i in range(n_cols) :
				new_line+='c|'
			new_line+='}\\hline\n'
			new_file_lines.append(new_line) #column format
			new_line = summary_table_tex_file.readline().rstrip('\\\\\n')+'& Fit ('+self.run_name+')  \\\\\n' #new run column
			new_file_lines.append(summary_table_tex_file.readline()) #last line before data
			for i in range(len(self.pars)) :
				new_line = summary_table_tex_file.readline().rstrip('\\\\\n')
				if self.par_errs[i] == 0 :
					new_line+='&	'+self.par_names_tex[i]+'	&	%.3f (fixed)	\\\\\n'%(self.pars[i])
				else :
					new_line+='&	'+self.par_names_tex[i]+'	&	%.3f$\\pm$%.3f	\\\\\n'%(self.pars[i],self.par_errs[i])
				new_file_lines.append(new_line)
			for i in range(4) :
				new_file_lines.append(summary_table_tex_file.readline()) #four footer lines
			#write the lines of the new file back on
			summary_table_tex_file.close()
			summary_table_tex_file = open(refdir+'summary_table.tex','w')
			summary_table_tex_file.writelines(new_file_lines)
		else :
			os.system('touch '+refdir+'summary_table.tex')
			summary_table_tex_file = open(refdir+'summary_table.tex','w')
			summary_table_tex_file.write('\\begin{table}[hbt]\n\\begin{center}\n')
			summary_table_tex_file.write('\\caption{\\small \\label{tab:summary_table} Summary table! Woo-hoo!})\n')
			summary_table_tex_file.write('\\vspace{3pt}\n')
			summary_table_tex_file.write('\\begin{tabular}{|c|c|}\\hline\n')
			summary_table_tex_file.write('Parameter                 & Fit ('+self.run_name+')  \\\\\n')
			summary_table_tex_file.write('\\hline\n')
			for i in range(len(self.pars)) :
				if self.par_errs[i] == 0 :
					summary_table_tex_file.write('&	%.3f (fixed)	\\\\\n'%(self.pars[i]))
				else :
					summary_table_tex_file.write('&	%.3f$\\pm$%.3f	\\\\\n'%(self.pars[i],self.par_errs[i]))
			summary_table_tex_file.write('\\hline\n\\end{tabular}\n\\end{center}\n\\end{table}\n')
		if 'summary_table.txt' in os.listdir(refdir) :
			summary_table_txt_file = open(refdir+'summary_table.txt','r')
			new_file_lines = []
			new_line = summary_table_txt_file.readline().rstrip('\n')+'& Fit ('+self.run_name+')  \n' #new run column
			new_file_lines.append(summary_table_txt_file.readline().rstrip('\n')+'--------------\n') #last line before data
			for i in range(len(self.pars)) :
				new_line = summary_table_txt_file.readline().rstrip('\n')
				if self.par_errs[i] == 0 :
					new_line+='	%.3f (fixed)	\n'%(self.pars[i])
				else :
					new_line+='	%.3f+/-%.3f	\n'%(self.pars[i],self.par_errs[i])
				new_file_lines.append(new_line)
			#write the lines of the new file back on
			summary_table_txt_file.close()
			summary_table_txt_file = open(refdir+'summary_table.txt','w')
			summary_table_txt_file.writelines(new_file_lines)
		else :
			os.system('touch '+refdir+'summary_table.txt')
			summary_table_txt_file = open(refdir+'summary_table.txt','w')
			summary_table_txt_file.write('Parameter                  Fit ('+self.run_name+')  \n')
			summary_table_txt_file.write('------------------------------------------------------------------\n')
			for i in range(len(self.pars)) :
				if self.par_errs[i] == 0 :
					summary_table_txt_file.write(self.par_names[i]+'		%.3f (fixed)	\n'%(self.pars[i]))
				else :
					summary_table_txt_file.write(self.par_names[i]+'		%.3f+/-%.3f	\n'%(self.pars[i],self.par_errs[i]))
		if 'summary.csv' not in os.listdir(refdir) :
			os.system('touch '+refdir+'summary_table.txt')
			summary_csv_file = open(refdir+'summary.csv','w')
			first_line = ''
			for i in range(len(self.pars)) :
				first_line+=par_names[i]+','+par_names[i]+'_err'
				if i<len(self.pars)-1 :
					first_line+=','
			first_line+='\n'
			summary_csv_file.write(first_line)
			summary_csv_file.close()
		summary_csv_file = open(refdir+'summary.csv','a')
		new_line = ''
		for i in range(len(self.pars)) :
			new_line+='%.3f,%.3f'%(pars[i],par_errs[i])
			if i<len(self.pars)-1 :
				new_line+=','
		new_line+='\n'
		summary_csv_file.write(new_line)
		summary_csv_file.close()
			

	#__mergeAllDataTTrees__ function merges the data TTrees by opening the input text file and pulling the trees
	#from each of the data files
	def __mergeAllDataTTrees__(self,input_filename) :
		chain = ROOT.TChain('tree')
		inputfile = open(input_filename,'r')
		for line in inputfile :
			if line.startswith('#') :
				continue
			chain.Add(line.rstrip())
		chain.Merge(TOTAL_DATA_FILENAME)
		self.data_tree = TOTAL_DATA_FILENAME.Get('tree')

	#__del__ function
	def __del__(self) :
		self.f.cd()
		for histo in self.histo_list :
			histo.Write()
		for canv in canv_list :
			canv.Write()
		self.f.Write()
		self.f.Close()