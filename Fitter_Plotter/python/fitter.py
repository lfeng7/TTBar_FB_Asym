#fitter.py is the fitter and plotter workhorse code that fits the MC templates to the data
#and makes and saves a bunch of plots too I guess
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

from ROOT import *
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

##############################         Fitter Class         ##############################

class fitter :
	
	#docstring
	"""Fitter class, fits MC templates to data TTrees, makes plots"""

	#__init__ function
	def __init__(self,runName,onGrid,templates_filename,output_name) :
		self.run_name = runName
		self.template_file = TFile(templates_filename)
		self.f = TFile(output_name,'recreate')
		self.on_grid = onGrid.lower()
		self.histo_list = []
		self.canv_list = []

	#fcn is the fitting function for MC to data
	def fcn(self, npar, deriv, f, par, flag) :
		lnL = 0.
		#loop over all the events in the tree
		#parameters, unambiguously
		ev_Rbck = par[0]; ev_Rqqbar = par[1]; ev_Afb = par[2]; ev_xi = par[3]; ev_delta = par[4];
	#	print 'in function, tree = '+str(tree) #DEBUG
		for entry in range(self.data_tree.GetEntriesFast()) :
			self.data_tree.GetEntry(entry)
			#get the event likelihood depending on whether the analysis is charge summed or separated
			ev_fqqs = 0.; ev_fqqa = 0.; ev_fqqs_delta = 0.; ev_fqqs_xi = 0.; ev_fqqa_delta = 0.; ev_fqqa_xi = 0.; ev_fgg = 0.; ev_fbck = 0.
			if self.chargeSummed :
				bin = self.fqqs.FindFixBin(self.cstar[0],self.x_F[0],self.M[0])
				ev_fqqs = self.fqqs.GetBinContent(bin)
				ev_fqqa = self.fqqa.GetBinContent(bin)
				ev_fqqs_delta = self.fqqs_delta.GetBinContent(bin)
				ev_fqqs_xi = self.fqqs_xi.GetBinContent(bin)
				ev_fqqa_delta = self.fqqa_delta.GetBinContent(bin)
				ev_fqqa_xi = self.fqqa_xi.GetBinContent(bin)
				ev_fgg = self.fgg.GetBinContent(bin)
				ev_fbck = self.fbck.GetBinContent(bin)
			else :
				bin = self.fqqs_plus.FindFixBin(self.cstar[0],self.x_F[0],self.M[0])
				if self.Q_l[0] > 0 :
					ev_fqqs = self.fqqs_plus.GetBinContent(bin)
					ev_fqqa = self.fqqa_plus.GetBinContent(bin)
					ev_fqqs_delta = self.fqqs_delta_plus.GetBinContent(bin)
					ev_fqqs_xi = self.fqqs_xi_plus.GetBinContent(bin)
					ev_fqqa_delta = self.fqqa_delta_plus.GetBinContent(bin)
					ev_fqqa_xi = self.fqqa_xi_plus.GetBinContent(bin)
					ev_fgg = self.fgg_plus.GetBinContent(bin)
					ev_fbck = self.fbck_plus.GetBinContent(bin)
				elif self.Q_l[0] < 0 :
					ev_fqqs = self.fqqs_minus.GetBinContent(bin)
					ev_fqqa = self.fqqa_minus.GetBinContent(bin)
					ev_fqqs_delta = self.fqqs_delta_minus.GetBinContent(bin)
					ev_fqqs_xi = self.fqqs_xi_minus.GetBinContent(bin)
					ev_fqqa_delta = self.fqqa_delta_minus.GetBinContent(bin)
					ev_fqqa_xi = self.fqqa_xi_minus.GetBinContent(bin)
					ev_fgg = self.fgg_minus.GetBinContent(bin)
					ev_fbck = self.fbck_minus.GetBinContent(bin)
				else :
					print 'WARNING: EVENT SKIPPED BECAUSE LEPTON CHARGE INFORMATION IS UNAVAILABLE'
			print 'bin = %d'%(bin)
			#skip events in places where any pos-def template is unpopulated
			if ev_fbck == 0. or ev_fgg == 0. or ev_fqqs == 0. :
				print 'bad phase space point: (%.4f,%.4f,%.4f) in bin %d'%(self.cstar[0],self.x_F[0],self.M[0],bin)
				continue
			ev_L = ev_Rbck*ev_fbck																					#background
			ev_L+= (1.0-ev_Rbck)*(1.0-ev_Rqqbar)*ev_fgg 															#glu-glu
			xi_delta_factor = (1.0/(1.0+ev_delta*self.F_delta+ev_xi*self.F_xi))
			ev_L+= (1.0-ev_Rbck)*ev_Rqqbar*xi_delta_factor*(ev_fqqs+ev_delta*ev_fqqs_delta+ev_xi*ev_fqqs_xi) 		#symmetric qqbar
			ev_L+= (1.0-ev_Rbck)*ev_Rqqbar*xi_delta_factor*ev_Afb*(ev_fqqa+ev_delta*ev_fqqa_delta+ev_xi*ev_fqqa_xi) #antisymmetric qqbar
			#make sure event likelihood is positive, then add log to lnL to return
			if not ev_L > 0. :
				print 'PARAMETERS ARE IN A BAD SPOT: NEGATIVE EVENT LIKELIHOOD!'
				ev_L = sys.float_info.epsilon #smallest possible positive value
			lnL+=log(ev_L)
		#set the return value
	#	print 'chi2 at this iteration: %.4f'%(lnL) #DEBUG
		f[0] = lnL

	#fit function fits MC templates to data
	def fit(self) :
		#Get the MC templates
		#charge summed case
		self.template_list = []
		if self.template_file.GetListOfKeys().Contains('fqqs') :
			self.chargeSummed = True
			self.fqqs 	   = self.template_file.Get('fqqs'); 	   self.template_list.append(self.fqqs)
			self.fqqs_xi    = self.template_file.Get('fqqs_xi');	   self.template_list.append(self.fqqs_xi)
			self.fqqs_delta = self.template_file.Get('fqqs_delta'); self.template_list.append(self.fqqs_delta)
			self.fqqa 	   = self.template_file.Get('fqqa'); 	   self.template_list.append(self.fqqa)
			self.fqqa_xi    = self.template_file.Get('fqqa_xi');	   self.template_list.append(self.fqqa_xi)
			self.fqqa_delta = self.template_file.Get('fqqa_delta'); self.template_list.append(self.fqqa_delta)
			self.fgg 	   = self.template_file.Get('fgg'); 	   self.template_list.append(self.fgg)
			self.fbck 	   = self.template_file.Get('fbck'); 	   self.template_list.append(self.fbck)
			self.fntmj 	   = self.template_file.Get('fntmj'); 	   self.template_list.append(self.fntmj)
			self.F_xi 	= self.fqqs_xi.Integral()
			self.F_delta = self.fqqs_delta.Integral()
		#charge separated case
		if self.template_file.GetListOfKeys().Contains('fqqs_plus') and self.template_file.GetListOfKeys().Contains('fqqs_minus') :
			self.chargeSummed = False
			self.fqqs_plus 		 = self.template_file.Get('fqqs_plus');		   self.template_list.append(self.fqqs_plus)
			self.fqqs_xi_plus 	 = self.template_file.Get('fqqs_xi_plus');	   self.template_list.append(self.fqqs_xi_plus)
			self.fqqs_delta_plus  = self.template_file.Get('fqqs_delta_plus');  self.template_list.append(self.fqqs_delta_plus)
			self.fqqa_plus 		 = self.template_file.Get('fqqa_plus');		   self.template_list.append(self.fqqa_plus)
			self.fqqa_xi_plus 	 = self.template_file.Get('fqqa_xi_plus');	   self.template_list.append(self.fqqa_xi_plus)
			self.fqqa_delta_plus  = self.template_file.Get('fqqa_delta_plus');  self.template_list.append(self.fqqa_delta_plus)
			self.fgg_plus 		 = self.template_file.Get('fgg_plus');		   self.template_list.append(self.fgg_plus)
			self.fbck_plus 		 = self.template_file.Get('fbck_plus');		   self.template_list.append(self.fbck_plus)
			self.fqqs_minus 		 = self.template_file.Get('fqqs_minus');	   self.template_list.append(self.fqqs_minus)
			self.fqqs_xi_minus 	 = self.template_file.Get('fqqs_xi_minus');	   self.template_list.append(self.fqqs_xi_minus)
			self.fqqs_delta_minus = self.template_file.Get('fqqs_delta_minus'); self.template_list.append(self.fqqs_delta_minus)
			self.fqqa_minus 		 = self.template_file.Get('fqqa_minus');	   self.template_list.append(self.fqqa_minus)
			self.fqqa_xi_minus 	 = self.template_file.Get('fqqa_xi_minus');	   self.template_list.append(self.fqqa_xi_minus)
			self.fqqa_delta_minus = self.template_file.Get('fqqa_delta_minus'); self.template_list.append(self.fqqa_delta_minus)
			self.fgg_minus 		 = self.template_file.Get('fgg_minus');		   self.template_list.append(self.fgg_minus)
			self.fbck_minus 		 = self.template_file.Get('fbck_minus');	   self.template_list.append(self.fbck_minus)
			self.F_xi 	= self.fqqs_xi_plus.Integral()+self.fqqs_xi_minus.Integral()
			self.F_delta = self.fqqs_delta_plus.Integral()+self.fqqs_delta_minus.Integral()
		#open the total data ttree to use outside of the class
		self.data_tree = self.template_file.Get('data_tree')
#		print 'self.data_tree = '+str(self.data_tree)+', tree = '+str(tree) #DEBUG
		self.__loadDataTreeBranches__()
		#fit parameters and errors
		self.par_names = []; self.par_names_tex = []; self.pars = []; self.par_errs = []
		self.par_ini_vals = []; self.par_mins = []; self.par_maxs = []; self.par_fixs = []
		self.newParam('Rbck',  '$R_{\mathrm{bk}}$',INITIAL_RBCK,  0.5,INITIAL_RBCK,  0.0, 1.0,False)
		self.newParam('Rqqbar','$R_q\\bar q}$',	  INITIAL_RQQBAR,0.5,INITIAL_RQQBAR,0.0, 1.0,False)
		self.newParam('Afb',   '$A^{(1)}_{FB}$',	  INITIAL_AFB, 	 0.5,INITIAL_AFB, 	-1.0,1.0,False)
		self.newParam('xi', 	  '$\\xi$',			  INITIAL_XI, 	 0.5,INITIAL_XI, 	-1.0,1.0,True)
		self.newParam('delta', '$\\delta$',		  INITIAL_DELTA, 0.5,INITIAL_DELTA, -1.0,1.0,True)
		#set up Minuit for fitting
		minuit = TMinuit(len(self.pars))
		minuit.SetFCN(self.fcn)
		ierflag = Long(1)
		arglist = array('d', [-1.0])
		arglist[0] = 100000.
		for i in range(len(self.pars)) :
			minuit.mnparm(i,self.par_names[i],self.par_ini_vals[i],self.par_errs[i],self.par_mins[i],self.par_maxs[i],ierflag)
			if self.par_fixs[i] :
				#minuit.mnfixp(i,ierflag)
				minuit.FixParameter(i)
		#minimize
		minuit.mnexcm('MIGRAD', arglist, 1,ierflag)
		if ierflag != 0 :
			print 'PROBLEM IN FIT: error flag = '+str(ierflag)+''
		#get final parameters values and errors
		for i in range(len(self.pars)) :
			tmp = Double(1.0);		tmp_err = Double(1.0)
			minuit.GetParameter(0,tmp,tmp_err)
			self.pars[i] = tmp;		self.par_errs[i] = tmp_err
		#print out final values
		print 'MIGRAD COMPLETE! FINAL PARAMETER VALUES: '
		for i in range(len(self.pars)) :
			print '	'+self.par_names[i]+' = %.4f +/- %.4f'%(self.pars[i],self.par_errs[i])
		#save final values to normal and LaTeX-formatted summary table files in the directory immediately above
		self.writeToSummaryFile()

	#makeComparisonPlots function makes the four projection plots and saves them
	def makeComparisonPlots(self,save_pdfs) :
		x_stack = THStack('x_stack','Fit Comparison, X Projection')
		y_stack = THStack('y_stack','Fit Comparison, Y Projection')
		z_stack = THStack('z_stack','Fit Comparison, Z Projection')
		q_stack = THStack('q_stack','Fit Comparison, Lepton Charge')
		#oh my goodness so much annoying code needs to go here

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

	#newParam function adds a parameter to a whole bunch of lists and stuff
	def newParam(self,name,tex_name,val,err,ini_val,min_val,max_val,fix) :
		self.par_names.append(name)
		self.par_names_tex.append(tex_name)
		self.pars.append(val)
		self.par_errs.append(err)
		self.par_ini_vals.append(ini_val)
		self.par_mins.append(min_val)
		self.par_maxs.append(max_val)
		self.par_fixs.append(fix)

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
			new_file_lines.append(new_line)
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
			os.system('touch '+refdir+'summary.csv')
			summary_csv_file = open(refdir+'summary.csv','w')
			first_line = ''
			for i in range(len(self.pars)) :
				first_line+=self.par_names[i]+','+self.par_names[i]+'_err'
				if i<len(self.pars)-1 :
					first_line+=','
			first_line+='\n'
			summary_csv_file.write(first_line)
			summary_csv_file.close()
		summary_csv_file = open(refdir+'summary.csv','a')
		new_line = ''
		for i in range(len(self.pars)) :
			new_line+='%.3f,%.3f'%(self.pars[i],self.par_errs[i])
			if i<len(self.pars)-1 :
				new_line+=','
		new_line+='\n'
		summary_csv_file.write(new_line)
		summary_csv_file.close()

	#__del__ function
	def __del__(self) :
		self.f.cd()
		for histo in self.histo_list :
			histo.Write()
		for canv in self.canv_list :
			canv.Write()
		self.f.Write()
		self.f.Close()