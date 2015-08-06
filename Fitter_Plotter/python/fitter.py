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
INITIAL_D 	   = 0.00
INITIAL_MU 	   = 0.00

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
		ev_Rbck = par[0]; ev_Rntmj = par[1]; ev_Rqqbar = par[2]; ev_Afb = par[3]; ev_d = par[4]; ev_mu = par[5];
	#	print 'in function, tree = '+str(tree) #DEBUG
		for entry in range(self.data_tree.GetEntriesFast()) :
			self.data_tree.GetEntry(entry)
			#get the event likelihood depending on whether the analysis is charge summed or separated
			ev_fg0 = 0.; ev_fg1 = 0.; ev_fg2 = 0.; ev_fg3 = 0.; ev_fg4 = 0.; ev_fqs0 = 0.; ev_fqs1 = 0.; ev_fqs2 = 0.
			ev_fqa0 = 0.; ev_fqa1 = 0.; ev_fqa2 = 0.; ev_fbck = 0.; ev_fntmj = 0.
			if self.chargeSummed :
				bin = self.fg0.FindFixBin(self.cstar[0],self.x_F[0],self.M[0])
				ev_fg0 = self.fg0.GetBinContent(bin)
				ev_fg1 = self.fg1.GetBinContent(bin)
				ev_fg2 = self.fg2.GetBinContent(bin)
				ev_fg3 = self.fg3.GetBinContent(bin)
				ev_fg4 = self.fg4.GetBinContent(bin)
				ev_fqs0 = self.fqs0.GetBinContent(bin)
				ev_fqs1 = self.fqs1.GetBinContent(bin)
				ev_fqs2 = self.fqs2.GetBinContent(bin)
				ev_fqa0 = self.fqa0.GetBinContent(bin)
				ev_fqa1 = self.fqa1.GetBinContent(bin)
				ev_fqa2 = self.fqa2.GetBinContent(bin)
				ev_fbck = self.fbck.GetBinContent(bin)
				ev_fntmj = self.fntmj.GetBinContent(bin)
			else :
				bin = self.fg0_plus.FindFixBin(self.cstar[0],self.x_F[0],self.M[0])
				if self.Q_l[0] > 0 :
					ev_fg0 = self.fg0_plus.GetBinContent(bin)
					ev_fg1 = self.fg1_plus.GetBinContent(bin)
					ev_fg2 = self.fg2_plus.GetBinContent(bin)
					ev_fg3 = self.fg3_plus.GetBinContent(bin)
					ev_fg4 = self.fg4_plus.GetBinContent(bin)
					ev_fqs0 = self.fqs0_plus.GetBinContent(bin)
					ev_fqs1 = self.fqs1_plus.GetBinContent(bin)
					ev_fqs2 = self.fqs2_plus.GetBinContent(bin)
					ev_fqa0 = self.fqa0_plus.GetBinContent(bin)
					ev_fqa1 = self.fqa1_plus.GetBinContent(bin)
					ev_fqa2 = self.fqa2_plus.GetBinContent(bin)
					ev_fbck = self.fbck_plus.GetBinContent(bin)
					ev_fntmj = self.fntmj_plus.GetBinContent(bin)
				elif self.Q_l[0] < 0 :
					ev_fg0 = self.fg0_minus.GetBinContent(bin)
					ev_fg1 = self.fg1_minus.GetBinContent(bin)
					ev_fg2 = self.fg2_minus.GetBinContent(bin)
					ev_fg3 = self.fg3_minus.GetBinContent(bin)
					ev_fg4 = self.fg4_minus.GetBinContent(bin)
					ev_fqs0 = self.fqs0_minus.GetBinContent(bin)
					ev_fqs1 = self.fqs1_minus.GetBinContent(bin)
					ev_fqs2 = self.fqs2_minus.GetBinContent(bin)
					ev_fqa0 = self.fqa0_minus.GetBinContent(bin)
					ev_fqa1 = self.fqa1_minus.GetBinContent(bin)
					ev_fqa2 = self.fqa2_minus.GetBinContent(bin)
					ev_fbck = self.fbck_minus.GetBinContent(bin)
					ev_fntmj = self.fntmj_minus.GetBinContent(bin)
				else :
					print 'WARNING: EVENT SKIPPED BECAUSE LEPTON CHARGE INFORMATION IS UNAVAILABLE'
			print 'bin = %d'%(bin)
			#skip events in places where any pos-def template is unpopulated
			if ev_fbck == 0. or ev_fg0 == 0. or ev_fqs0 == 0. :
				print 'bad phase space point: (%.4f,%.4f,%.4f) in bin %d'%(self.cstar[0],self.x_F[0],self.M[0],bin)
				continue
			m2pd2 = (ev_mu*ev_mu+ev_d*ev_d)
			FGG = 1. + ev_mu*(1.+ev_mu)*self.F_g1 + m2pd2*(1.+ev_mu)*self.F_g2
			FGG+= m2pd2*(1.-5.*ev_mu)*self.F_g3 + m2pd2*m2pd2*self.F_g4
			FQQBAR = 1. + (2*ev_mu+ev_mu*ev_mu-ev_d*ev_d)*self.F_qs1 + m2pd2*m2pd2*self.F_qs2
			ev_L_func_gg = ev_fg0 + ev_mu*(1.+ev_mu)*ev_fg1 + m2pd2*(1.+ev_mu)*ev_fg2 + m2pd2*(1.-5*ev_mu)*ev_fg3 + m2pd2*m2pd2*ev_fg4
			ev_L_func_qs = ev_fqs0 + (2.*ev_mu+ev_mu*ev_mu-ev_d*ev_d)*ev_fqs1 + m2pd2*ev_fqs2
			ev_L_func_qa = ev_fqa0 + (2.*ev_mu+ev_mu*ev_mu-ev_d*ev_d)*ev_fqa1 + m2pd2*ev_fqa2
			ev_L = ev_Rbck*ev_fbck																					#background
			ev_L+= ev_Rntmj*ev_fntmj																				#NTMJ background
			ev_L+= (1.0-ev_Rbck-ev_Rntmj)*(1.0-ev_Rqqbar)*(1./FGG)*ev_L_func_gg 									#glu-glu
			ev_L+= (1.0-ev_Rbck-ev_Rntmj)*ev_Rqqbar*(1./FQQBAR)*ev_L_func_qs 										#symmetric qqbar
			ev_L+= (1.0-ev_Rbck-ev_Rntmj)*ev_Rqqbar*(1./FQQBAR)*ev_Afb*ev_L_func_qa 								#antisymmetric qqbar
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
		if self.template_file.GetListOfKeys().Contains('fg0') :
			self.chargeSummed = True
			self.fg0 	   = self.template_file.Get('fg0'); 	   self.template_list.append(self.fg0)
			self.fg1 	   = self.template_file.Get('fg1'); 	   self.template_list.append(self.fg1)
			self.fg2 	   = self.template_file.Get('fg2'); 	   self.template_list.append(self.fg2)
			self.fg3 	   = self.template_file.Get('fg3'); 	   self.template_list.append(self.fg3)
			self.fg4 	   = self.template_file.Get('fg4'); 	   self.template_list.append(self.fg4)

			self.fqs0 	   = self.template_file.Get('fqs0'); 	   self.template_list.append(self.fqs0)
			self.fqs1 	   = self.template_file.Get('fqs1'); 	   self.template_list.append(self.fqs1)
			self.fqs2 	   = self.template_file.Get('fqs2'); 	   self.template_list.append(self.fqs2)

			self.fqa0 	   = self.template_file.Get('fqa0'); 	   self.template_list.append(self.fqa0)
			self.fqa1 	   = self.template_file.Get('fqa1'); 	   self.template_list.append(self.fqa1)
			self.fqa2 	   = self.template_file.Get('fqa2'); 	   self.template_list.append(self.fqa2)

			self.fbck 	   = self.template_file.Get('fbck'); 	   self.template_list.append(self.fbck)
			self.fntmj 	   = self.template_file.Get('fntmj'); 	   self.template_list.append(self.fntmj)

			self.F_g1 = self.fg1.Integral()
			self.F_g2 = self.fg2.Integral() 
			self.F_g3 = self.fg3.Integral() 
			self.F_g4 = self.fg4.Integral() 

			self.F_qs1 = self.fqs1.Integral()
			self.F_qs2 = self.fqs2.Integral()

		#charge separated case
		if self.template_file.GetListOfKeys().Contains('fg0_plus') and self.template_file.GetListOfKeys().Contains('fg0_minus') :
			self.chargeSummed = False
			self.fg0_plus 	   = self.template_file.Get('fg0_plus'); 	   self.template_list.append(self.fg0_plus)
			self.fg1_plus 	   = self.template_file.Get('fg1_plus'); 	   self.template_list.append(self.fg1_plus)
			self.fg2_plus 	   = self.template_file.Get('fg2_plus'); 	   self.template_list.append(self.fg2_plus)
			self.fg3_plus 	   = self.template_file.Get('fg3_plus'); 	   self.template_list.append(self.fg3_plus)
			self.fg4_plus 	   = self.template_file.Get('fg4_plus'); 	   self.template_list.append(self.fg4_plus)

			self.fqs0_plus 	   = self.template_file.Get('fqs0_plus'); 	   self.template_list.append(self.fqs0_plus)
			self.fqs1_plus 	   = self.template_file.Get('fqs1_plus'); 	   self.template_list.append(self.fqs1_plus)
			self.fqs2_plus 	   = self.template_file.Get('fqs2_plus'); 	   self.template_list.append(self.fqs2_plus)

			self.fqa0_plus 	   = self.template_file.Get('fqa0_plus'); 	   self.template_list.append(self.fqa0_plus)
			self.fqa1_plus 	   = self.template_file.Get('fqa1_plus'); 	   self.template_list.append(self.fqa1_plus)
			self.fqa2_plus 	   = self.template_file.Get('fqa2_plus'); 	   self.template_list.append(self.fqa2_plus)

			self.fbck_plus 	   = self.template_file.Get('fbck_plus'); 	   self.template_list.append(self.fbck_plus)
			self.fntmj_plus    = self.template_file.Get('fntmj_plus'); 	   self.template_list.append(self.fntmj_plus)

			self.fg0_minus 	   = self.template_file.Get('fg0_minus'); 	   self.template_list.append(self.fg0_minus)
			self.fg1_minus 	   = self.template_file.Get('fg1_minus'); 	   self.template_list.append(self.fg1_minus)
			self.fg2_minus 	   = self.template_file.Get('fg2_minus'); 	   self.template_list.append(self.fg2_minus)
			self.fg3_minus 	   = self.template_file.Get('fg3_minus'); 	   self.template_list.append(self.fg3_minus)
			self.fg4_minus 	   = self.template_file.Get('fg4_minus'); 	   self.template_list.append(self.fg4_minus)

			self.fqs0_minus    = self.template_file.Get('fqs0_minus'); 	   self.template_list.append(self.fqs0_minus)
			self.fqs1_minus    = self.template_file.Get('fqs1_minus'); 	   self.template_list.append(self.fqs1_minus)
			self.fqs2_minus    = self.template_file.Get('fqs2_minus'); 	   self.template_list.append(self.fqs2_minus)

			self.fqa0_minus    = self.template_file.Get('fqa0_minus'); 	   self.template_list.append(self.fqa0_minus)
			self.fqa1_minus    = self.template_file.Get('fqa1_minus'); 	   self.template_list.append(self.fqa1_minus)
			self.fqa2_minus    = self.template_file.Get('fqa2_minus'); 	   self.template_list.append(self.fqa2_minus)

			self.fbck_minus    = self.template_file.Get('fbck_minus'); 	   self.template_list.append(self.fbck_minus)
			self.fntmj_minus   = self.template_file.Get('fntmj_minus');    self.template_list.append(self.fntmj_minus)

			self.F_g1 = self.fg1_plus.Integral()+self.fg1_minus.Integral()
			self.F_g2 = self.fg2_plus.Integral()+self.fg2_minus.Integral() 
			self.F_g3 = self.fg3_plus.Integral()+self.fg3_minus.Integral() 
			self.F_g4 = self.fg4_plus.Integral()+self.fg4_minus.Integral() 

			self.F_qs1 = self.fqs1_plus.Integral()+self.fqs1_minus.Integral()
			self.F_qs2 = self.fqs2_plus.Integral()+self.fqs2_minus.Integral()
		#open the total data ttree to use outside of the class
		self.data_tree = self.template_file.Get('data_tree')
#		print 'self.data_tree = '+str(self.data_tree)+', tree = '+str(tree) #DEBUG
		self.__loadDataTreeBranches__()
		#fit parameters and errors
		self.par_names = []; self.par_names_tex = []; self.pars = []; self.par_errs = []
		self.par_ini_vals = []; self.par_mins = []; self.par_maxs = []; self.par_fixs = []
		self.newParam('Rbck',  '$R_{\mathrm{bk}}$',   INITIAL_RBCK,  0.5,INITIAL_RBCK,  0.0, 1.0,False)
		self.newParam('Rntmj', '$R_{\mathrm{NTMJ}}$', INITIAL_RNTMJ, 0.5,INITIAL_RNTMJ, 0.0, 1.0,False)
		self.newParam('Rqqbar','$R_q\\bar q}$', 	  INITIAL_RQQBAR,0.5,INITIAL_RQQBAR,0.0, 1.0,False)
		self.newParam('Afb',   '$A^{(1)}_{FB}$', 	  INITIAL_AFB,   0.5,INITIAL_AFB, 	-1.0,1.0,False)
		self.newParam('d', 	   '$d$', 				  INITIAL_D, 	 0.5,INITIAL_D, 	-1.0,1.0,True)
		self.newParam('mu',    '$\\mu$', 			  INITIAL_MU, 	 0.5,INITIAL_MU, 	-1.0,1.0,True)
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