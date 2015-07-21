#template_maker workhorse code used in making templates from ttrees with a bunch of options
#contains three classes, one for background templates and one for signal templates, and one 
#for the total template file
#signal templates produce many histograms at once because they have antisymmetric, xi, and
#delta reweighting factors to keep track of
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

from ROOT import *
import glob
from array import array

#global variables
#histogram limits
XBINS = 20;		XMIN = -1.0;	XMAX = 1.0
YBINS = 20;		YMIN = 0.0;		YMAX = 0.6
ZBINS = 20;		ZMIN = 500.;	ZMAX = 2500.
#luminosity
LUMINOSITY = 19748.

##############################		   Template File Class  		##############################

class template_file :
	#docstring
	"""template_file class; holds all the details of the final template file with the summed distributions"""
	
	#__init__function
	def __init__(self,filename,sumCharges,leptons) :
		#set the input variables
		self.f = TFile(filename,'recreate')
		self.sum_charge = sumCharges == 'yes'
		self.leptype = 'none'
		if 'mu' in leptons :
			self.leptype = 'muons'
		if 'el' in leptons :
			self.leptype = 'electrons'
		#final distributions
		self.all_histo_names = []; self.all_histos = []
		if self.sum_charge :
			self.__addDistribution__('fqqs',	  'Symmetric q#bar{q} distribution')
			self.__addDistribution__('fqqs_xi',	  'Symmetric q#bar{q} distribution, xi reweighted')
			self.__addDistribution__('fqqs_delta','Symmetric q#bar{q} distribution, delta reweighted')
			self.__addDistribution__('fqqa',	  'Antisymmetric q#bar{q} distribution')
			self.__addDistribution__('fqqa_xi',	  'Antisymmetric q#bar{q} distribution, xi reweighted')
			self.__addDistribution__('fqqa_delta','Antisymmetric q#bar{q} distribution, delta reweighted')
			self.__addDistribution__('fgg',	  	  'gg (qg,q_{i}q_{j},etc.) distribution')
			self.__addDistribution__('fbck',	  'background distribution')
			self.__addDistribution__('fntmj',	  'NTMJ background distribution')
			self.__addNTMJDistribution__('fntmj_low_pass',	  'NTMJ background distribution, low mass, passing #tau_{32} cut')
			self.__addNTMJDistribution__('fntmj_low_fail',	  'NTMJ background distribution, low mass, failing #tau_{32} cut')
			self.__addNTMJDistribution__('fntmj_hi_pass',	  'NTMJ background distribution, high mass, passing #tau_{32} cut')
			self.__addNTMJDistribution__('fntmj_hi_fail',	  'NTMJ background distribution, high mass, failing #tau_{32} cut')
		else :
			self.__addDistribution__('fqqs_plus',	   'Symmetric q#bar{q} distribution, Q_{l}>0')
			self.__addDistribution__('fqqs_xi_plus',   'Symmetric q#bar{q} distribution, xi reweighted, Q_{l}>0')
			self.__addDistribution__('fqqs_delta_plus','Symmetric q#bar{q} distribution, delta reweighted, Q_{l}>0')
			self.__addDistribution__('fqqa_plus',	   'Antisymmetric q#bar{q} distribution, Q_{l}>0')
			self.__addDistribution__('fqqa_xi_plus',   'Antisymmetric q#bar{q} distribution, xi reweighted, Q_{l}>0')
			self.__addDistribution__('fqqa_delta_plus','Antisymmetric q#bar{q} distribution, delta reweighted, Q_{l}>0')
			self.__addDistribution__('fgg_plus',	   'gg (qg,q_{i}q_{j},etc.) distribution, Q_{l}>0')
			self.__addDistribution__('fbck_plus',	   'background distribution, Q_{l}>0')
			self.__addDistribution__('fntmj_plus',	   'NTMJ background distribution, Q_{l}>0')
			self.__addNTMJDistribution__('fntmj_plus_low_pass',	  'NTMJ background distribution, low mass, passing #tau_{32} cut, Q_{l}>0')
			self.__addNTMJDistribution__('fntmj_plus_low_fail',	  'NTMJ background distribution, low mass, failing #tau_{32} cut, Q_{l}>0')
			self.__addNTMJDistribution__('fntmj_plus_hi_pass',	  'NTMJ background distribution, high mass, passing #tau_{32} cut, Q_{l}>0')
			self.__addNTMJDistribution__('fntmj_plus_hi_fail',	  'NTMJ background distribution, high mass, failing #tau_{32} cut, Q_{l}>0')
			self.__addDistribution__('fqqs_minus',	    'Symmetric q#bar{q} distribution, Q_{l}<0')
			self.__addDistribution__('fqqs_xi_minus',   'Symmetric q#bar{q} distribution, xi reweighted, Q_{l}<0')
			self.__addDistribution__('fqqs_delta_minus','Symmetric q#bar{q} distribution, delta reweighted, Q_{l}<0')
			self.__addDistribution__('fqqa_minus',	    'Antisymmetric q#bar{q} distribution, Q_{l}<0')
			self.__addDistribution__('fqqa_xi_minus',   'Antisymmetric q#bar{q} distribution, xi reweighted, Q_{l}<0')
			self.__addDistribution__('fqqa_delta_minus','Antisymmetric q#bar{q} distribution, delta reweighted, Q_{l}<0')
			self.__addDistribution__('fgg_minus',	    'gg (qg,q_{i}q_{j},etc.) distribution, Q_{l}<0')
			self.__addDistribution__('fbck_minus',	    'background distribution, Q_{l}<0')
			self.__addDistribution__('fntmj_minus',	    'NTMJ background distribution, Q_{l}<0')
			self.__addNTMJDistribution__('fntmj_minus_low_pass',	  'NTMJ background distribution, low mass, passing #tau_{32} cut, Q_{l}<0')
			self.__addNTMJDistribution__('fntmj_minus_low_fail',	  'NTMJ background distribution, low mass, failing #tau_{32} cut, Q_{l}<0')
			self.__addNTMJDistribution__('fntmj_minus_hi_pass',	  'NTMJ background distribution, high mass, passing #tau_{32} cut, Q_{l}<0')
			self.__addNTMJDistribution__('fntmj_minus_hi_fail',	  'NTMJ background distribution, high mass, failing #tau_{32} cut, Q_{l}<0')
		self.NTMJ_tree = TTree('NTMJ_tree','NTMJ_tree')
		self.NTMJ_weight = array('d',[1.0]); self.NTMJ_tree.Branch('NTMJ_weight',self.NTMJ_weight,'NTMJ_weight/D')
		self.NTMJ_cstar  = array('d',[100.]); self.NTMJ_tree.Branch('NTMJ_cstar',self.NTMJ_cstar,'NTMJ_cstar/D')
		self.NTMJ_x_F  = array('d',[100.]); self.NTMJ_tree.Branch('NTMJ_x_F',self.NTMJ_x_F,'NTMJ_x_F/D')
		self.NTMJ_M  = array('d',[-1.0]); self.NTMJ_tree.Branch('NTMJ_M',self.NTMJ_M,'NTMJ_M/D')
		self.NTMJ_hadt_M  = array('d',[-1.0]); self.NTMJ_tree.Branch('NTMJ_hadt_M',self.NTMJ_hadt_M,'NTMJ_hadt_M/D')
		self.NTMJ_Q_l  = array('i',[0]); self.NTMJ_tree.Branch('NTMJ_Q_l',self.NTMJ_Q_l,'NTMJ_Q_l/I')

	#addTemplate function adds a new set of histograms for the given file, and sums the new template into the appropriate distribution
	def addToTemplate(self,ttree_dir_path,template_name,template_ifd) :
		#make a new 3D histogram and projections for this sample file
		new_3D_histo = TH3D(template_name,template_name+' distribution; c*; x_{F}; M (GeV)',XBINS,XMIN,XMAX,YBINS,YMIN,YMAX,ZBINS,ZMIN,ZMAX)
		new_histo_x  = TH1D(template_name+'_x',template_name+' distribution X Projection; c*',	  XBINS,XMIN,XMAX)
		new_histo_y  = TH1D(template_name+'_y',template_name+' distribution Y Projection; x_{F}',			  YBINS,YMIN,YMAX)
		new_histo_z  = TH1D(template_name+'_z',template_name+' distribution Z Projection; M (GeV)',		  ZBINS,ZMIN,ZMAX)
		#chain up the ttree files
		chain = TChain('tree')
		filenamelist = glob.glob(ttree_dir_path+'/*_skim_tree.root')
		for filename in filenamelist :
			chain.Add(filename)
		#copy the chain into a tree after making selection cuts
		muon_preselection = 'muon1_pt>ele1_pt && lepW_pt>50. && hadt_pt>300. && hadt_M>100.'
		muon_kinematics = 'muon1_pt>40. && abs(muon1_eta)<2.4'
		muon_ID = 'muon1_isLoose==1'
		muon_2D = '(muon1_relPt>25. || muon1_dR>0.5)'
		ele_preselection = 'ele1_pt>muon1_pt && lepW_pt>50. && hadt_pt>300. && hadt_M>100.'
		ele_kinematics = 'ele1_pt>40. && abs(ele1_eta)<2.4'
		ele_ID = 'ele1_isLoose==1'
		ele_2D = '(ele1_relPt>25. || ele1_dR>0.5)'
		lep_top_mass = 'lept_M>140. && lept_M<250.'
		muon_full_leptonic = muon_preselection+' && '+muon_kinematics+' && '+muon_ID+' && '+muon_2D+' && '+lep_top_mass
		muon_hadronic_pretag = muon_full_leptonic+' && hadt_tau21>0.1'
		ele_full_leptonic = ele_preselection+' && '+ele_kinematics+' && '+ele_ID+' && '+ele_2D+' && '+lep_top_mass
		ele_hadronic_pretag = ele_full_leptonic+' && hadt_tau21>0.1'
		signal_mass = 'hadt_M>140. && hadt_M<250.'
		signal_tau32 = 'hadt_tau32<0.55' 
		cuts = ''
		if self.leptype=='muons' :
			cuts+=muon_hadronic_pretag
			#cuts+=' && (('+muon_hadronic_pretag+') || ('+ele_hadronic_pretag+'))'
		if self.leptype=='electrons' :
			cuts+=ele_hadronic_pretag
		tree = chain.CopyTree(cuts)
		#define where to put all the branch variables we need
		self.__initializeBranchesToRead__(tree)
		#loop over all events in the tree
		nEntries = tree.GetEntriesFast()
		print '	# of entries: '+str(nEntries)
		for entry in range(nEntries) :
			percent_done = 100.*entry/nEntries
			if percent_done%10 < 100./nEntries :
				print '	'+str((int)(percent_done))+'%'
			tree.GetEntry(entry)
			#calculate the reweighting factor
			eventweight = self.weight[0]*self.sf_pileup[0] #this will eventually be much more complicated
			NTMJ_weight = -1.0*LUMINOSITY*eventweight
			if (template_ifd == 'mudata' and self.leptype=='muons') or (template_ifd == 'eledata' and self.leptype=='electrons') :
				NTMJ_weight = 1.0
			if (template_ifd == 'eledata' and self.leptype=='muons') or (template_ifd == 'mudata' and self.leptype=='electrons') :
				NTMJ_weight = 0.0
			#add to the new histograms (only once, these are just for double-checks)
			new_3D_histo.Fill(self.cstar[0],self.x_F[0],self.M[0],eventweight)
			new_histo_x.Fill(self.cstar[0],eventweight)
			new_histo_y.Fill(self.x_F[0],eventweight)
			new_histo_z.Fill(self.M[0],eventweight)
			#add antitagged events to the NTMJ tree
			if self.hadt_tau32[0]>0.55 and self.hadt_M[0]>140. and self.hadt_M[0]<250. :
				self.NTMJ_weight[0]=NTMJ_weight
				self.NTMJ_cstar[0]=self.cstar[0]
				self.NTMJ_x_F[0]=self.x_F[0] 
				self.NTMJ_M[0]=self.M[0] 
				self.NTMJ_hadt_M[0]=self.hadt_M[0] 
				self.NTMJ_Q_l[0]=self.Q_l[0]
				self.NTMJ_tree.Fill() 
			#add to the appropriate distributions
			for i in range(len(self.all_histos)/4) :
				dist_name = self.all_histo_names[4*i]
				#First add to the NTMJ background templates because they're sort of complex and require fewer cuts overall
				if ( 'fntmj' in dist_name and (('pass' in dist_name and self.hadt_tau32[0]<0.55) or ('pass' not in dist_name and self.hadt_tau32[0]>0.55))
					and (('low' in dist_name and self.hadt_M[0]<140.) or ('hi' in dist_name and self.hadt_M[0]>250.)) ) :
					if self.sum_charge :
						self.__Fill__(i,self.cstar[0],self.x_F[0],self.hadt_M[0],NTMJ_weight)
					else :
						if (self.Q_l[0]>0 and 'plus' in dist_name) or (self.Q_l[0]<0 and 'minus' in dist_name) :
							self.__Fill__(i,self.cstar[0],self.x_F[0],self.hadt_M[0],NTMJ_weight)
				#continue if the event isn't a signal event from this point forward
				if self.hadt_tau32[0]>0.55 or self.hadt_M[0]<140. or self.hadt_M[0]>250. :
					continue
				#if it had a symmetric initial state, half the event weight and add it twice
				if self.addTwice[0] == 1 :
					eventweight = eventweight/2.
				#The background and gg cases are relatively simple
				if ( ((template_ifd == 'background' or template_ifd == 'bck' or template_ifd=='bkg' or template_ifd=='bg') and 'fbck' in dist_name)
					 or (template_ifd == 'gg' and 'fgg' in dist_name) ) :
					if self.sum_charge :
						self.__Fill__(i,self.cstar[0],self.x_F[0],self.M[0],eventweight)
						if self.addTwice[0] == 1 :
							self.__Fill__(i,-1.0*self.cstar[0],self.x_F[0],self.M[0],eventweight)
					else :
						if (self.Q_l[0]>0 and 'plus' in dist_name) or (self.Q_l[0]<0 and 'minus' in dist_name) :
							self.__Fill__(i,self.cstar[0],self.x_F[0],self.M[0],eventweight)
						if self.addTwice[0] == 1 and ((self.Q_l[0]>0 and 'minus' in dist_name) or (self.Q_l[0]<0 and 'plus' in dist_name)) :
							self.__Fill__(i,-1.0*self.cstar[0],self.x_F[0],self.M[0],eventweight)
				#The qqbar templates are reweighted in a bunch of different ways (and also automatically entered twice)
				elif template_ifd == 'qq' or template_ifd == 'qqbar' :
					if self.sum_charge :
						if 'fqqs' in dist_name and not ('xi' in dist_name or 'delta' in dist_name) :
							self.__Fill__(i,self.cstar[0],self.x_F[0],self.M[0],eventweight)
							self.__Fill__(i,-1.0*self.cstar[0],self.x_F[0],self.M[0],eventweight)
						elif 'fqqs_xi' in dist_name :
							self.__Fill__(i,self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_s_xi[0])
							self.__Fill__(i,-1.0*self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_s_xi_opp[0])
						elif 'fqqs_delta' in dist_name :
							self.__Fill__(i,self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_s_delta[0])
							self.__Fill__(i,-1.0*self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_s_delta_opp[0])
						elif 'fqqa' in dist_name and not ('xi' in dist_name or 'delta' in dist_name) :
							self.__Fill__(i,self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_a[0])
							self.__Fill__(i,-1.0*self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_a_opp[0])
						elif 'fqqa_xi' in dist_name :
							self.__Fill__(i,self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_a_xi[0])
							self.__Fill__(i,-1.0*self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_a_xi_opp[0])
						elif 'fqqa_delta' in dist_name :
							self.__Fill__(i,self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_a_delta[0])
							self.__Fill__(i,-1.0*self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_a_delta_opp[0])
					else :
						if (self.Q_l[0]>0 and 'plus' in dist_name) or (self.Q_l[0]<0 and 'minus' in dist_name) :
							if 'fqqs' in dist_name and not ('xi' in dist_name or 'delta' in dist_name) :
								self.__Fill__(i,self.cstar[0],self.x_F[0],self.M[0],eventweight)
							elif 'fqqs_xi' in dist_name :
								self.__Fill__(i,self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_s_xi[0])
							elif 'fqqs_delta' in dist_name :
								self.__Fill__(i,self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_s_delta[0])
							elif 'fqqa' in dist_name and not ('xi' in dist_name or 'delta' in dist_name) :
								self.__Fill__(i,self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_a[0])
							elif 'fqqa_xi' in dist_name :
								self.__Fill__(i,self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_a_xi[0])
							elif 'fqqa_delta' in dist_name :
								self.__Fill__(i,self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_a_delta[0])
						if (self.Q_l[0]>0 and 'minus' in dist_name) or (self.Q_l[0]<0 and 'plus' in dist_name) :
							if 'fqqs' in dist_name and not ('xi' in dist_name or 'delta' in dist_name) :
								self.__Fill__(i,-1.0*self.cstar[0],self.x_F[0],self.M[0],eventweight)
							elif 'fqqs_xi' in dist_name :
								self.__Fill__(i,-1.0*self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_s_xi_opp[0])
							elif 'fqqs_delta' in dist_name :
								self.__Fill__(i,-1.0*self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_s_delta_opp[0])
							elif 'fqqa' in dist_name and not ('xi' in dist_name or 'delta' in dist_name) :
								self.__Fill__(i,-1.0*self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_a_opp[0])
							elif 'fqqa_xi' in dist_name :
								self.__Fill__(i,-1.0*self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_a_xi_opp[0])
							elif 'fqqa_delta' in dist_name :
								self.__Fill__(i,-1.0*self.cstar[0],self.x_F[0],self.M[0],eventweight*self.w_a_delta_opp[0])
		#write each of the new histograms to the total template file
		self.f.cd(); new_3D_histo.Write(); new_histo_x.Write(); new_histo_y.Write(); new_histo_z.Write()

	#normalizeDistributions function used to renormalize the final distributions in preparation for fitting
	def normalizeDistributions(self) :
		#three separate normalizations factors for the simpler distributions
		qq_int  = 0.;	gg_int  = 0.;	bck_int = 0.
		#for each distribution, add its integral to the appropriate factor
		for i in range(len(self.all_histos)/4) :
			dist_name = self.all_histo_names[4*i]
			if 'fqqs' in dist_name and not 'xi' in dist_name and not 'delta' in dist_name :
				qq_int+=self.all_histos[4*i].Integral()
			elif 'fgg' in dist_name :
				gg_int+=self.all_histos[4*i].Integral()
			elif 'fbck' in dist_name :
				bck_int+=self.all_histos[4*i].Integral()
		#renormalize all the distributions, including their projections, accordingly
		for i in range(len(self.all_histos)/4) :
			dist_name = self.all_histo_names[4*i]
			factor = 1.
			if 'fqqs' in dist_name or 'fqqa' in dist_name :
				factor = 1.0/qq_int
			elif 'fgg' in dist_name :
				factor = 1.0/gg_int
			elif 'fbck' in dist_name :
				factor = 1.0/bck_int
			self.all_histos[4*i+0].Scale(factor)
			self.all_histos[4*i+1].Scale(factor)
			self.all_histos[4*i+2].Scale(factor)
			self.all_histos[4*i+3].Scale(factor)
		#conversion factors for both mass sidebands and linear function
		function_plus = TF1('function_plus','[0]*x+[1]',100.,500.)
		function_minus = TF1('function_minus','[0]*x+[1]',100.,500.)
		if self.sum_charge :
			low_passed_i = 0; low_failed_i = 0; hi_passed_i = 0; hi_failed_i = 0
			for i in range(len(self.all_histos)/4) :
				dist_name = self.all_histo_names[4*i]
				if not 'fntmj' in dist_name :
					continue
				if 'low' in dist_name and 'pass' in dist_name :
					low_passed_i = 4*i+3
				elif 'low' in dist_name and 'fail' in dist_name :
					low_failed_i = 4*i+3
				elif 'hi' in dist_name and 'pass' in dist_name :
					hi_passed_i = 4*i+3
				elif 'hi' in dist_name and 'fail' in dist_name :
					hi_failed_i = 4*i+3
			all_low = all_histos[low_passed_i]+all_histos[low_failed_i]
			all_hi  = all_histos[hi_passed_i]+all_histos[hi_failed_i]
			conv_factor_low_x = all_low.GetMean(); conv_factor_hi_x = all_hi.GetMean()
			print 'mean X value for low mass bin = %f'%(conv_factor_low_x)
			print 'mean X value for high mass bin = %f'%(conv_factor_hi_x)
			conv_factor_low = all_histos[low_passed_i].Integral()/all_histos[low_failed_i].Integral()
			conv_factor_hi = all_histos[hi_passed_i].Integral()/all_histos[hi_failed_i].Integral()
			slope = (conv_factor_hi-conv_factor_low)/(conv_factor_hi_x-conv_factor_low_x)
			intercept = conv_factor_low-(slope*conv_factor_low_x)
			print 'linear function: y = %f*x+%f'%(slope,intercept)
			function_plus.SetParameter(0,slope); function_plus.SetParameter(1,intercept)
			function_minus.SetParameter(0,slope); function_minus.SetParameter(1,intercept)
		else :
			low_passed_i_plus = 0; low_failed_i_plus = 0; hi_passed_i_plus = 0; hi_failed_i_plus = 0
			low_passed_i_minus = 0; low_failed_i_minus = 0; hi_passed_i_minus = 0; hi_failed_i_minus = 0
			for i in range(len(self.all_histos)/4) :
				dist_name = self.all_histo_names[4*i]
				if not 'fntmj' in dist_name :
					continue
				if 'low' in dist_name and 'pass' in dist_name and 'plus' in dist_name :
					low_passed_i_plus = 4*i+3
				elif 'low' in dist_name and 'fail' in dist_name and 'plus' in dist_name :
					low_failed_i_plus = 4*i+3
				elif 'hi' in dist_name and 'pass' in dist_name and 'plus' in dist_name :
					hi_passed_i_plus = 4*i+3
				elif 'hi' in dist_name and 'fail' in dist_name and 'plus' in dist_name :
					hi_failed_i_plus = 4*i+3
				elif 'low' in dist_name and 'pass' in dist_name and 'minus' in dist_name :
					low_passed_i_minus = 4*i+3
				elif 'low' in dist_name and 'fail' in dist_name and 'minus' in dist_name :
					low_failed_i_minus = 4*i+3
				elif 'hi' in dist_name and 'pass' in dist_name and 'minus' in dist_name :
					hi_passed_i_minus = 4*i+3
				elif 'hi' in dist_name and 'fail' in dist_name and 'minus' in dist_name :
					hi_failed_i_minus = 4*i+3
			all_low_plus = all_histos[low_passed_i_plus]+all_histos[low_failed_i_plus]
			all_hi_plus  = all_histos[hi_passed_i_plus]+all_histos[hi_failed_i_plus]
			all_low_minus = all_histos[low_passed_i_minus]+all_histos[low_failed_i_minus]
			all_hi_minus  = all_histos[hi_passed_i_minus]+all_histos[hi_failed_i_minus]
			conv_factor_low_x_plus = all_low_plus.GetMean(); conv_factor_hi_x_plus = all_hi_plus.GetMean()
			conv_factor_low_x_minus = all_low_minus.GetMean(); conv_factor_hi_x_minus = all_hi_minus.GetMean()
			print 'mean X value for low mass bin (positive leptons) = %f'%(conv_factor_low_x_plus)
			print 'mean X value for high mass bin (positive leptons) = %f'%(conv_factor_hi_x_plus)
			print 'mean X value for low mass bin (negative leptons) = %f'%(conv_factor_low_x_minus)
			print 'mean X value for high mass bin (negative leptons) = %f'%(conv_factor_hi_x_minus)
			conv_factor_low_plus = all_histos[low_passed_i_plus].Integral()/all_histos[low_failed_i_plus].Integral()
			conv_factor_hi_plus = all_histos[hi_passed_i_plus].Integral()/all_histos[hi_failed_i_plus].Integral()
			conv_factor_low_minus = all_histos[low_passed_i_minus].Integral()/all_histos[low_failed_i_minus].Integral()
			conv_factor_hi_minus = all_histos[hi_passed_i_minus].Integral()/all_histos[hi_failed_i_minus].Integral()
			slope_plus = (conv_factor_hi_plus-conv_factor_low_plus)/(conv_factor_hi_x_plus-conv_factor_low_x_plus)
			intercept_plus = conv_factor_low_plus-(slope_plus*conv_factor_low_x_plus)
			slope_minus = (conv_factor_hi_minus-conv_factor_low_minus)/(conv_factor_hi_x_minus-conv_factor_low_x_minus)
			intercept_minus = conv_factor_low_minus-(slope_minus*conv_factor_low_x_minus)
			print 'linear function (positive leptons): y = %f*x+%f'%(slope_plus,intercept_plus)
			print 'linear function (negative leptons): y = %f*x+%f'%(slope_minus,intercept_minus)
			function_plus.SetParameter(0,slope_plus); function_plus.SetParameter(1,intercept_plus)
			function_minus.SetParameter(0,slope_minus); function_minus.SetParameter(1,intercept_minus)
		self.__build_NTMJ_template__(function_plus,function_minus)

	#__addDistribution__ function adds a 3D histogram and 1D projections to the file and to the lists
	def __addDistribution__(self,name,formatted_name) :
		histo_3D = TH3D(name,	   formatted_name+'; c*; x_{F}; M (GeV)',XBINS,XMIN,XMAX,YBINS,YMIN,YMAX,ZBINS,ZMIN,ZMAX)
		histo_x  = TH1D(name+'_x',formatted_name+' X Projection; c*',	  XBINS,XMIN,XMAX)
		histo_y  = TH1D(name+'_y',formatted_name+' Y Projection; x_{F}',		  YBINS,YMIN,YMAX)
		histo_z  = TH1D(name+'_z',formatted_name+' Z Projection; M (GeV)',		  ZBINS,ZMIN,ZMAX)
		self.all_histo_names.append(name); self.all_histo_names.append(name+'_x') 
		self.all_histo_names.append(name+'_y'); self.all_histo_names.append(name+'_z')
		self.all_histos.append(histo_3D); self.all_histos.append(histo_x); self.all_histos.append(histo_y); self.all_histos.append(histo_z)

	#__addNTMJDistribution__ function adds a 3D histogram and 1D projections to the file and to the lists 
	#but with hadronic t mass instead of ttbar mass on the z-axis
	def __addNTMJDistribution__(self,name,formatted_name) :
		histo_3D = TH3D(name,	   formatted_name+'; c*; x_{F}; M_{hadt} (GeV)',XBINS,XMIN,XMAX,YBINS,YMIN,YMAX,ZBINS,ZMIN,ZMAX)
		histo_x  = TH1D(name+'_x',formatted_name+' X Projection; c*',	  XBINS,XMIN,XMAX)
		histo_y  = TH1D(name+'_y',formatted_name+' Y Projection; x_{F}',		  YBINS,YMIN,YMAX)
		histo_z  = TH1D(name+'_z',formatted_name+' Z Projection; M_{hadt} (GeV)',		  40,100.,500.)
		self.all_histo_names.append(name); self.all_histo_names.append(name+'_x') 
		self.all_histo_names.append(name+'_y'); self.all_histo_names.append(name+'_z')
		self.all_histos.append(histo_3D); self.all_histos.append(histo_x); self.all_histos.append(histo_y); self.all_histos.append(histo_z)

	#__initializeBranchesToRead__ function just puts all the variables from the inputted TTree into variables we can use
	def __initializeBranchesToRead__(self,tree) :
		#hadronic top jet mass and tau32
		self.hadt_M = array('d',[-1.0]); tree.SetBranchAddress('hadt_M',self.hadt_M)
		self.hadt_tau32 = array('d',[-1.0]); tree.SetBranchAddress('hadt_tau32',self.hadt_tau32)
		#weights (renormalization, scale factors, analysis)
		self.weight    = array('d',[1.0]); tree.SetBranchAddress('weight',self.weight)
		self.w_a 	   	   = array('d',[1.0]); tree.SetBranchAddress('w_a',		  self.w_a)
		self.w_s_xi    	   = array('d',[1.0]); tree.SetBranchAddress('w_s_xi',	      self.w_s_xi)
		self.w_a_xi    	   = array('d',[1.0]); tree.SetBranchAddress('w_a_xi',	      self.w_a_xi)
		self.w_s_delta 	   = array('d',[1.0]); tree.SetBranchAddress('w_s_delta',    self.w_s_delta)
		self.w_a_delta 	   = array('d',[1.0]); tree.SetBranchAddress('w_a_delta',    self.w_a_delta)
		self.w_a_opp 	   = array('d',[1.0]); tree.SetBranchAddress('w_a_opp',	  self.w_a_opp)
		self.w_s_xi_opp    = array('d',[1.0]); tree.SetBranchAddress('w_s_xi_opp',	  self.w_s_xi_opp)
		self.w_a_xi_opp    = array('d',[1.0]); tree.SetBranchAddress('w_a_xi_opp',	  self.w_a_xi_opp)
		self.w_s_delta_opp = array('d',[1.0]); tree.SetBranchAddress('w_s_delta_opp',self.w_s_delta_opp)
		self.w_a_delta_opp = array('d',[1.0]); tree.SetBranchAddress('w_a_delta_opp',self.w_a_delta_opp)
		self.sf_pileup 	   = array('d',[1.0]); tree.SetBranchAddress('sf_pileup', self.sf_pileup)
		#lepton charge
		self.Q_l = array('i',[0]); tree.SetBranchAddress('Q_l',self.Q_l)
		#kinematic fit chi2
		self.chi2 = array('d',[0.0]); tree.SetBranchAddress('chi2',self.chi2)
		#whether or not this event should be added twice and have its weight halved based on whether its initial state
		#was symmetric (this will only be nonzero for qqbar and some gg events)
		self.addTwice = array('I',[0]); tree.SetBranchAddress('addTwice',self.addTwice)
		#cosine(theta)
		self.cstar = array('d',[100.0]); tree.SetBranchAddress('cstar_scaled',self.cstar)
		#Feynman x
		self.x_F = array('d',[100.0]); tree.SetBranchAddress('x_F_scaled',self.x_F)
		#ttbar invariant mass
		self.M = array('d',[-1.0]); tree.SetBranchAddress('M_scaled',self.M)

	#__Fill__ function just fills the 3D histo and its 1D projections
	def __Fill__(self,i,c,x,m,w) :
		self.all_histos[4*i].Fill(c,x,m,w)
		self.all_histos[4*i+1].Fill(c,w)
		self.all_histos[4*i+2].Fill(x,w)
		self.all_histos[4*i+3].Fill(m,w)

	def __build_NTMJ_template__(self,plus_func,minus_func) :
		#find the distributions to fill
		plus_index = 0; minus_index = 0;
		for i in range(len(self.all_histos)/4) :
			dist_name = self.all_histo_names[4*i]
			if 'fntmj' in dist_name and ('plus' in dist_name or self.sum_charge) :
				plus_index = 4*i
			if 'fntmj' in dist_name and ('minus' in dist_name or self.sum_charge) :
				minus_index = 4*i
		#initialize the branches to read from
		w = array('d',[1.0]); self.NTMJ_tree.SetBranchAddress('NTMJ_weight',w)
		c = array('d',[100.]); self.NTMJ_tree.SetBranchAddress('NTMJ_cstar',c)
		x = array('d',[100.]); self.NTMJ_tree.SetBranchAddress('NTMJ_x_F',x)
		M = array('d',[-1.0]); self.NTMJ_tree.SetBranchAddress('NTMJ_M',M)
		t_M = array('d',[-1.0]); self.NTMJ_tree.SetBranchAddress('NTMJ_hadt_M',t_M)
		Q = array('i',[0]); self.NTMJ_tree.SetBranchAddress('NTMJ_Q_l',Q)
		#read all the entries
		nEntries = self.NTMJ_tree.GetEntriesFast()
		print '	# of entries in NTMJ tree: '+str(nEntries)
		for entry in range(nEntries) :
			percent_done = 100.*entry/nEntries
			if percent_done%10 < 100./nEntries :
				print '	'+str((int)(percent_done))+'%'
			self.NTMJ_tree.GetEntry(entry)
			#recalculate weights and fill histograms
			if Q[0]>0 :
				new_weight = plus_func.Eval(t_M[0])*w[0]
				self.__Fill__(plus_index,c[0],x[0],M[0],new_weight)
			elif Q[0]<0 :
				new_weight = minus_func.Eval(t_M[0])*w[0]
				self.__Fill__(minus_index,c[0],x[0],M[0],new_weight)



	#__del__ function
	def __del__(self) :
		#write to the file
		self.f.cd()
		for histo in self.all_histos :
			histo.Write()
		self.f.Write()
		self.f.Close()