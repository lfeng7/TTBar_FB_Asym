from ROOT import *
import glob
from array import array

#global variables
#histogram limits
XBINS = 20;		XMIN = -1.0;	XMAX = 1.0
YBINS = 6;		YMIN = 0.0;		YMAX = 0.6
ZBINS = 10;		ZMIN = 500.;	ZMAX = 2500.
#luminosity
LUMINOSITY = 19748.

##############################		   Template File Class  		##############################

class template_file :
	#docstring
	"""template_file class; holds all the details of the final template file for use with theta and the auxilliary file for checking"""
	
	#__init__function
	def __init__(self,filename,sumCharges,leptons) :
		#set the input variables
		self.f = TFile(filename,'recreate')
		self.f_aux = TFile(filename.rstrip('.root')+'_aux.root')
		self.sum_charge = sumCharges == 'yes'
		self.leptype = 'none'
		lepprefix = ''
		if 'mu' in leptons :
			self.leptype = 'muons'
		if 'el' in leptons :
			self.leptype = 'electrons'
		#final distributions
		self.dists = []
		self.__addAllDistributions__()

	#addTemplate function adds a new set of histograms for the given file, and sums the new template into the appropriate distribution
	def addToTemplate(self,ttree_dir_path,file_type,file_ifd) :
		#make a new 3D histogram and projections for this sample file
		new_3D_histo = TH3D(file_type,file_type+' distribution; c*; x_{F}; M (GeV)', 	   XBINS,XMIN,XMAX,YBINS,YMIN,YMAX,ZBINS,ZMIN,ZMAX)
		new_histo_x  = TH1D(file_type+'_x',file_type+' distribution X Projection; c*', 	   XBINS,XMIN,XMAX)
		new_histo_y  = TH1D(file_type+'_y',file_type+' distribution Y Projection; x_{F}',  YBINS,YMIN,YMAX)
		new_histo_z  = TH1D(file_type+'_z',file_type+' distribution Z Projection; M (GeV)',ZBINS,ZMIN,ZMAX)
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
		if self.leptype=='electrons' :
			cuts+=ele_hadronic_pretag
		tree = chain.CopyTree(cuts)
		tree.SetDirectory(0)
		#define where to put all the branch variables we need
		self.__initializeBranchesToRead__(te)
		#loop over all events in the tree
		nEntries = tree.GetEntriesFast()
		print '	# of entries: '+str(nEntries)
		for entry in range(nEntries) :
			percent_done = 100.*entry/nEntries
			if percent_done%10 < 100./nEntries :
				print '	'+str((int)(percent_done))+'%'
			tree.GetEntry(entry)
			#calculate the reweighting factor
			eventweight = LUMINOSITY*self.weight[0]*self.sf_pileup[0]*self.sf_top_pT[0] #this will eventually be much more complicated
			NTMJ_weight = -1.0*eventweight
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
				self.NTMJ_w[0]=NTMJ_weight
				self.NTMJ_cstar[0]=self.cstar[0]
				self.NTMJ_x_F[0]=self.x_F[0] 
				self.NTMJ_M[0]=self.M[0] 
				self.NTMJ_hadt_M[0]=self.hadt_M[0] 
				self.NTMJ_Q_l[0]=self.Q_l[0]
				self.NTMJ_tree.Fill() 
			#If it's a data file and a signal event, add it to the final data tree
			if (template_ifd == 'mudata' and self.leptype=='muons') or (template_ifd == 'eledata' and self.leptype=='electrons') :
				if self.hadt_tau32[0]<0.55 and self.hadt_M[0]>140. and self.hadt_M[0]<250. :
					self.data_chi2[0] = self.chi2[0]
					self.data_cstar[0] = self.cstar[0]
					self.data_x_F[0] = self.x_F[0]
					self.data_M[0] = self.M[0]
					self.data_Q_l[0] = self.Q_l[0]
					self.data_tree.Fill()
			#add to the appropriate distributions
			for i in range(len(self.all_histos)/4) :
				dist_name = self.all_histo_names[4*i]
				#First add to the NTMJ background templates because they're sort of complex and require fewer cuts overall
				if ( 'fntmj' in dist_name and (('pass' in dist_name and self.hadt_tau32[0]<0.55) or ('fail' in dist_name and self.hadt_tau32[0]>0.55))
					and (('low' in dist_name and self.hadt_M[0]<140.) or ('hi' in dist_name and self.hadt_M[0]>250.)) ) :
					if self.sum_charge :
						self.__FillNTMJ__(i,self.cstar[0],self.x_F[0],self.M[0],self.hadt_M[0],NTMJ_weight)
					else :
						if (self.Q_l[0]>0 and 'plus' in dist_name) or (self.Q_l[0]<0 and 'minus' in dist_name) :
							self.__FillNTMJ__(i,self.cstar[0],self.x_F[0],self.M[0],self.hadt_M[0],NTMJ_weight)
				#continue if the event isn't a signal event from this point forward
				if self.hadt_tau32[0]>0.55 or self.hadt_M[0]<140. or self.hadt_M[0]>250. :
					continue
				#if it had a symmetric initial state, half the event weight and add it twice
				if self.addTwice[0] == 1 :
					eventweight = eventweight/2.
				#The background case is relatively simple
				if ((template_ifd == 'background' or template_ifd == 'bck' or template_ifd=='bkg' or template_ifd=='bg') and 'fbck' in dist_name) :
					if self.sum_charge :
						self.__Fill__(i,self.cstar[0],self.x_F[0],self.M[0],eventweight)
					else :
						if (self.Q_l[0]>0 and 'plus' in dist_name) or (self.Q_l[0]<0 and 'minus' in dist_name) :
							self.__Fill__(i,self.cstar[0],self.x_F[0],self.M[0],eventweight)
				#The qqbar and gg templates are reweighted in a bunch of different ways (and also automatically entered twice)
				elif ((template_ifd == 'qq' or template_ifd == 'qqbar') and 'fq' in dist_name) or (template_ifd == 'gg' and 'fg' in dist_name) :
					thisweight = eventweight
					thisweight_opp = eventweight
					if 'fg1' in dist_name :
						thisweight*=self.wg1[0]; thisweight_opp*=self.wg1_opp[0]
					elif 'fg2' in dist_name :
						thisweight*=self.wg2[0]; thisweight_opp*=self.wg2_opp[0]
					elif 'fg3' in dist_name :
						thisweight*=self.wg3[0]; thisweight_opp*=self.wg3_opp[0]
					elif 'fg4' in dist_name :
						thisweight*=self.wg4[0]; thisweight_opp*=self.wg4_opp[0]
					elif 'fqs1' in dist_name :
						thisweight*=self.wqs1[0]; thisweight_opp*=self.wqs1_opp[0]
					elif 'fqs2' in dist_name :
						thisweight*=self.wqs2[0]; thisweight_opp*=self.wqs2_opp[0]
					elif 'fqa0' in dist_name :
						thisweight*=self.wqa0[0]; thisweight_opp*=self.wqa0_opp[0]
					elif 'fqa1' in dist_name :
						thisweight*=self.wqa1[0]; thisweight_opp*=self.wqa1_opp[0]
					elif 'fqa2' in dist_name :
						thisweight*=self.wqa2[0]; thisweight_opp*=self.wqa2_opp[0]
					if self.sum_charge or (self.Q_l[0]>0 and 'plus' in dist_name) or (self.Q_l[0]<0 and 'minus' in dist_name) :
						self.__Fill__(i,self.cstar[0],self.x_F[0],self.M[0],thisweight)
					if self.addTwice[0] == 1 and (self.sum_charge or (self.Q_l[0]>0 and 'minus' in dist_name) or (self.Q_l[0]<0 and 'plus' in dist_name)) :
						#print 'adding an event twice to template %s with weight %f. Template is for distribution %s.'%(dist_name,thisweight_opp,template_ifd)
						self.__Fill__(i,-1.0*self.cstar[0],self.x_F[0],self.M[0],thisweight_opp)
		#write each of the new histograms to the total template file
		self.f.cd(); new_3D_histo.Write(); new_histo_x.Write(); new_histo_y.Write(); new_histo_z.Write()

	def build_NTMJ_template(self) :
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
					low_passed_i = 4*i
				elif 'low' in dist_name and 'fail' in dist_name :
					low_failed_i = 4*i
				elif 'hi' in dist_name and 'pass' in dist_name :
					hi_passed_i = 4*i
				elif 'hi' in dist_name and 'fail' in dist_name :
					hi_failed_i = 4*i
			all_low = self.all_histos[low_passed_i+3]+self.all_histos[low_failed_i+3]
			all_hi  = self.all_histos[hi_passed_i+3]+self.all_histos[hi_failed_i+3]
			conv_factor_low_x = all_low.GetMean(); conv_factor_hi_x = all_hi.GetMean()
			print 'mean X value for low mass bin = %f +/- %f'%(conv_factor_low_x,all_low.GetMeanError())
			print 'mean X value for high mass bin = %f +/- %f'%(conv_factor_hi_x,all_hi.GetMeanError())
			conv_factor_low = self.all_histos[low_passed_i].Integral()/self.all_histos[low_failed_i].Integral()
			conv_factor_low_err = conv_factor_low*sqrt(1./abs(self.all_histos[low_passed_i].Integral())+1./abs(self.all_histos[low_failed_i].Integral()))
			conv_factor_hi = self.all_histos[hi_passed_i].Integral()/self.all_histos[hi_failed_i].Integral()
			conv_factor_hi_err = conv_factor_hi*sqrt(1./abs(self.all_histos[hi_passed_i].Integral())+1./abs(self.all_histos[hi_failed_i].Integral()))
			print 'conversion factor for low mass bin = %f +/- %f'%(conv_factor_low,conv_factor_low_err)
			print 'conversion factor for high mass bin = %f +/- %f'%(conv_factor_hi,conv_factor_hi_err)
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
					low_passed_i_plus = 4*i
				elif 'low' in dist_name and 'fail' in dist_name and 'plus' in dist_name :
					low_failed_i_plus = 4*i
				elif 'hi' in dist_name and 'pass' in dist_name and 'plus' in dist_name :
					hi_passed_i_plus = 4*i
				elif 'hi' in dist_name and 'fail' in dist_name and 'plus' in dist_name :
					hi_failed_i_plus = 4*i
				elif 'low' in dist_name and 'pass' in dist_name and 'minus' in dist_name :
					low_passed_i_minus = 4*i
				elif 'low' in dist_name and 'fail' in dist_name and 'minus' in dist_name :
					low_failed_i_minus = 4*i
				elif 'hi' in dist_name and 'pass' in dist_name and 'minus' in dist_name :
					hi_passed_i_minus = 4*i
				elif 'hi' in dist_name and 'fail' in dist_name and 'minus' in dist_name :
					hi_failed_i_minus = 4*i
			all_low_plus = self.all_histos[low_passed_i_plus+3]+self.all_histos[low_failed_i_plus+3]
			all_hi_plus  = self.all_histos[hi_passed_i_plus+3]+self.all_histos[hi_failed_i_plus+3]
			all_low_minus = self.all_histos[low_passed_i_minus+3]+self.all_histos[low_failed_i_minus+3]
			all_hi_minus  = self.all_histos[hi_passed_i_minus+3]+self.all_histos[hi_failed_i_minus+3]
			conv_factor_low_x_plus = all_low_plus.GetMean(); conv_factor_hi_x_plus = all_hi_plus.GetMean()
			conv_factor_low_x_minus = all_low_minus.GetMean(); conv_factor_hi_x_minus = all_hi_minus.GetMean()
			print 'mean X value for low mass bin (positive leptons) = %f +/- %f'%(conv_factor_low_x_plus,all_low_plus.GetMeanError())
			print 'mean X value for high mass bin (positive leptons) = %f +/- %f'%(conv_factor_hi_x_plus,all_hi_plus.GetMeanError())
			print 'mean X value for low mass bin (negative leptons) = %f +/- %f'%(conv_factor_low_x_minus,all_low_minus.GetMeanError())
			print 'mean X value for high mass bin (negative leptons) = %f +/- %f'%(conv_factor_hi_x_minus,all_hi_minus.GetMeanError())
			conv_factor_low_plus = self.all_histos[low_passed_i_plus].Integral()/self.all_histos[low_failed_i_plus].Integral()
			conv_factor_low_plus_err = conv_factor_low_plus*sqrt(1./abs(self.all_histos[low_passed_i_plus].Integral())+1./abs(self.all_histos[low_failed_i_plus].Integral()))
			conv_factor_hi_plus = self.all_histos[hi_passed_i_plus].Integral()/self.all_histos[hi_failed_i_plus].Integral()
			conv_factor_hi_plus_err = conv_factor_hi_plus*sqrt(1./abs(self.all_histos[hi_passed_i_plus].Integral())+1./abs(self.all_histos[hi_failed_i_plus].Integral()))
			conv_factor_low_minus = self.all_histos[low_passed_i_minus].Integral()/self.all_histos[low_failed_i_minus].Integral()
			conv_factor_low_minus_err = conv_factor_low_minus*sqrt(1./abs(self.all_histos[low_passed_i_minus].Integral())+1./abs(self.all_histos[low_failed_i_minus].Integral()))
			conv_factor_hi_minus = self.all_histos[hi_passed_i_minus].Integral()/self.all_histos[hi_failed_i_minus].Integral()
			conv_factor_hi_minus_err = conv_factor_hi_minus*sqrt(1./abs(self.all_histos[hi_passed_i_minus].Integral())+1./abs(self.all_histos[hi_failed_i_minus].Integral()))
			print 'conversion factor for low mass bin (positive leptons) = %f +/- %f'%(conv_factor_low_plus,conv_factor_low_plus_err)
			print 'conversion factor for high mass bin (positive leptons) = %f +/- %f'%(conv_factor_hi_plus,conv_factor_hi_plus_err)
			print 'conversion factor for low mass bin (negative leptons) = %f +/- %f'%(conv_factor_low_minus,conv_factor_low_minus_err)
			print 'conversion factor for high mass bin (negative leptons) = %f +/- %f'%(conv_factor_hi_minus,conv_factor_hi_minus_err)
			slope_plus = (conv_factor_hi_plus-conv_factor_low_plus)/(conv_factor_hi_x_plus-conv_factor_low_x_plus)
			intercept_plus = conv_factor_low_plus-(slope_plus*conv_factor_low_x_plus)
			slope_minus = (conv_factor_hi_minus-conv_factor_low_minus)/(conv_factor_hi_x_minus-conv_factor_low_x_minus)
			intercept_minus = conv_factor_low_minus-(slope_minus*conv_factor_low_x_minus)
			print 'linear function (positive leptons): y = %f*x+%f'%(slope_plus,intercept_plus)
			print 'linear function (negative leptons): y = %f*x+%f'%(slope_minus,intercept_minus)
			function_plus.SetParameter(0,slope_plus); function_plus.SetParameter(1,intercept_plus)
			function_minus.SetParameter(0,slope_minus); function_minus.SetParameter(1,intercept_minus)
		#find the distributions to fill
		plus_index = 0; minus_index = 0;
		for i in range(len(self.all_histos)/4) :
			dist_name = self.all_histo_names[4*i]
			if 'fntmj' in dist_name and ('plus' in dist_name or self.sum_charge) and not 'pass' in dist_name and not 'fail' in dist_name :
				plus_index = i
			if 'fntmj' in dist_name and ('minus' in dist_name or self.sum_charge) and not 'pass' in dist_name and not 'fail' in dist_name :
				minus_index = i
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
				#new_weight = w[0]
				new_weight = function_plus.Eval(t_M[0])*w[0]
				self.__Fill__(plus_index,c[0],x[0],M[0],new_weight)
			elif Q[0]<0 :
				#new_weight = w[0]
				new_weight = function_minus.Eval(t_M[0])*w[0]
				self.__Fill__(minus_index,c[0],x[0],M[0],new_weight)

	#__addAllDistributions__ sets up all of the final distributions depending on whether we want the charges summed
	def __addAllDistributions__(self) :
		lepprefix = 'none'
		if self.leptype == 'muons' :
			lepprefix = 'mu'
		elif self.leptype == 'electrons' :
			lepprefix = 'el'
		std_reweights = ['weight','sf_top_pT','sf_pileup']
		std_systematics = ['sf_lep_id']
		if self.sum_charge :
			(self,name,formatted_name,sample_reweight,reweights,systematics,function)
			self.dists.append(distribution(lepprefix+'__fg0','0th gg (qg,q_{i}q_{j},etc.) distribution',None,std_reweights,std_systematics,'(3.-#Rbck#-#Rntmj#)*(2.-#Rqqbar#)'))
			self.dists.append(distribution(lepprefix+'__fqs0','0th Symmetric q#bar{q distribution',None,std_reweights,std_systematics,'(3.-#Rbck#-#Rntmj#)*#Rqqbar#'))
			self.dists.append(distribution(lepprefix+'__fqa0','0th Antisymmetric q#bar{q distribution','wqa0',std_reweights,std_systematics,'(3.-#Rbck#-#Rntmj#)*#Rqqbar#*#Afb#'))
			self.dists.append(distribution(lepprefix+'__fbck','background distribution',None,std_reweights,std_systematics,'#Rbck#'))
		else :
			print "You're SOL because I haven't written this yet"

	#__addDistribution__ function adds a 3D histogram and 1D projections to the file and to the lists
	def __addDistribution__(self,name,formatted_name) :
		histo_3D = TH3D(name,	   formatted_name+'; c*; x_{F}; M (GeV)',XBINS,XMIN,XMAX,YBINS,YMIN,YMAX,ZBINS,ZMIN,ZMAX)
		histo_x  = TH1D(name+'_x',formatted_name+' X Projection; c*',	  XBINS,XMIN,XMAX)
		histo_y  = TH1D(name+'_y',formatted_name+' Y Projection; x_{F}',		  YBINS,YMIN,YMAX)
		histo_z  = TH1D(name+'_z',formatted_name+' Z Projection; M (GeV)',		  ZBINS,ZMIN,ZMAX)
		self.all_histo_names.append(name); self.all_histo_names.append(name+'_x') 
		self.all_histo_names.append(name+'_y'); self.all_histo_names.append(name+'_z')
		self.all_histos.append(histo_3D); self.all_histos.append(histo_x); self.all_histos.append(histo_y); self.all_histos.append(histo_z)
		histo_3D.SetDirectory(0); histo_x.SetDirectory(0); histo_y.SetDirectory(0); histo_z.SetDirectory(0)

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
		histo_3D.SetDirectory(0); histo_x.SetDirectory(0); histo_y.SetDirectory(0); histo_z.SetDirectory(0)

	#__initializeBranchesToRead__ function just puts all the variables from the inputted TTree into variables we can use
	def __initializeBranchesToRead__(self,tree) :
		#hadronic top jet mass and tau32
		self.hadt_M = array('d',[-1.0]); tree.SetBranchAddress('hadt_M',self.hadt_M)
		self.hadt_tau32 = array('d',[-1.0]); tree.SetBranchAddress('hadt_tau32',self.hadt_tau32)
		#weights (renormalization, scale factors, analysis)
		self.weight    = array('d',[1.0]); tree.SetBranchAddress('weight',self.weight)
		self.wg1 	   	   = array('d',[1.0]); tree.SetBranchAddress('wg1',		  self.wg1)
		self.wg2 	   	   = array('d',[1.0]); tree.SetBranchAddress('wg2',		  self.wg2)
		self.wg3 	   	   = array('d',[1.0]); tree.SetBranchAddress('wg3',		  self.wg3)
		self.wg4 	   	   = array('d',[1.0]); tree.SetBranchAddress('wg4',		  self.wg4)
		self.wqs1 	   	   = array('d',[1.0]); tree.SetBranchAddress('wqs1', 	  self.wqs1)
		self.wqs2 	   	   = array('d',[1.0]); tree.SetBranchAddress('wqs2', 	  self.wqs2)
		self.wqa0 	   	   = array('d',[1.0]); tree.SetBranchAddress('wqa0', 	  self.wqa0)
		self.wqa1 	   	   = array('d',[1.0]); tree.SetBranchAddress('wqa1', 	  self.wqa1)
		self.wqa2 	   	   = array('d',[1.0]); tree.SetBranchAddress('wqa2', 	  self.wqa2)
		self.wg1_opp 	   = array('d',[1.0]); tree.SetBranchAddress('wg1_opp',   self.wg1_opp)
		self.wg2_opp 	   = array('d',[1.0]); tree.SetBranchAddress('wg2_opp',   self.wg2_opp)
		self.wg3_opp 	   = array('d',[1.0]); tree.SetBranchAddress('wg3_opp',   self.wg3_opp)
		self.wg4_opp 	   = array('d',[1.0]); tree.SetBranchAddress('wg4_opp',   self.wg4_opp)
		self.wqs1_opp 	   = array('d',[1.0]); tree.SetBranchAddress('wqs1_opp',  self.wqs1_opp)
		self.wqs2_opp 	   = array('d',[1.0]); tree.SetBranchAddress('wqs2_opp',  self.wqs2_opp)
		self.wqa0_opp 	   = array('d',[1.0]); tree.SetBranchAddress('wqa0_opp',  self.wqa0_opp)
		self.wqa1_opp 	   = array('d',[1.0]); tree.SetBranchAddress('wqa1_opp',  self.wqa1_opp)
		self.wqa2_opp 	   = array('d',[1.0]); tree.SetBranchAddress('wqa2_opp',  self.wqa2_opp)
		self.sf_pileup 	   = array('d',[1.0]); tree.SetBranchAddress('sf_pileup', self.sf_pileup)
		self.sf_top_pT 	   = array('d',[1.0]); tree.SetBranchAddress('sf_top_pT', self.sf_top_pT)
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

	#__FillNTMJ__ function just fills the 3D histo, some 1-D projections, and the hadronic top mass distribution
	def __FillNTMJ__(self,i,c,x,m,tm,w) :
		self.all_histos[4*i].Fill(c,x,m,w)
		self.all_histos[4*i+1].Fill(c,w)
		self.all_histos[4*i+2].Fill(x,w)
		self.all_histos[4*i+3].Fill(tm,w)

	#__del__ function
	def __del__(self) :
		#write to the file
		self.f.cd()
		for histo in self.all_histos :
			histo.Write()
		self.data_tree.Write()
		self.f.Write()
		self.f.Close()

class distribution :
	#docstring
	"""distribution class"""
	
	#__init__function
	def __init__(self,name,formatted_name,sample_reweight,reweights,systematics,function) :
		self.name = name
		self.formatted_name = formatted_name
		self.function = function
		self.tree = TTree(self.name+'_tree',self.name+'_tree')
		self.tree.SetDirectory(0)
		self.all_branches = []
		self.cstar 			  = array('d',[100.]); self.tree.Branch('cstar',self.cstar,'cstar/D'); 								 self.all_branches.append(('cstar_scaled',self.cstar))
		self.x_F 			  = array('d',[100.]); self.tree.Branch('x_F',self.x_F,'x_F/D'); 									 self.all_branches.append(('x_F_scaled',self.x_F))
		self.M  			  = array('d',[-1.0]); self.tree.Branch('M',self.M,'M/D'); 											 self.all_branches.append(('M_scaled',self.M))
		self.hadt_M  		  = array('d',[-1.0]); self.tree.Branch('hadt_M',self.hadt_M,'hadt_M/D'); 							 self.all_branches.append(('scaled_hadt_M',self.hadt_M))
		self.hadt_tau32  	  = array('d',[-1.0]); self.tree.Branch('hadt_tau32',self.hadt_tau32,'hadt_tau32/D'); 				 self.all_branches.append(('hadt_tau32',self.hadt_tau32))
		self.Q_l  			  = array('i',[0]);    self.tree.Branch('Q_l',self.Q_l,'Q_l/I'); 									 self.all_branches.append(('Q_l',self.Q_l))
		self.addTwice 		  = array('I',[2]);    self.tree.Branch('addTwice',self.addTwice,'addTwice/i'); 					 self.all_branches.append(('addTwice',self.addTwice))
		self.sample_reweight  = array('d',[1.0]);  self.tree.Branch('sample_reweight',self.sample_reweight,'sample_reweight/D'); self.all_branches((sample_reweight,self.sample_reweight))
		self.reweights_arrays = []
		for i in range(len(reweights)) :
			self.reweights_arrays.append(array('d',[1.0]))
			self.tree.Branch(reweights[i],self.reweights_arrays[i],reweights[i]+'/D')
			self.all_branches.append((reweights[i],self.reweights_arrays[i]))
		self.systematics_arrays 	 = []
		self.systematics_arrays_up   = []
		self.systematics_arrays_down = []
		for i in range(len(systematics)) :
			self.systematics_arrays.append(array('d',[1.0]))
			self.tree.Branch(systematics[i],self.systematics_arrays[len(self.systematics_arrays)-1],systematics[i]+'/D')
			self.all_branches.append((systematics[i],self.systematics_arrays[len(self.systematics_arrays)-1]))
			self.systematics_arrays_up.append(array('d',[1.0]))
			self.tree.Branch(systematics[i]+'_up',self.systematics_arrays[len(self.systematics_arrays)-1],systematics[i]+'_up/D')
			self.all_branches.append((systematics[i]+'_up',self.systematics_arrays[len(self.systematics_arrays)-1]))
			self.systematics_arrays_down.append(array('d',[1.0]))
			self.tree.Branch(systematics[i]+'_down',self.systematics_arrays[len(self.systematics_arrays)-1],systematics[i]+'_down/D')
			self.all_branches.append((systematics[i]+'_down',self.systematics_arrays[len(self.systematics_arrays)-1]))

	#__del__ function
	def __del__(self) :
		print 'lol'
