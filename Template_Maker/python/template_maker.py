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
	"""template_file class"""
	
	#__init__function
	def __init__(self,filename,sumCharges,leptons) :
		#set the input variables
		self.f = TFile(filename,'recreate')
		self.f_aux = TFile(filename.rstrip('.root')+'_aux.root','recreate')
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
		#histograms
		self.histos = []

	#addTemplate function adds a new set of histograms for the given file, and sums the new template into the appropriate distribution
	def addToDistributions(self,ttree_dir_path,file_type,file_ifd) :
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
		#find the distributions this tree will contribute to
		for i in range(len(self.dists)) :
			distname = self.dists[i].name
			#if this tree should be added to this distribution
			contribution = self.__contributesToDist__(file_ifd,i)
			if contribution != 0. :
				#set all of the branch addresses
				for branch in self.dists[i].all_branches :
					tree.SetBranchAddress(branch[0],branch[2])
				#copy over the relevant parts of the tree
				nEntries = tree.GetEntriesFast()
				print 'Adding trees from '+ttree_dir_path+' to distribution '+distname+' ('+str(nEntries)+' entries)'
				for entry in range(nEntries) :
					percent_done = 100.*entry/nEntries
					if percent_done%10 < 100./nEntries :
						print '	'+str((int)(percent_done))+'%'
					tree.GetEntry(entry)
					#change the event weight to add or subtract as appropriate
					self.dists[i].contrib_weight[0]=contribution
					self.dists[i].tree.Fill()

	#build_templates builds the template histograms from the distributions
	def build_templates(self) :
		for dist in self.dists :
			if not 'fntmj' in dist.name :
				print 'Building templates for distribution '+dist.name+'. . .'
				self.histos+=dist.build_templates(self.f_aux,self.sum_charge)
				print 'Done'

	#build_NTMJ_templates automatically generates templates for all of the data-driven NTMJ background
	#automatically accounting for systematics
	def build_NTMJ_templates(self) :
		ntmjdists = []
		for dist in self.dists :
			if 'ntmj' in dist.name :
				ntmjdists.append(dist)
		for ntmjdist in ntmjdists :
			#first make room for all of the templates we'll eventually have
			#Nominal and due to fitting function parameters
			ntmjdist.__organizeFunctionTemplates__()
			#Due to systematics
			ntmjdist.__organizeSystematicsTemplates__()
			#Break the MC-subtracted data events into regions
			antitag_tree  = ntmjdist.tree.CopyTree('hadt_M > 140. && hadt_M < 250. && hadt_tau32 > 0.55')
			sideband_trees = []; sideband_counts = []; sideband_xs = []
			low_pass_tree = ntmjdist.tree.CopyTree('hadt_M < 140. && hadt_tau32 < 0.55'); sideband_trees.append(low_pass_tree)
			low_fail_tree = ntmjdist.tree.CopyTree('hadt_M < 140. && hadt_tau32 > 0.55'); sideband_trees.append(low_fail_tree)
			hi_pass_tree  = ntmjdist.tree.CopyTree('hadt_M > 250. && hadt_tau32 < 0.55'); sideband_trees.append(hi_pass_tree)
			hi_fail_tree  = ntmjdist.tree.CopyTree('hadt_M > 250. && hadt_tau32 > 0.55'); sideband_trees.append(hi_fail_tree)
			#for each of the sideband regions get the total event count
			for i in range(len(sideband_trees)) :
				sideband_counts.append([]); sideband_xs.append([])
				#get ready to read all of the branches
				for branch in ntmjdist.all_branches :
					sideband_trees[i].SetBranchAddress(branch[1],branch[2])
				#Find the overall number of events in the region in each morphing case
				nEntries = sideband_trees[i].GetEntriesFast()
				for entry in range(nEntries) :
					sideband_trees[i].GetEntry(entry)
					#for each template, find the total number of events in this sideband
					for j in range(len(ntmjdist.all_histos)/4) :
						sideband_counts[i].append(0.); sideband_xs[i].append(0.)
					for j in range(len(ntmjdist.all_histos)/4) :
						template_name = ntmjdist.all_histo_names[4*i]
						#Build the reweighting factors
						eventweight = ntmjdist.contrib_weight[0]*ntmjdist.sample_reweight[0]
						eventweight_opp = ntmjdist.contrib_weight[0]*ntmjdist.sample_reweight_opp[0]
						for reweight in ntmjdist.reweights_arrays :
							eventweight*=reweight[0]; eventweight_opp*=reweight[0]
						#apply all the necessary systematic factors
						if ntmjdist.systematics != None :
							for k in range(len(ntmjdist.systematics)) :
								if ntmjdist.systematics[k]+'__down' in template_name :
									eventweight*=ntmjdist.systematics_arrays_down[k][0]; eventweight_opp*=ntmjdist.systematics_arrays_down[k][0]
								elif ntmjdist.systematics[k]+'__up' in template_name : 
									eventweight*=ntmjdist.systematics_arrays_up[k][0]; eventweight_opp*=ntmjdist.systematics_arrays_up[k][0]
								else :
									eventweight*=ntmjdist.systematics_arrays[k][0]; eventweight_opp*=ntmjdist.systematics_arrays[k][0]
						#apply the function of the fitting parameters
						eventweight*=ntmjdist.factors[j]; eventweight_opp*=ntmjdist.factors[j]
						#if we want to add the event twice half the weight
						if ntmjdist.addTwice[0]==1 :
							eventweight*=0.5; eventweight_opp*=0.5
						#add to the running total
						if self.sum_charge or (ntmjdist.Q_l[0]>0 and 'plus' in template_name) or (ntmjdist.Q_l[0]<0 and 'minus' in template_name) :
							sideband_counts[i][j]+=eventweight; sideband_xs[i][j]+=eventweight*ntmjdist.hadt_M[0]
						if ntmjdist.addTwice[0]==1 and (self.sum_charge or (ntmjdist.Q_l[0]<0 and 'plus' in template_name) or (ntmjdist.Q_l[0]>0 and 'minus' in template_name)) :
							sideband_counts[i][j]+=eventweight_opp; sideband_xs[i][j]+=eventweight_opp*ntmjdist.hadt_M[0]
			#Calculate the conversion factor functions for each template
			conversion_functions = []
			for i in range(len(ntmjdist.all_histos)/4) :
				conversion_functions.append(TF1('conv_func_'+ntmjdist.all_histo_names[4*i],'[0]*x+[1]',100.,500.))
				x_low = (sideband_xs[0][i]+sideband_xs[1][i])/(sideband_counts[0][i]+sideband_counts[1][i])
				x_hi  = (sideband_xs[2][i]+sideband_xs[3][i])/(sideband_counts[2][i]+sideband_counts[3][i])
				y_low = sideband_counts[0][i]/sideband_counts[1][i]
				y_hi  = sideband_counts[2][i]/sideband_counts[3][i]
				slope = (y_hi-y_low)/(x_hi-x_low)
				intercept = y_low-(slope*x_low)
				conversion_functions[i].SetParameter(0,slope); conversion_functions[i].SetParameter(1,intercept)
			#Now build the final NTMJ templates by morphing the antitagged regions into the signal regions
			for branch in ntmjdist.all_branches :
				antitag_tree.SetBranchAddress(branch[1],branch[2])
			nEntries = antitag_tree.GetEntriesFast()
			for entry in range(nEntries) :
				antitag_tree.GetEntry(entry)
				#add to all of the templates
				for j in range(len(ntmjdist.all_histos)/4) :
					template_name = ntmjdist.all_histo_names[4*i]
					#Build the reweighting factors
					eventweight = ntmjdist.contrib_weight[0]*ntmjdist.sample_reweight[0]
					eventweight_opp = ntmjdist.contrib_weight[0]*ntmjdist.sample_reweight_opp[0]
					for reweight in ntmjdist.reweights_arrays :
						eventweight*=reweight[0]; eventweight_opp*=reweight[0]
					#apply all the necessary systematic factors
					if ntmjdist.systematics != None :
						for k in range(len(ntmjdist.systematics)) :
							if ntmjdist.systematics[k]+'__down' in template_name :
								eventweight*=ntmjdist.systematics_arrays_down[k][0]; eventweight_opp*=ntmjdist.systematics_arrays_down[k][0]
							elif ntmjdist.systematics[k]+'__up' in template_name : 
								eventweight*=ntmjdist.systematics_arrays_up[k][0]; eventweight_opp*=ntmjdist.systematics_arrays_up[k][0]
							else :
								eventweight*=ntmjdist.systematics_arrays[k][0]; eventweight_opp*=ntmjdist.systematics_arrays[k][0]
					#apply the function of the fitting parameters
					eventweight*=ntmjdist.factors[j]; eventweight_opp*=ntmjdist.factors[j]
					#apply the conversion function from antitagged to signal region
					eventweight*=conversion_functions[j].Eval(ntmjdist.hadt_M[0]); eventweight_opp*=conversion_functions[j].Eval(ntmjdist.hadt_M[0])
					#if we want to add the event twice half the weight
					if ntmjdist.addTwice[0]==1 :
						eventweight*=0.5; eventweight_opp*=0.5
					#add to the running total
					if self.sum_charge or (ntmjdist.Q_l[0]>0 and 'plus' in template_name) or (ntmjdist.Q_l[0]<0 and 'minus' in template_name) :
						ntmjdist.__Fill__(j,ntmjdist.cstar[0],ntmjdist.x_F[0],ntmjdist.M[0],eventweight)
					if ntmjdist.addTwice[0]==1 and (self.sum_charge or (ntmjdist.Q_l[0]<0 and 'plus' in template_name) or (ntmjdist.Q_l[0]>0 and 'minus' in template_name)) :
						ntmjdist.__Fill__(j,-1.0*ntmjdist.cstar[0],ntmjdist.x_F[0],ntmjdist.M[0],eventweight_opp)
			#Save the templates to the auxiliary file
			self.f_aux.cd()
			returnhistos = []
			for histo in ntmjdist.all_histos :
				histo.Write()
			for i in range(len(ntmjdist.all_histos)/4) :
				returnhistos.append(convertTo1D(ntmjdist.all_histos[4*i]))
			#add to the list of new 1D templates
			self.histos+=returnhistos

	#__addAllDistributions__ sets up all of the final distributions depending on whether we want the charges summed
	def __addAllDistributions__(self) :
		lepprefix = 'none'
		if self.leptype == 'muons' :
			lepprefix = 'mu'
		elif self.leptype == 'electrons' :
			lepprefix = 'el'
		std_reweights = ['weight','sf_top_pT','sf_pileup']
		std_systematics = ['sf_lep_ID']
		if self.sum_charge :
			self.dists.append(distribution(lepprefix+'__DATA','data distribution',None,None,None,None))
			self.dists.append(distribution(lepprefix+'__fg0','0th gg (qg,q_{i}q_{j},etc.) distribution',None,std_reweights,std_systematics,'(3.-#Rbck#-#Rntmj#)*(2.-#Rqqbar#)'))
			self.dists.append(distribution(lepprefix+'__fqs0','0th Symmetric q#bar{q} distribution',None,std_reweights,std_systematics,'(3.-#Rbck#-#Rntmj#)*#Rqqbar#'))
			self.dists.append(distribution(lepprefix+'__fqa0','0th Antisymmetric q#bar{q} distribution','wqa0',std_reweights,std_systematics,'(3.-#Rbck#-#Rntmj#)*#Rqqbar#*#Afb#'))
			self.dists.append(distribution(lepprefix+'__fbck','background distribution',None,std_reweights,std_systematics,'#Rbck#'))
			self.dists.append(distribution(lepprefix+'__fntmj','NTMJ background distribution',None,None,std_systematics,'#Rntmj#'))
		else :
			print "You're SOL because I haven't written this yet"

	#__contributesToDist__ finds out whether the tree coming in should be added to the ith distribution
	def __contributesToDist__(self,ifd,i) :
		distname = self.dists[i].name
		if (self.leptype=='muons' and ifd=='mudata' and ('DATA' in distname or 'fntmj' in distname)) :
			return 1.0
		elif (self.leptype=='electrons' and ifd=='eledata' and ('DATA' in distname or 'fntmj' in distname)) :
			return 1.0
		elif (ifd == 'qq' and ('fqs' in distname or 'fqa' in distname)) :
			return LUMINOSITY
		elif (ifd == 'gg' and 'fg' in distname) :
			return LUMINOSITY
		elif (ifd == 'bck' and 'fbck' in distname) :
			return LUMINOSITY
		elif (ifd == 'qq' or ifd == 'gg' or ifd == 'bck') and 'fntmj' in distname :
			return -1.0*LUMINOSITY
		else :
			return 0.0

	#__del__ function
	def __del__(self) :
		#write to the file
		self.f_aux.cd()
		for dist in self.dists :
			dist.tree.Write()
		self.f_aux.Close()
		self.f.cd()
		for histo in self.histos :
			histo.Write()
		self.f.Close()

class distribution :
	#docstring
	"""distribution class"""
	
	#__init__function
	def __init__(self,name,formatted_name,sample_reweight,reweights,systematics,function) :
		self.name = name
		self.formatted_name = formatted_name
		self.systematics = systematics
		self.function = function
		self.tree = TTree(self.name+'_tree',self.name+'_tree')
		self.tree.SetDirectory(0)
		self.all_branches = []
		self.cstar 			  = array('d',[100.]); self.tree.Branch('cstar',self.cstar,'cstar/D'); 								 self.all_branches.append(('cstar_scaled','cstar',self.cstar))
		self.x_F 			  = array('d',[100.]); self.tree.Branch('x_F',self.x_F,'x_F/D'); 									 self.all_branches.append(('x_F_scaled','x_F',self.x_F))
		self.M  			  = array('d',[-1.0]); self.tree.Branch('M',self.M,'M/D'); 											 self.all_branches.append(('M_scaled','M',self.M))
		self.hadt_M  		  = array('d',[-1.0]); self.tree.Branch('hadt_M',self.hadt_M,'hadt_M/D'); 							 self.all_branches.append(('scaled_hadt_M','hadt_M',self.hadt_M))
		self.hadt_tau32  	  = array('d',[-1.0]); self.tree.Branch('hadt_tau32',self.hadt_tau32,'hadt_tau32/D'); 				 self.all_branches.append(('hadt_tau32','hadt_tau32',self.hadt_tau32))
		self.Q_l  			  = array('i',[0]);    self.tree.Branch('Q_l',self.Q_l,'Q_l/I'); 									 self.all_branches.append(('Q_l','Q_l',self.Q_l))
		self.addTwice 		  = array('I',[2]);    self.tree.Branch('addTwice',self.addTwice,'addTwice/i'); 					 self.all_branches.append(('addTwice','addTwice',self.addTwice))
		self.contrib_weight   = array('d',[1.0]);  self.tree.Branch('contrib_weight',self.contrib_weight,'contrib_weight/D');
		self.sample_reweight  = array('d',[1.0]);  self.tree.Branch('sample_reweight',self.sample_reweight,'sample_reweight/D');
		self.sample_reweight_opp  = array('d',[1.0]);  self.tree.Branch('sample_reweight_opp',self.sample_reweight_opp,'sample_reweight_opp/D'); 
		if sample_reweight != None :
			self.all_branches.append((sample_reweight,'sample_reweight',self.sample_reweight))
			self.all_branches.append((sample_reweight+'_opp','sample_reweight_opp',self.sample_reweight_opp))
		self.reweights_arrays = []
		if reweights != None :
			for i in range(len(reweights)) :
				self.reweights_arrays.append(array('d',[1.0]))
				self.tree.Branch(reweights[i],self.reweights_arrays[i],reweights[i]+'/D')
				self.all_branches.append((reweights[i],reweights[i],self.reweights_arrays[i]))
		self.systematics_arrays 	 = []
		self.systematics_arrays_up   = []
		self.systematics_arrays_down = []
		if systematics != None :
			for i in range(len(systematics)) :
				self.systematics_arrays.append(array('d',[1.0]))
				self.tree.Branch(systematics[i],self.systematics_arrays[len(self.systematics_arrays)-1],systematics[i]+'/D')
				self.all_branches.append((systematics[i],systematics[i],self.systematics_arrays[len(self.systematics_arrays)-1]))
				self.systematics_arrays_up.append(array('d',[1.0]))
				self.tree.Branch(systematics[i]+'_hi',self.systematics_arrays_up[len(self.systematics_arrays_up)-1],systematics[i]+'_hi/D')
				self.all_branches.append((systematics[i]+'_hi',systematics[i]+'_hi',self.systematics_arrays_up[len(self.systematics_arrays_up)-1]))
				self.systematics_arrays_down.append(array('d',[1.0]))
				self.tree.Branch(systematics[i]+'_low',self.systematics_arrays_down[len(self.systematics_arrays_down)-1],systematics[i]+'_low/D')
				self.all_branches.append((systematics[i]+'_low',systematics[i]+'_low',self.systematics_arrays_down[len(self.systematics_arrays_down)-1]))
		self.all_histo_names = []; self.all_histos = []; self.all_templates = []

	#build_templates builds the 3D templates for a distribution, and returns a list of 1D templates
	def build_templates(self,aux_file,sumcharge) :
		#Nominal and due to fitting function parameters
		self.__organizeFunctionTemplates__()
		#Due to systematics
		self.__organizeSystematicsTemplates__()
		#look only at the signal events
		signal_tree = self.tree.CopyTree('hadt_M > 140. && hadt_M < 250. && hadt_tau32 < 0.55')
		#set the branches to read
		for branch in self.all_branches :
			signal_tree.SetBranchAddress(branch[1],branch[2])
		#read and add events from the tree to all of the derived templates
		nEntries = signal_tree.GetEntriesFast()
		for entry in range(nEntries) :
			signal_tree.GetEntry(entry)
			for i in range(len(self.all_histos)/4) :
				#the common reweighting factors
				eventweight = self.contrib_weight[0]*self.sample_reweight[0]
				eventweight_opp = self.contrib_weight[0]*self.sample_reweight_opp[0]
				for reweight in self.reweights_arrays :
					eventweight*=reweight[0]; eventweight_opp*=reweight[0]
				template_name = self.all_histo_names[4*i]
				#apply all the necessary systematic factors
				if self.systematics != None :
					for j in range(len(self.systematics)) :
						if self.systematics[j]+'__down' in template_name :
							eventweight*=self.systematics_arrays_down[j][0]; eventweight_opp*=self.systematics_arrays_down[j][0]
						elif self.systematics[j]+'__up' in template_name : 
							eventweight*=self.systematics_arrays_up[j][0]; eventweight_opp*=self.systematics_arrays_up[j][0]
						else :
							eventweight*=self.systematics_arrays[j][0]; eventweight_opp*=self.systematics_arrays[j][0]
				#apply the function of the fitting parameters
				eventweight*=self.factors[i]; eventweight_opp*=self.factors[i]
				#if we want to add the event twice half the weight
				if self.addTwice[0]==1 :
					eventweight*=0.5; eventweight_opp*=0.5
				#add to the template
				if sumcharge or (self.Q_l[0]>0 and 'plus' in template_name) or (self.Q_l[0]<0 and 'minus' in template_name) :
					self.__Fill__(i,self.cstar[0],self.x_F[0],self.M[0],eventweight)
				if self.addTwice[0]==1 and (sumcharge or (self.Q_l[0]<0 and 'plus' in template_name) or (self.Q_l[0]>0 and 'minus' in template_name)) :
					self.__Fill__(i,-1.0*self.cstar[0],self.x_F[0],self.M[0],eventweight_opp)
		#Save the templates to the auxiliary file
		aux_file.cd()
		returnhistos = []
		for histo in self.all_histos :
			histo.Write()
		for i in range(len(self.all_histos)/4) :
			returnhistos.append(convertTo1D(self.all_histos[4*i]))
		#return the list of new 1D templates
		return returnhistos

	#__organizeFunctionTemplates__ adds all of the necessary up/down templates for the fitting parameters
	#and returns the nominal factor for scaling the templates that are morphed for systematics
	def __organizeFunctionTemplates__(self) :
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
		#lists of new template prefactors and names
		self.factors = []; args = []; argispar = []
		if self.function != None :
			args = self.function.split('#')
		for j in range(len(args)) :
			argispar.append(False)
		for j in range(len(args)) : #replace the parameters within the arguments
			for par in pars :
				if par[0] == args[j] :
					args[j] = par
					argispar[j] = True
		#nominal distribution
		factorstring = ''
		for j in range(len(args)) :
			if not argispar[j] :
				factorstring+=args[j]
			else :
				factorstring+=str(args[j][1])
		if factorstring != '' :
			self.factors.append(eval(factorstring))
		else :
			self.factors.append(1.0)
		self.__addTemplate__(self.name,self.formatted_name)
		print 'adding template with name %s and factor %s=%.4f'%(self.name,factorstring,self.factors[len(self.factors)-1])
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
					self.__addTemplate__(self.name+'__'+args[j][0]+'__'+'down',self.formatted_name+' '+args[j][0]+' down')
					self.__addTemplate__(self.name+'__'+args[j][0]+'__'+'up',self.formatted_name+' '+args[j][0]+' up')
					print 'adding template with name %s and factor %s=%.4f'%(self.name+'__'+args[j][0]+'__'+'down',factorstring_down,self.factors[len(self.factors)-2])
					print 'adding template with name %s and factor %s=%.4f'%(self.name+'__'+args[j][0]+'__'+'up',factorstring_up,  self.factors[len(self.factors)-1])
					countpars+=1
			self.factors.append(eval(factorstring_down)); self.factors.append(eval(factorstring_up))

	#__organizeSystematicsTemplates__ adds all of the templates for the up/down systematics distributions
	def __organizeSystematicsTemplates__(self) :
		if self.systematics == None :
			return
		for sys in self.systematics :
			self.factors.append(1.0)
			self.__addTemplate__(self.name+'__'+sys+'__down',self.formatted_name+' '+sys+' down')
			self.factors.append(1.0)
			self.__addTemplate__(self.name+'__'+sys+'__up',self.formatted_name+' '+sys+' up')
			print 'adding template with name '+self.name+'__'+sys+'__down'
			print 'adding template with name '+self.name+'__'+sys+'__up'

	#__addDistribution__ function adds a 3D histogram and 1D projections to the file and to the lists
	def __addTemplate__(self,name,formatted_name) :
		histo_3D = TH3D(name,	   formatted_name+'; c*; |x_{F}|; M (GeV)',XBINS,XMIN,XMAX,YBINS,YMIN,YMAX,ZBINS,ZMIN,ZMAX)
		histo_x  = TH1D(name+'_x',formatted_name+' X Projection; c*',	  XBINS,XMIN,XMAX)
		histo_y  = TH1D(name+'_y',formatted_name+' Y Projection; |x_{F}|',		  YBINS,YMIN,YMAX)
		histo_z  = TH1D(name+'_z',formatted_name+' Z Projection; M (GeV)',		  ZBINS,ZMIN,ZMAX)
		self.all_histo_names.append(name); self.all_histo_names.append(name+'_x') 
		self.all_histo_names.append(name+'_y'); self.all_histo_names.append(name+'_z')
		self.all_histos.append(histo_3D); self.all_histos.append(histo_x); self.all_histos.append(histo_y); self.all_histos.append(histo_z)
		histo_3D.SetDirectory(0); histo_x.SetDirectory(0); histo_y.SetDirectory(0); histo_z.SetDirectory(0)

	#__Fill__ function just fills the 3D histo and its 1D projections
	def __Fill__(self,i,c,x,m,w) :
		self.all_histos[4*i].Fill(c,x,m,w)
		self.all_histos[4*i+1].Fill(c,w)
		self.all_histos[4*i+2].Fill(x,w)
		self.all_histos[4*i+3].Fill(m,w)

	#__del__ function
	def __del__(self) :
		print 'lol'

#convertTo1D takes a 3D distribution and makes it 1D for use with theta
def convertTo1D(original) :
	nBins = original.GetNbinsX()*original.GetNbinsY()*original.GetNbinsZ()
	newHisto = TH1F(original.GetName(),original.GetTitle(),nBins,0.,nBins-1.)
	newHisto.SetDirectory(0)
	for k in range(nBins) :
		newHisto.SetBinContent(k,original.GetBinContent(k))
	return newHisto