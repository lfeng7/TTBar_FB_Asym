from ROOT import *
from array import array
from template import template

#global variables
#list of constant reweights
const_reweights_trees = ['weight']
const_reweights_dists = ['cs_weight']
#list of systematic reweights
simple_systematics_trees = ['sf_pileup']#,     'sf_top_pT',     'sf_lep_ID',     'sf_trig_eff', 	  'luminosity']
simple_systematics_dists = ['pileup_weight']#, 'top_pT_weight', 'lep_ID_weight', 'trig_eff_weight', 'luminosity']
#PDF reweight
pdf_reweight_vector_trees = 'CT10_weights'
pdf_reweight_vector_dists = 'CT10_weights'
#Luminosity min bias
LUMI_MIN_BIAS = 0.026
#Cuts for conversion function calculation
LM1_LOW = 50.
LM1_HI  = 95.
LM2_LOW = 95.
LM2_HI  = 140.
SM_LOW  = 140.
SM_HI   = 250.
HM_LOW  = 250.
HM_HI   = 350.
SUBS_CUT = 0.55


#TDR Style
#gROOT.Macro('rootlogon.C')

##############################		   Distribution Class  		##############################

class distribution :
	#docstring
	"""distribution class"""
	
	#__init__function
	def __init__(self,parfile_name,includepdf,name,color,formatted_name,function) :
		#cosmetic stuff
		print '		Creating distribution named '+name+''
		self.parfile_name = parfile_name
		self.include_PDF = includepdf
		self.step = parfile_name.split('_')[0]
		self.name = name
		self.color = color
		self.formatted_name = formatted_name
		#fitting function string
		self.function = function
		#trees and branches
		self.signal_tree     = TTree(self.name+'_signal_tree',self.name+'_signal_tree')
		self.antitagged_tree = TTree(self.name+'_antitagged_tree',self.name+'_antitagged_tree')
		self.sidebands_tree  = TTree(self.name+'_sidebands_tree',self.name+'_sidebands_tree')
		self.signal_tree.SetDirectory(0)
		self.antitagged_tree.SetDirectory(0)
		self.sidebands_tree.SetDirectory(0)
		self.all_branches = []
		print '			Adding branches'
		self.__add_all_branches__()
		#automatically generate the list of templates that will be made from this distribution
		self.all_templates = []
		print '			Adding templates'
		self.__add_all_templates__()
		print '		Done'

	def filltree(self) :
		hadtM = self.hadt_M[0]; hadttau32 = self.hadt_tau32[0]
		if hadtM>SM_LOW and hadtM<SM_HI :
			if hadttau32<SUBS_CUT :
				self.signal_tree.Fill()
			else :
				self.antitagged_tree.Fill()
		else:
			self.sidebands_tree.Fill()

	#build_templates builds the 3D templates for a distribution, and returns a list of 1D templates
	def build_templates(self,distlist) :
		#set the branches to read
		for branch in self.all_branches :
			self.signal_tree.SetBranchAddress(branch[1],branch[2])
		nEntries = self.signal_tree.GetEntries()
		#read and add events from the tree to all of the derived templates
		for i in range(len(self.all_templates)) :
			template_name = self.all_templates[i].name
			print '		Doing template '+template_name
			#Get the total numbers of events for 
			NGG,NG1,NG2,NG3,NG4,NQQ,NQ1,NQ2,NBCK,NNTMJ = self.__get_total_event_numbers__(distlist,'sm_ps',template_name)
#			print '		NGG, NG1, NG2, NG3, NG4, NQQ, NQ1, NQ2, NBCK, NNTMJ = %f, %f, %f, %f, %f, %f, %f, %f, %f, %f'%(NGG,NG1,NG2,NG3,NG4,NQQ,NQ1,NQ2,NBCK,NNTMJ) #DEBUG
			for entry in range(nEntries) :
				self.signal_tree.GetEntry(entry)
				eventweight, eventweight_opp = self.__get_event_weights__(template_name)
				funcweight, funcweight_opp = self.__get_func_weights__(template_name,NGG,NG1,NG2,NG3,NG4,NQQ,NQ1,NQ2,NBCK,NNTMJ)
				eventweight*=funcweight; eventweight_opp*=funcweight_opp
				#if we want to add the event twice half the weight
				if self.addTwice[0]==1 :
					eventweight*=0.5; eventweight_opp*=0.5
				#add to the template
				self.all_templates[i].Fill(self.cstar[0],abs(self.x_F[0]),self.M[0],eventweight)
				if self.addTwice[0]==1 :
					self.all_templates[i].Fill(-1.0*self.cstar[0],abs(self.x_F[0]),self.M[0],eventweight_opp)

	def add_to_NTMJ_dictionary(self,distlist,event_counts,tree_cuts) :
		#Set Branches
		for branch in self.all_branches :
			self.sidebands_tree.SetBranchAddress(branch[1],branch[2])
		nEntries = self.sidebands_tree.GetEntries()
		#For all of the lepton types
		for ltype in event_counts :
			#check that this distribution contributes to templates with leptons of this type
			if self.name.split('__')[0].split('plus')[0] != ltype and self.name.split('__')[0].split('minus')[0] != ltype :
				continue
			#for all of the template types
			for ttype in event_counts[ltype] :
				if (self.name.find('JES__up')!=-1 and not ttype.find('JES__up')!=-1) or (self.name.find('JER__up')!=-1 and not ttype.find('JER__up')!=-1) :
					continue
				if (ttype.find('JES__up')!=-1 and not self.name.find('JES__up')!=-1) or (ttype.find('JER__up')!=-1 and not self.name.find('JER__up')!=-1) :
					continue
				if (self.name.find('JES__down')!=-1 and not ttype.find('JES__down')!=-1) or (self.name.find('JER__down')!=-1 and not ttype.find('JER__down')!=-1) :
					continue
				if (ttype.find('JES__down')!=-1 and not self.name.find('JES__down')!=-1) or (ttype.find('JER__down')!=-1 and not self.name.find('JER__down')!=-1) :
					continue
				template_name = ltype+'__fntmj__'+ttype
				print '			Adding with lepton type '+ltype+' and template type '+ttype #DEBUG	
				#Get the total event numbers based on the template type and set of cuts
				NGGS = []; NG1S = []; NG2S = []; NG3S = []; NG4S = []; NQQS = []; NQ1S = []; NQ2S = []; NBCKS = []; NNTMJS = []
				for j in range(len(event_counts[ltype][ttype])) :
					first, second = tree_cuts[j].split('_')[0], tree_cuts[j].split('_')[1]
					if not (first == 'lm1' or first == 'lm2' or first =='sm' or first == 'hm') or not (second == 'ps' or second == 'fs') :
						print 'WARNING: CUTSTRING NOT RECOGNIZED: '+tree_cuts[j]
						return
					NGG,NG1,NG2,NG3,NG4,NQQ,NQ1,NQ2,NBCK,NNTMJ = self.__get_total_event_numbers__(distlist,tree_cuts[j],template_name)
					NNTMJ = 0. #because we're going to find what it is
					NGGS.append(NGG); NG1S.append(NG1); NG2S.append(NG2); NG3S.append(NG3); NG4S.append(NG4) 
					NQQS.append(NQQ); NQ1S.append(NQ1); NQ2S.append(NQ2) 
					NBCKS.append(NBCK) 
					NNTMJS.append(NNTMJ)
				#some loop optimization
				getentry = self.sidebands_tree.GetEntry
				geteventweights = self.__get_event_weights__
				getfuncweights = self.__get_func_weights__
				for entry in range(nEntries) :
					getentry(entry)
					hadtM = self.hadt_M[0]; hadttau32 = self.hadt_tau32[0]
					lm1 = hadtM>LM1_LOW and hadtM<=LM1_HI
					lm2 = hadtM>LM2_LOW and hadtM<=LM2_HI
					hm  = hadtM>HM_LOW and hadtM<=HM_HI
					sm  = hadtM>SM_LOW and hadtM<=SM_HI
					ps  = hadttau32<SUBS_CUT
					fs  = hadttau32>SUBS_CUT
					eventweight_event, eventweight_opp_event = geteventweights(template_name)
					if self.addTwice[0]==1 :
						eventweight_event*=0.5; eventweight_opp_event*=0.5
					#Add to all distributions
					for j in range(len(event_counts[ltype][ttype])) :
						eventweight, eventweight_opp = eventweight_event, eventweight_opp_event
						if tree_cuts[j] == 'sm_ps' and not (sm and ps) :
							continue
						elif tree_cuts[j] == 'lm1_ps' and not (lm1 and ps) :
							continue
						elif tree_cuts[j] == 'lm2_ps' and not (lm2 and ps) :
							continue
						elif tree_cuts[j] == 'hm_ps' and not (hm and ps) :
							continue
						elif tree_cuts[j] == 'lm1_fs' and not (lm1 and fs) :
							continue
						elif tree_cuts[j] == 'lm2_fs' and not (lm2 and fs) :
							continue
						elif tree_cuts[j] == 'hm_fs' and not (hm and fs) :
							continue
						funcweight, funcweight_opp = getfuncweights(template_name,NGGS[j],NG1S[j],NG2S[j],NG3S[j],NG4S[j],NQQS[j],NQ1S[j],NQ2S[j],NBCKS[j],NNTMJS[j])
						eventweight*=funcweight; eventweight_opp*=funcweight_opp
						event_counts[ltype][ttype][j]+=eventweight
						if self.addTwice[0]==1 :
							event_counts[ltype][ttype][j]+=eventweight_opp

	def apply_conversion_functions(self,distlist,conv_functions) :
		#Set Branches
		for branch in self.all_branches :
			self.antitagged_tree.SetBranchAddress(branch[1],branch[2])
		nEntries = self.antitagged_tree.GetEntries()
		#For all of the NTMJ templates
		for i in range(len(self.all_templates)) :
			template_name = self.all_templates[i].name
			print '			Doing template '+template_name
			for ltype in conv_functions :
				if self.name.split('__')[0].split('plus')[0] != ltype and self.name.split('__')[0].split('minus')[0] != ltype :
					continue
				for ttype in conv_functions[ltype] :
					if (len(self.name.split('__'))==2 and ttype=='nominal') or (len(self.name.split('__'))>2 and ttype==self.name.split('__')[2]+'__'+self.name.split('__')[3]) :
						#Get the event numbers from other distributions
						NGG,NG1,NG2,NG3,NG4,NQQ,NQ1,NQ2,NBCK,NNTMJ = self.__get_total_event_numbers__(distlist,'sm_fs',template_name)
						NNTMJ = 0. #because we're going to find what it is
						#Loop optimization
						getentry_here = self.antitagged_tree.GetEntry
						geteventweights_here = self.__get_event_weights__
						getfuncweights_here = self.__get_func_weights__
						for entry in range(nEntries) :
							getentry_here(entry)
							eventweight, eventweight_opp = geteventweights_here(template_name)
							funcweight, funcweight_opp = getfuncweights_here(template_name,NGG,NG1,NG2,NG3,NG4,NQQ,NQ1,NQ2,NBCK,NNTMJ)
							eventweight*=funcweight; eventweight_opp*=funcweight_opp
							#apply the conversion function
							conv_value = 1.0
							if template_name.find('__fit__up')!=-1 :
								conv_value = conv_functions[ltype][ttype][1].Eval(self.hadt_M[0])
							elif template_name.find('__fit__down')!=-1 :
								conv_value = conv_functions[ltype][ttype][2].Eval(self.hadt_M[0])
							elif template_name.find('__fit__')==-1 :
								conv_value = conv_functions[ltype][ttype][0].Eval(self.hadt_M[0])
							eventweight*=conv_value; eventweight_opp*=conv_value
							if self.addTwice[0]==1 :
								eventweight*=0.5; eventweight_opp*=0.5
							self.all_templates[i].Fill(self.cstar[0],abs(self.x_F[0]),self.M[0],eventweight)
							if self.addTwice[0]==1 :
								self.all_templates[i].Fill(-1.0*self.cstar[0],abs(self.x_F[0]),self.M[0],eventweight_opp)

	def __add_all_branches__(self) :
		print '				Adding observables branches'
		self.cstar 				  = array('d',[100.]);  self.all_branches.append(('cstar_scaled','cstar',self.cstar,'/D'))
		self.x_F 				  = array('d',[100.]);  self.all_branches.append(('x_F_scaled','x_F',self.x_F,'/D'))
		self.M 					  = array('d',[-1.0]);  self.all_branches.append(('M_scaled','M',self.M,'/D'))
		self.hadt_M 			  = array('d',[-1.0]);  self.all_branches.append(('scaled_hadt_M','hadt_M',self.hadt_M,'/D'))
		self.hadt_tau32 		  = array('d',[-1.0]);  self.all_branches.append(('hadt_tau32','hadt_tau32',self.hadt_tau32,'/D'))
		self.Q_l 				  = array('i',[0]);     self.all_branches.append(('Q_l','Q_l',self.Q_l,'/I'))
		self.addTwice 			  = array('I',[2]);     self.all_branches.append(('addTwice','addTwice',self.addTwice,'/i'))
		self.contrib_weight 	  = array('d',[1.0]);   self.all_branches.append((None,'contrib_weight',self.contrib_weight,'/D'))
		self.dist_reweight_names = []
		self.dist_reweight_arrays = [] 
		self.dist_reweight_opp_arrays = []
		self.const_reweights_arrays = []
		self.simple_systematic_reweights_arrays = []
		self.pdf_reweight_array = array('d',[1.0])
		self.fit_function_reweight_arrays = []
		self.fit_parameter_names = []
		#only MC distribtions need other branches 
		if self.name.find('DATA') == -1 :
			#constant reweights
			for i in range(len(const_reweights_trees)) :
				print '				Adding branch for constant reweight '+const_reweights_dists[i]
				self.const_reweights_arrays.append(array('d',[1.0])) 
				self.all_branches.append((const_reweights_trees[i],const_reweights_dists[i],self.const_reweights_arrays[i],'/D'))
			#simple systematics
			for i in range(len(simple_systematics_trees)) :
				print '				Adding branches for simple systematic reweight '+simple_systematics_dists[i]
				if simple_systematics_trees[i]=='luminosity' :
					self.simple_systematic_reweights_arrays.append(array('d',[1.0]))
					self.all_branches.append((None,simple_systematics_dists[i],self.simple_systematic_reweights_arrays[3*i],'/D'))
					self.simple_systematic_reweights_arrays.append(array('d',[1.0+LUMI_MIN_BIAS]))
					self.all_branches.append((None,simple_systematics_dists[i]+'_up',self.simple_systematic_reweights_arrays[3*i+1],'/D'))
					self.simple_systematic_reweights_arrays.append(array('d',[1.0-LUMI_MIN_BIAS]))
					self.all_branches.append((None,simple_systematics_dists[i]+'_down',self.simple_systematic_reweights_arrays[3*i+2],'/D'))
				else :
					self.simple_systematic_reweights_arrays.append(array('d',[1.0]))
					self.all_branches.append((simple_systematics_trees[i],simple_systematics_dists[i],self.simple_systematic_reweights_arrays[3*i],'/D'))
					self.simple_systematic_reweights_arrays.append(array('d',[1.0]))
					self.all_branches.append((simple_systematics_trees[i]+'_hi',simple_systematics_dists[i]+'_up',self.simple_systematic_reweights_arrays[3*i+1],'/D'))
					self.simple_systematic_reweights_arrays.append(array('d',[1.0]))
					self.all_branches.append((simple_systematics_trees[i]+'_low',simple_systematics_dists[i]+'_down',self.simple_systematic_reweights_arrays[3*i+2],'/D'))
			if self.include_PDF :
				#pdf systematics vectors
				print '				Adding branch for pdf systematic reweight vector '+pdf_reweight_vector_dists
				self.pdf_reweight_array = array('d',53*[1.0])
				self.all_branches.append((pdf_reweight_vector_trees,pdf_reweight_vector_dists,self.pdf_reweight_array,'[53]/D'))
			#fit parameters and distribution reweights
			#read the numerical values in from the parameter file
			pars = []
			parfile = open(self.parfile_name)
			for line in parfile :
				if line.startswith('#') :
					continue
				a = line.rstrip().split()
				if len(a) == 4 :
					[parname,initialvalue,downvalue,upvalue] = a
				elif len(a) == 5 :
					[parname,initialvalue,downvalue,upvalue,sigma] = a
				pars.append([parname,float(initialvalue),float(downvalue),float(upvalue)])
			#find which portions of the function string are fit parameters or distribution reweights
			args = []; argispar = []
			if self.function != None :
				args = self.function.split('#')
			for j in range(len(args)) : #replace the parameters within the arguments
				argispar.append(False)
				if args[j].startswith('wg') or args[j].startswith('wq') :
					self.dist_reweight_names.append(args[j])
					self.dist_reweight_arrays.append(array('d',[1.0])); self.dist_reweight_opp_arrays.append(array('d',[1.0]))
					self.all_branches.append((args[j],args[j],self.dist_reweight_arrays[len(self.dist_reweight_arrays)-1],'/D'))
					self.all_branches.append((args[j]+'_opp',args[j]+'_opp',self.dist_reweight_opp_arrays[len(self.dist_reweight_opp_arrays)-1],'/D'))
				else :
					for par in pars :
						if par[0] == args[j] :
							args[j] = par
							argispar[j] = True
							if par[0] not in self.fit_parameter_names :
								self.fit_parameter_names.append(par[0])
							break
			#nominal distribution
			factorstring = ''
			for j in range(len(args)) :
				if not argispar[j] :
					if args[j].startswith('N') or args[j].startswith('wg') or args[j].startswith('wq') :
						factorstring+='#'+args[j]+'#'
					else :
						factorstring+=args[j]
				else :
					factorstring+='('+str(args[j][1])+')'
			if factorstring != '' :
				print '				Adding nominal fit parameter function: '+factorstring
			self.fit_function_reweight_arrays.append(factorstring)
			if self.name.find('JES')==-1 and self.name.find('JER')==-1 :
				#all the up and down distributions
				for k in range(len(self.fit_parameter_names)) :
					factorstring_down = ''; factorstring_up = ''
					thisparname = ''
					for j in range(len(args)) :
						if not argispar[j] :
							if args[j].startswith('N') or args[j].startswith('wg') or args[j].startswith('wq') :
								factorstring_down+='#'+args[j]+'#'; factorstring_up+='#'+args[j]+'#'
							else :
								factorstring_down+=args[j]; factorstring_up+=args[j]
						elif argispar[j] and self.fit_parameter_names[k]!=args[j][0] :
							factorstring_down+='('+str(args[j][1])+')'; factorstring_up+='('+str(args[j][1])+')'
						elif argispar[j] and self.fit_parameter_names[k]==args[j][0] :
							thisparname = args[j][0]
							factorstring_down+='('+str(args[j][2])+')'; factorstring_up+='('+str(args[j][3])+')'
					print '				Adding '+self.fit_parameter_names[k]+' up function: '+factorstring_up
					self.fit_function_reweight_arrays.append(factorstring_up)
					print '				Adding '+self.fit_parameter_names[k]+' down function: '+factorstring_down
					self.fit_function_reweight_arrays.append(factorstring_down)
		#finally add all the new branches to the tree
		for branch in self.all_branches :
			self.signal_tree.Branch(branch[1],branch[2],branch[1]+branch[3])
			self.antitagged_tree.Branch(branch[1],branch[2],branch[1]+branch[3])
			self.sidebands_tree.Branch(branch[1],branch[2],branch[1]+branch[3])

	def __add_all_templates__(self) :
		#nominal
		if self.name.find('DATA')!=-1 or self.name.find('JES') != -1 or self.name.find('JER') != -1 :
			self.all_templates.append(template(self.name,self.formatted_name))
		else :
			self.all_templates.append(template(self.name,self.formatted_name+', nominal'))
		if self.name.find('JES') == -1 and self.name.find('JER') == -1 :
			if self.step == 'initial' :
				#simple systematics up/down
				for i in range(len(self.simple_systematic_reweights_arrays)/3) :
					sysname = simple_systematics_dists[i]
					self.all_templates.append(template(self.name+'__'+sysname+'__up',self.formatted_name+', '+sysname+' up'))
					self.all_templates.append(template(self.name+'__'+sysname+'__down',self.formatted_name+', '+sysname+' down'))
				#PDF systematics up/down
				for i in range(1,(len(self.pdf_reweight_array)-1)/2) :
					self.all_templates.append(template(self.name+'__pdf_lambda_'+str(i)+'__up',self.formatted_name+', PDF lambda'+str(i)+' up'))
					self.all_templates.append(template(self.name+'__pdf_lambda_'+str(i)+'__down',self.formatted_name+', PDF lambda'+str(i)+' down'))
			if self.step != 'final' :
				#fit parameters up/down
				for i in range(len(self.fit_parameter_names)) :
					self.all_templates.append(template(self.name+'__par_'+self.fit_parameter_names[i]+'__up',self.formatted_name+', '+self.fit_parameter_names[i]+' up'))
					self.all_templates.append(template(self.name+'__par_'+self.fit_parameter_names[i]+'__down',self.formatted_name+', '+self.fit_parameter_names[i]+' down'))
				#NTMJ fit parameters up/down
				if self.name.find('ntmj')!=-1 :
					self.all_templates.append(template(self.name+'__fit__up',self.formatted_name+', NTMJ fit error up'))
					self.all_templates.append(template(self.name+'__fit__down',self.formatted_name+', NTMJ fit error down'))

	def  __get_total_event_numbers__(self,ds,cuts,templatename) :
		ns = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
		names = ['fgg_0','fgg_1','fgg_2','fgg_3','fgg_4','fqq_0','fqq_1','fqq_2','fbck','fntmj']
		rws   = [None,   'wg1',  'wg2',  'wg3',  'wg4',  None,   'wqs1', 'wqs2', None,  None]
		#check the cut string
		f, s = cuts.split('_')[0], cuts.split('_')[1]
		if not (f == 'lm1' or f == 'lm2' or f =='sm' or f == 'hm') or not (s == 'ps' or s == 'fs') :
			print 'WARNING: CUTSTRING NOT RECOGNIZED: '+cuts
			return
		#find the right trees
		for i in range(len(ds)) :
			#check that the distribution is in the right channel
			thisdistsplit = self.name.split('__')
			if ds[i].name.split('__')[0]!=thisdistsplit[0] :
				continue
			if len(thisdistsplit)>2 and (len(ds[i].name.split('__'))<4 or ds[i].name.split('__')[2]!=thisdistsplit[2] or ds[i].name.split('__')[3]!=thisdistsplit[3]) :
				continue
			#get each type of number
			for j in range(len(names)) :
				if names[j]=='fntmj' and ds[i].name.split('__')[1].find(names[j])!=-1 :
					for k in range(len(ds[i].all_templates)) :
						thistempnamesplit = ds[i].all_templates[k].name.split('__')
						origtempnamesplit = templatename.split('__')
						if (len(thistempnamesplit)==4 and len(origtempnamesplit)==4 and thistempnamesplit[2]==origtempnamesplit[2] and thistempnamesplit[3]==origtempnamesplit[3]) :
							ns[j] = ds[i].all_templates[k].histo_3D.Integral()
							break
						elif len(thistempnamesplit)<4 : 
							ns[j] = ds[i].all_templates[k].histo_3D.Integral()
				elif ds[i].name.split('__')[1].find(names[j].split('_')[0])!=-1 :
					#Set branches
					tree = ds[i].sidebands_tree
					if f == 'sm' and s == 'ps' :
						tree = ds[i].signal_tree
					elif f == 'sm' and s == 'fs' :
						tree = ds[i].antitagged_tree
					for branch in ds[i].all_branches :
						tree.SetBranchAddress(branch[1],branch[2])
					dist_rw_branch = array('d',[1.0])
					dist_rw_opp_branch = array('d',[1.0])
					for k in range(len(ds[i].dist_reweight_names)) :
						if ds[i].dist_reweight_names[k] == rws[j] :
							dist_rw_branch = ds[i].dist_reweight_arrays[k]
							dist_rw_opp_branch = ds[i].dist_reweight_opp_arrays[k]
							break
					#loop over events
					nEntries = tree.GetEntries()
					#some loop optimization
					getentry_local = tree.GetEntry
					geteventweights_local = ds[i].__get_event_weights__
					for entry in range(nEntries) :
						getentry_local(entry)
						hadtM_local = ds[i].hadt_M[0]; hadttau32_local = ds[i].hadt_tau32[0]
						lm1_local = hadtM_local>LM1_LOW and hadtM_local<=LM1_HI
						lm2_local = hadtM_local>LM2_LOW and hadtM_local<=LM2_HI
						hm_local  = hadtM_local>HM_LOW and hadtM_local<=HM_HI
						sm_local  = hadtM_local>SM_LOW and hadtM_local<=SM_HI
						ps_local  = hadttau32_local<SUBS_CUT
						fs_local  = hadttau32_local>SUBS_CUT
						if f == 'sm' :
							if cuts == 'sm_ps' and not (sm_local and ps_local) :
								continue
							elif cuts == 'sm_fs' and not (sm_local and fs_local) :
								continue
						else :
							if cuts == 'lm1_ps' and not (lm1_local and ps_local) :
								continue
							elif cuts == 'lm2_ps' and not (lm2_local and ps_local) :
								continue
							elif cuts == 'hm_ps' and not (hm_local and ps_local) :
								continue
							elif cuts == 'lm1_fs' and not (lm1_local and fs_local) :
								continue
							elif cuts == 'lm2_fs' and not (lm2_local and fs_local) :
								continue
							elif cuts == 'hm_fs' and not (hm_local and fs_local) :
								continue
						#get event weights
						eventweight, eventweight_opp = geteventweights_local(templatename)
						eventweight*=dist_rw_branch[0]; eventweight_opp*=dist_rw_opp_branch[0]
						if ds[i].addTwice[0]==1 :
							eventweight*=0.5; eventweight_opp*=0.5
						ns[j]+=eventweight
						if ds[i].addTwice[0]==1 :
							ns[j]+=eventweight_opp
		return ns[0],ns[1],ns[2],ns[3],ns[4],ns[5],ns[6],ns[7],ns[8],ns[9]

	def __get_event_weights__(self,templatename) :
		#Lumi/cross section reweight
		eweight = self.contrib_weight[0]; eweight_opp = self.contrib_weight[0]
#		print '---------------------------------NEW EVENT--------------------------------' #DEBUG
#		print '	eweight = contrib weight = '+str(self.contrib_weight[0]) #DEBUG
		#data doesn't need anything else
		if templatename.find('DATA')!=-1 :
			return eweight, eweight_opp
		#Constant reweights
		for reweight in self.const_reweights_arrays :
			eweight*=reweight[0]; eweight_opp*=reweight[0]
#			print '	eweight *= constreweight ('+str(reweight[0])+') = '+str(eweight) #DEBUG
		#Simple systematics
		ssras = self.simple_systematic_reweights_arrays
		for j in range(len(self.simple_systematic_reweights_arrays)/3) :
			thissys = simple_systematics_dists[j]
			if templatename.find(thissys+'__up')!=-1 : 
				eweight*=ssras[3*j+1][0]
#				print '	eweight *= ssras['+str(3*j+1)+'][0] ('+str(ssras[3*j+1][0])+') = '+str(eweight) #DEBUG 
				eweight_opp*=ssras[3*j+1][0]
			elif templatename.find(thissys+'__down')!=-1 :
				eweight*=ssras[3*j+2][0]
#				print '	eweight *= ssras['+str(3*j+2)+'][0] ('+str(ssras[3*j+2][0])+') = '+str(eweight) #DEBUG 
				eweight_opp*=ssras[3*j+2][0]
			else :
				eweight*=ssras[3*j][0]
#				print '	eweight *= ssras['+str(3*j)+'][0] ('+str(ssras[3*j][0])+') = '+str(eweight) #DEBUG 
				eweight_opp*=ssras[3*j][0]
		#PDF systematics
		pra = self.pdf_reweight_array
		tnamesplit = templatename.split('pdf_lambda_')
		if len(tnamesplit) > 1 :
			ipdf = eval(tnamesplit[len(tnamesplit)-1].split('__')[0])
			upordown = tnamesplit[len(tnamesplit)-1].split('__')[1]
			if upordown=='up' :
				eweight*=pra[2*(ipdf-1)+1]
#				print '	eweight *= pra['+str(2*(ipdf-1)+1)+'] ('+str(pra[2*(ipdf-1)+1])+') = '+str(eweight) #DEBUG 
				eweight_opp*=pra[2*(ipdf-1)+1]
			elif upordown=='down' :
				eweight*=pra[2*(ipdf-1)+2]
#				print '	eweight *= pra['+str(2*(ipdf-1)+2)+'] ('+str(pra[2*(ipdf-1)+2])+') = '+str(eweight) #DEBUG 
				eweight_opp*=pra[2*(ipdf-1)+2]
		else :
			eweight*=pra[0]
#			print '	eweight *= pra[0] ('+str(pra[0])+') = '+str(eweight) #DEBUG
			eweight_opp*=pra[0]
		return eweight,eweight_opp

	def __get_func_weights__(self,templatename,ngg,ng1,ng2,ng3,ng4,nqq,nq1,nq2,nbck,nntmj) :
		eweight = 1.0; eweight_opp = 1.0
		#fit parameter function weights
		#replace the total event numbers in the function strings
		real_rw_arrays = []; real_rw_opp_arrays = []
		ffras = self.fit_function_reweight_arrays
		for i in range(len(ffras)) :
			newfunction = ''; newfunction_opp = ''
			arraysplit = ffras[i].split('#')
			for j in range(len(arraysplit)) :
				check = arraysplit[j]
				if check == 'NTOT' :
					newfunction+='('+str(ngg+nqq+nbck+nntmj)+')'; newfunction_opp+='('+str(ngg+nqq+nbck+nntmj)+')'
				elif check == 'NBCK' :
					newfunction+='('+str(nbck)+')'; newfunction_opp+='('+str(nbck)+')'
				elif check == 'NNTMJ' :
					newfunction+='('+str(nntmj)+')'; newfunction_opp+='('+str(nntmj)+')'
				elif check == 'NTTBAR' :
					newfunction+='('+str(ngg+nqq)+')'; newfunction_opp+='('+str(ngg+nqq)+')'
				elif check == 'NQQBAR' :
					newfunction+='('+str(nqq)+')'; newfunction_opp+='('+str(nqq)+')'
				elif check == 'NG1' :
					newfunction+='('+str(ng1)+')'; newfunction_opp+='('+str(ng1)+')'
				elif check == 'NG2' :
					newfunction+='('+str(ng2)+')'; newfunction_opp+='('+str(ng2)+')'
				elif check == 'NG3' :
					newfunction+='('+str(ng3)+')'; newfunction_opp+='('+str(ng3)+')'
				elif check == 'NG4' :
					newfunction+='('+str(ng4)+')'; newfunction_opp+='('+str(ng4)+')'
				elif check == 'NQ1' :
					newfunction+='('+str(nq1)+')'; newfunction_opp+='('+str(nq1)+')'
				elif check == 'NQ2' :
					newfunction+='('+str(nq2)+')'; newfunction_opp+='('+str(nq2)+')'
				elif check.startswith('wg') or check.startswith('wq') :
					for k in range(len(self.dist_reweight_names)) :
						if self.dist_reweight_names[k]==arraysplit[j] :
							newfunction+='('+str(self.dist_reweight_arrays[k][0])+')'; newfunction_opp+='('+str(self.dist_reweight_opp_arrays[k][0])+')'
				else :
					newfunction+=arraysplit[j]; newfunction_opp+=arraysplit[j]
			real_rw_arrays.append(eval(newfunction)); real_rw_opp_arrays.append(eval(newfunction_opp))
		for j in range(len(self.fit_parameter_names)) :
			if len(templatename.split('__'))>2 and templatename.split('__')[2].find(self.fit_parameter_names[j])!=-1 :
				if templatename.split('__')[3].find('up')!=-1 :
					eweight*=real_rw_arrays[2*j+1]
#					print '	eweight *= real_rw_arrays['+str(2*j+1)+'] ('+str(real_rw_arrays[2*j+1])+') = '+str(eweight) #DEBUG 
					eweight_opp*=real_rw_opp_arrays[2*j+1]
					break
				elif templatename.split('__')[3].find('down')!=-1 :
					eweight*=real_rw_arrays[2*j+2]
#					print '	eweight *= real_rw_arrays['+str(2*j+2)+'] ('+str(real_rw_arrays[2*j+2])+') = '+str(eweight) #DEBUG 
					eweight_opp*=real_rw_opp_arrays[2*j+2]
					break
			elif j == len(self.fit_parameter_names)-1 :
				eweight*=real_rw_arrays[0]
#				print '	eweight *= real_rw_arrays[0] ('+str(real_rw_arrays[0])+') = '+str(eweight) #DEBUG 
				eweight_opp*=real_rw_opp_arrays[0]
		return eweight,eweight_opp
