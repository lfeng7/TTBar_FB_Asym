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
		#pretagged tree and branches
		self.tree = TTree(self.name+'_tree',self.name+'_tree')
		self.tree.SetDirectory(0)
		self.all_branches = []
		print '			Adding branches'
		self.__add_all_branches__()
		#automatically generate the list of templates that will be made from this distribution
		self.all_templates = []
		print '			Adding templates'
		self.__add_all_templates__()
		print '		Done'

	#build_templates builds the 3D templates for a distribution, and returns a list of 1D templates
	def build_templates(self,distlist) :
		#set the branches to read
		for branch in self.all_branches :
			self.tree.SetBranchAddress(branch[1],branch[2])
		nEntries = self.tree.GetEntries()
		#read and add events from the tree to all of the derived templates
		for i in range(len(self.all_templates)) :
			template_name = self.all_templates[i].name
			print '		Doing template '+template_name
			#Get the total numbers of events for 
			NGG,NG1,NG2,NG3,NG4,NQQ,NQ1,NQ2,NBCK,NNTMJ = self.__get_total_event_numbers__(distlist,'hadt_M>140. && hadt_M<=250. && hadt_tau32<0.55',template_name)
#			print '		NGG, NG1, NG2, NG3, NG4, NQQ, NQ1, NQ2, NBCK, NNTMJ = %f, %f, %f, %f, %f, %f, %f, %f, %f, %f'%(NGG,NG1,NG2,NG3,NG4,NQQ,NQ1,NQ2,NBCK,NNTMJ) #DEBUG
			for entry in range(nEntries) :
				self.tree.GetEntry(entry)
				if self.hadt_M[0]<140. or self.hadt_M[0]>250. or self.hadt_tau32[0]>0.55 :
					continue
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

	def fix_NTMJ_templates(self,f_aux,distlist) :
		#start by splitting the tree into several
		all_tree_cuts = []
		all_tree_cuts.append('hadt_M>100. && hadt_M<=120. && hadt_tau32<0.55')
		all_tree_cuts.append('hadt_M>120. && hadt_M<=140. && hadt_tau32<0.55')
		all_tree_cuts.append('hadt_M>250. && hadt_M<=500. && hadt_tau32<0.55')
		all_tree_cuts.append('hadt_M>100. && hadt_M<=120. && hadt_tau32>0.55')
		all_tree_cuts.append('hadt_M>120. && hadt_M<=140. && hadt_tau32>0.55')
		all_tree_cuts.append('hadt_M>250. && hadt_M<=500. && hadt_tau32>0.55')
		all_tree_cuts.append('hadt_M>140. && hadt_M<=250. && hadt_tau32>0.55')
		all_event_counts = []
		for i in range(6) :
			all_event_counts.append([])
		#Set Branches
		for branch in self.all_branches :
			self.tree.SetBranchAddress(branch[1],branch[2])
		nEntries = self.tree.GetEntries()
		#For all of the NTMJ templates
		for i in range(len(self.all_templates)) :
			template_name = self.all_templates[i].name
			print '		Doing template '+template_name
			#Get the total number of events in each of the sideband trees
			for j in range(len(all_event_counts)) :
				#Get the event numbers from other distributions
				NGG,NG1,NG2,NG3,NG4,NQQ,NQ1,NQ2,NBCK,NNTMJ = self.__get_total_event_numbers__(distlist,all_tree_cuts[j],template_name)
				NNTMJ = 0. #because we're going to find what it is
				all_event_counts[j].append(0.)
				for entry in range(nEntries) :
					self.tree.GetEntry(entry)
					if all_tree_cuts[j] == 'hadt_M>140. && hadt_M<=250. && hadt_tau32<0.55' :
						if not (self.hadt_M[0] > 140. and self.hadt_M[0] <= 250. and self.hadt_tau32[0] < 0.55) :
							continue
					elif all_tree_cuts[j] == 'hadt_M>100. && hadt_M<=120. && hadt_tau32<0.55' :
						if not (self.hadt_M[0]>100. and self.hadt_M[0]<=120. and self.hadt_tau32[0] < 0.55) :
							continue
					elif all_tree_cuts[j] == 'hadt_M>120. && hadt_M<=140. && hadt_tau32<0.55' :
						if not (self.hadt_M[0]>120. and self.hadt_M[0]<=140. and self.hadt_tau32[0] < 0.55) :
							continue
					elif all_tree_cuts[j] == 'hadt_M>250. && hadt_M<=500. && hadt_tau32<0.55' :
						if not (self.hadt_M[0]>250. and self.hadt_M[0]<=500. and self.hadt_tau32[0] < 0.55) :
							continue
					elif all_tree_cuts[j] == 'hadt_M>100. && hadt_M<=120. && hadt_tau32>0.55' :
						if not (self.hadt_M[0]>100. and self.hadt_M[0]<=120. and self.hadt_tau32[0] > 0.55) :
							continue
					elif all_tree_cuts[j] == 'hadt_M>120. && hadt_M<=140. && hadt_tau32>0.55' :
						if not (self.hadt_M[0]>120. and self.hadt_M[0]<=140. and self.hadt_tau32[0] > 0.55) :
							continue
					elif all_tree_cuts[j] == 'hadt_M>250. && hadt_M<=500. && hadt_tau32>0.55' :
						if not (self.hadt_M[0]>250. and self.hadt_M[0]<=500. and self.hadt_tau32[0] > 0.55) :
							continue
					else :
						print 'CUTSTRING NOT RECOGNIZED ('+all_tree_cuts[j]+')'
						continue
					eventweight, eventweight_opp = self.__get_event_weights__(template_name)
					funcweight, funcweight_opp = self.__get_func_weights__(template_name,NGG,NG1,NG2,NG3,NG4,NQQ,NQ1,NQ2,NBCK,NNTMJ)
					eventweight*=funcweight; eventweight_opp*=funcweight_opp
					if self.addTwice[0]==1 :
						eventweight*=0.5; eventweight_opp*=0.5
					all_event_counts[j][i]+=eventweight
					if self.addTwice[0]==1 :
						all_event_counts[j][i]+=eventweight_opp
				if all_event_counts[j][i] <= 0. :
					all_event_counts[j][i] = 1.0
			print 'NUMBERS OF EVENTS: ' #DEBUG
			for j in range(len(all_event_counts)) : #DEBUG
				print '		'+str(all_event_counts[j][i]) #DEBUG
			#Find the y values and errors for the points on the graph
			l1y = all_event_counts[0][i]/all_event_counts[3][i]
			l2y = all_event_counts[1][i]/all_event_counts[4][i]
			hy  = all_event_counts[2][i]/all_event_counts[5][i]
			l1ye = l1y*sqrt((1./all_event_counts[3][i])+(1./all_event_counts[0][i]))
			l2ye = l2y*sqrt((1./all_event_counts[4][i])+(1./all_event_counts[1][i]))
			hye  = hy*sqrt((1./all_event_counts[5][i])+(1./all_event_counts[2][i]))
			print '		l1y = '+str(l1y) #DEBUG
			print '		l2y = '+str(l2y) #DEBUG
			print '		hy = '+str(hy) #DEBUG
			print '		l1ye = '+str(l1ye) #DEBUG
			print '		l2ye = '+str(l2ye) #DEBUG
			print '		hye = '+str(hye) #DEBUG
			#Build the TGraph to fit with the conversion function
			n=3
			xs  = array('d',[110.,130.,375.])
			xes = array('d',[10.,10.,125.])
			ys  = array('d',[l1y,l2y,hy])
			yes = array('d',[l1ye,l2ye,hye])
			gr = TGraphErrors(n,xs,ys,xes,yes)
			#Define the conversion function
			conv_func = TF1('conv_func','[0]*x+[1]',100.,500.)
			#Fit the graph
			gr.Fit('conv_func')
			#get the fitter
			fitter = TVirtualFitter.GetFitter()
			#Build the nominal, slope up/down, and intercept up/down functions
			nom_func = TF1('nom_func','[0]*x+[1]',100.,500.) 
			nom_func.SetParameter(0,conv_func.GetParameter(0)) 
			nom_func.SetParameter(1,conv_func.GetParameter(1))
			fit_up_func = TF1('fit_up_func','[0]*x+[1]+sqrt([2]*[2]+2*x*[3]*[3]+x*x*[4]*[4])',100.,500.) 
			fit_up_func.SetParameter(0,conv_func.GetParameter(0)) 
			fit_up_func.SetParameter(1,conv_func.GetParameter(1))
			fit_up_func.SetParameter(2,conv_func.GetParError(1))
			fit_up_func.SetParameter(3,fitter.GetCovarianceMatrixElement(0,1))
			fit_up_func.SetParameter(4,conv_func.GetParError(0))
			fit_down_func = TF1('fit_down_func','[0]*x+[1]-sqrt([2]*[2]+2*x*[3]*[3]+x*x*[4]*[4])',100.,500.) 
			fit_down_func.SetParameter(0,conv_func.GetParameter(0)) 
			fit_down_func.SetParameter(1,conv_func.GetParameter(1))
			fit_down_func.SetParameter(2,conv_func.GetParError(1))
			fit_down_func.SetParameter(3,fitter.GetCovarianceMatrixElement(0,1))
			fit_down_func.SetParameter(4,conv_func.GetParError(0))
			#Build the template from the converted antitagged tree
			plot_func = conv_func
			if template_name.find('__fit_')==-1 :
				plot_func = nom_func
			elif template_name.find('__fit__up')!=-1 :
				plot_func = fit_up_func
			elif template_name.find('__fit__down')!=-1 :
				plot_func = fit_down_func
			NGG,NG1,NG2,NG3,NG4,NQQ,NQ1,NQ2,NBCK,NNTMJ = self.__get_total_event_numbers__(distlist,all_tree_cuts[len(all_tree_cuts)-1],template_name)
			NNTMJ = 0. #because we're going to find what it is
#			print '		NGG, NG1, NG2, NG3, NG4, NQQ, NQ1, NQ2, NBCK, NNTMJ = %f, %f, %f, %f, %f, %f, %f, %f, %f, %f'%(NGG,NG1,NG2,NG3,NG4,NQQ,NQ1,NQ2,NBCK,NNTMJ) #DEBUG
			for entry in range(nEntries) :
				self.tree.GetEntry(entry)
				if self.hadt_M[0]<140. or self.hadt_M[0]>250. or self.hadt_tau32[0]<0.55 :
					continue
				eventweight, eventweight_opp = self.__get_event_weights__(template_name)
				funcweight, funcweight_opp = self.__get_func_weights__(template_name,NGG,NG1,NG2,NG3,NG4,NQQ,NQ1,NQ2,NBCK,NNTMJ)
				eventweight*=funcweight; eventweight_opp*=funcweight_opp
				#apply the conversion function
				conv_value = 1.0
				if template_name.find('__fit__up')!=-1 :
					conv_value = fit_up_func.Eval(self.hadt_M[0])
				elif template_name.find('__fit__down')!=-1 :
					conv_value = fit_down_func.Eval(self.hadt_M[0])
				elif template_name.find('__fit_')==-1 :
					conv_value = nom_func.Eval(self.hadt_M[0])
				eventweight*=conv_value; eventweight_opp*=conv_value
				if self.addTwice[0]==1 :
					eventweight*=0.5; eventweight_opp*=0.5
				self.all_templates[i].Fill(self.cstar[0],abs(self.x_F[0]),self.M[0],eventweight)
				if self.addTwice[0]==1 :
					self.all_templates[i].Fill(-1.0*self.cstar[0],abs(self.x_F[0]),self.M[0],eventweight_opp)
			#Make the plot of the fit
#			canv = TCanvas(template_name+'_conv_func_canv',template_name+' conversion function canvas',900,900)
#			gr.SetTitle(template_name+' conversion function fit')
#			gr.GetXaxis().SetTitle('hadronic top candidate mass (GeV)')
#			gr.GetYaxis().SetTitle('conversion rate (N_{passed}/N_{failed})')
#			gr.GetXaxis().SetRangeUser(100.,500.)
#			gr.SetMarkerStyle(21)
#			gr.Draw('AP')
#			plot_func.SetLineWidth(3)
#			plot_func.SetLineColor(kRed)
#			plot_func.Draw('L SAME')
#			leg = TLegend(0.62,0.67,0.9,0.9)
#			leg.AddEntry(gr,'measured rates','PE')
#			leg.AddEntry(plot_func,'linear fit','L')
#			leg.Draw()
#			f_aux.cd()
#			canv.Write()


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
				for par in pars :
					if par[0] == args[j] :
						args[j] = par
						argispar[j] = True
						if par[0] not in self.fit_parameter_names :
							self.fit_parameter_names.append(par[0])
						break
				if args[j].startswith('wg') or args[j].startswith('wq') :
					self.dist_reweight_names.append(args[j])
					self.dist_reweight_arrays.append(array('d',[1.0])); self.dist_reweight_opp_arrays.append(array('d',[1.0]))
					self.all_branches.append((args[j],args[j],self.dist_reweight_arrays[len(dist_reweight_arrays)-1],'/D'))
					self.all_branches.append((args[j]+'_opp',args[j]+'_opp',self.dist_reweight_opp_arrays[len(dist_reweight_opp_arrays)-1],'/D'))
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
					print '				Adding '+self.fit_parameter_names[k]+' up function: '+factorstring_up+' = '+factorstring_up
					self.fit_function_reweight_arrays.append(factorstring_up)
					print '				Adding '+self.fit_parameter_names[k]+' down function: '+factorstring_down+' = '+factorstring_down
					self.fit_function_reweight_arrays.append(factorstring_down)
		#finally add all the new branches to the tree
		for branch in self.all_branches :
			self.tree.Branch(branch[1],branch[2],branch[1]+branch[3])

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
					for branch in ds[i].all_branches :
						ds[i].tree.SetBranchAddress(branch[1],branch[2])
					dist_rw_branch = array('d',[1.0])
					dist_rw_opp_branch = array('d',[1.0])
					for k in range(len(ds[i].dist_reweight_names)) :
						if ds[i].dist_reweight_names[k] == rws[j] :
							dist_rw_branch = ds[i].dist_reweight_arrays[k]
							dist_rw_opp_branch = ds[i].dist_reweight_opp_arrays[k]
							break
					#loop over events
					nEntries = ds[i].tree.GetEntries()
					for entry in range(nEntries) :
						ds[i].tree.GetEntry(entry)
						if cuts == 'hadt_M>140. && hadt_M<=250. && hadt_tau32<0.55' :
							if not (ds[i].hadt_M[0] > 140. and ds[i].hadt_M[0] <= 250. and ds[i].hadt_tau32[0] < 0.55) :
								continue
						elif cuts == 'hadt_M>100. && hadt_M<=120. && hadt_tau32<0.55' :
							if not (ds[i].hadt_M[0]>100. and ds[i].hadt_M[0]<=120. and ds[i].hadt_tau32[0] < 0.55) :
								continue
						elif cuts == 'hadt_M>120. && hadt_M<=140. && hadt_tau32<0.55' :
							if not (ds[i].hadt_M[0]>120. and ds[i].hadt_M[0]<=140. and ds[i].hadt_tau32[0] < 0.55) :
								continue
						elif cuts == 'hadt_M>250. && hadt_M<=500. && hadt_tau32<0.55' :
							if not (ds[i].hadt_M[0]>250. and ds[i].hadt_M[0]<=500. and ds[i].hadt_tau32[0] < 0.55) :
								continue
						elif cuts == 'hadt_M>100. && hadt_M<=120. && hadt_tau32>0.55' :
							if not (ds[i].hadt_M[0]>100. and ds[i].hadt_M[0]<=120. and ds[i].hadt_tau32[0] > 0.55) :
								continue
						elif cuts == 'hadt_M>120. && hadt_M<=140. && hadt_tau32>0.55' :
							if not (ds[i].hadt_M[0]>120. and ds[i].hadt_M[0]<=140. and ds[i].hadt_tau32[0] > 0.55) :
								continue
						elif cuts == 'hadt_M>250. && hadt_M<=500. && hadt_tau32>0.55' :
							if not (ds[i].hadt_M[0]>250. and ds[i].hadt_M[0]<=500. and ds[i].hadt_tau32[0] > 0.55) :
								continue
						elif cuts == 'hadt_M>140. && hadt_M<=250. && hadt_tau32>0.55' :
							if not (ds[i].hadt_M[0]>140. and ds[i].hadt_M[0]<=250. and ds[i].hadt_tau32[0] > 0.55) :
								continue
						else :
							print 'CUTSTRING NOT RECOGNIZED ('+cuts+')'
						#get event weights
						eventweight, eventweight_opp = ds[i].__get_event_weights__(templatename)
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
		for j in range(len(self.simple_systematic_reweights_arrays)/3) :
			if templatename.find(simple_systematics_dists[j]+'__up')!=-1 : 
				eweight*=self.simple_systematic_reweights_arrays[3*j+1][0]
#				print '	eweight *= self.simple_systematic_reweights_arrays['+str(3*j+1)+'][0] ('+str(self.simple_systematic_reweights_arrays[3*j+1][0])+') = '+str(eweight) #DEBUG 
				eweight_opp*=self.simple_systematic_reweights_arrays[3*j+1][0]
			elif templatename.find(simple_systematics_dists[j]+'__down')!=-1 :
				eweight*=self.simple_systematic_reweights_arrays[3*j+2][0]
#				print '	eweight *= self.simple_systematic_reweights_arrays['+str(3*j+2)+'][0] ('+str(self.simple_systematic_reweights_arrays[3*j+2][0])+') = '+str(eweight) #DEBUG 
				eweight_opp*=self.simple_systematic_reweights_arrays[3*j+2][0]
			else :
				eweight*=self.simple_systematic_reweights_arrays[3*j][0]
#				print '	eweight *= self.simple_systematic_reweights_arrays['+str(3*j)+'][0] ('+str(self.simple_systematic_reweights_arrays[3*j][0])+') = '+str(eweight) #DEBUG 
				eweight_opp*=self.simple_systematic_reweights_arrays[3*j][0]
		#PDF systematics
		tnamesplit = templatename.split('pdf_lambda_')
		if len(tnamesplit) > 1 :
			ipdf = eval(tnamesplit[len(tnamesplit)-1].split('__')[0])
			upordown = tnamesplit[len(tnamesplit)-1].split('__')[1]
			if upordown=='up' :
				eweight*=self.pdf_reweight_array[2*(ipdf-1)+1]
#				print '	eweight *= self.pdf_reweight_array['+str(2*(ipdf-1)+1)+'] ('+str(self.pdf_reweight_array[2*(ipdf-1)+1])+') = '+str(eweight) #DEBUG 
				eweight_opp*=self.pdf_reweight_array[2*(ipdf-1)+1]
			elif upordown=='down' :
				eweight*=self.pdf_reweight_array[2*(ipdf-1)+2]
#				print '	eweight *= self.pdf_reweight_array['+str(2*(ipdf-1)+2)+'] ('+str(self.pdf_reweight_array[2*(ipdf-1)+2])+') = '+str(eweight) #DEBUG 
				eweight_opp*=self.pdf_reweight_array[2*(ipdf-1)+2]
		else :
			eweight*=self.pdf_reweight_array[0]
#			print '	eweight *= self.pdf_reweight_array[0] ('+str(self.pdf_reweight_array[0])+') = '+str(eweight) #DEBUG
			eweight_opp*=self.pdf_reweight_array[0]
		return eweight,eweight_opp

	def __get_func_weights__(self,templatename,ngg,ng1,ng2,ng3,ng4,nqq,nq1,nq2,nbck,nntmj) :
		eweight = 1.0; eweight_opp = 1.0
		#fit parameter function weights
		#replace the total event numbers in the function strings
		real_rw_arrays = []; real_rw_opp_arrays = []
		for i in range(len(self.fit_function_reweight_arrays)) :
			newfunction = ''; newfunction_opp = ''
			arraysplit = self.fit_function_reweight_arrays[i].split('#')
			for j in range(len(arraysplit)) :
				if arraysplit[j] == 'NTOT' :
					newfunction+='('+str(ngg+nqq+nbck+nntmj)+')'; newfunction_opp+='('+str(ngg+nqq+nbck+nntmj)+')'
				elif arraysplit[j] == 'NBCK' :
					newfunction+='('+str(nbck)+')'; newfunction_opp+='('+str(nbck)+')'
				elif arraysplit[j] == 'NNTMJ' :
					newfunction+='('+str(nntmj)+')'; newfunction_opp+='('+str(nntmj)+')'
				elif arraysplit[j] == 'NTTBAR' :
					newfunction+='('+str(ngg+nqq)+')'; newfunction_opp+='('+str(ngg+nqq)+')'
				elif arraysplit[j] == 'NQQBAR' :
					newfunction+='('+str(nqq)+')'; newfunction_opp+='('+str(nqq)+')'
				elif arraysplit[j] == 'NG1' :
					newfunction+='('+str(ng1)+')'; newfunction_opp+='('+str(ng1)+')'
				elif arraysplit[j] == 'NG2' :
					newfunction+='('+str(ng2)+')'; newfunction_opp+='('+str(ng2)+')'
				elif arraysplit[j] == 'NG3' :
					newfunction+='('+str(ng3)+')'; newfunction_opp+='('+str(ng3)+')'
				elif arraysplit[j] == 'NG4' :
					newfunction+='('+str(ng4)+')'; newfunction_opp+='('+str(ng4)+')'
				elif arraysplit[j] == 'NQ1' :
					newfunction+='('+str(nq1)+')'; newfunction_opp+='('+str(nq1)+')'
				elif arraysplit[j] == 'NQ2' :
					newfunction+='('+str(nq2)+')'; newfunction_opp+='('+str(nq2)+')'
				elif arraysplit[j].startswith('wg') or arraysplit[j].startswith('wq') :
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


#convertTo1D takes a 3D distribution and makes it 1D for use with theta
def convertTo1D(original) :
	nBins = original.GetNbinsX()*original.GetNbinsY()*original.GetNbinsZ()
	newHisto = TH1F(original.GetName(),original.GetTitle(),nBins,0.,nBins-1.)
	newHisto.SetDirectory(0)
	for k in range(nBins) :
		newHisto.SetBinContent(k,original.GetBinContent(k))
	return newHisto