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

#TDR Style
gROOT.Macro('rootlogon.C')

##############################		   Template File Class  		##############################

class template_file :
	#docstring
	"""template_file class"""
	
	#__init__function
	def __init__(self,filename,parfilename,sumCharges) :
		#set the input variables
		self.f = TFile(filename,'recreate')
		self.f_aux = TFile(filename.rstrip('.root')+'_aux.root','recreate')
		self.parfile = parfilename
		self.sum_charge = sumCharges == 'yes'
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
		muon1_pt 	  = array('d',[-1.0]);  chain.SetBranchAddress('muon1_pt',muon1_pt)
		ele1_pt 	  = array('d',[-1.0]);  chain.SetBranchAddress('ele1_pt',ele1_pt)
		lepW_pt 	  = array('d',[-1.0]);  chain.SetBranchAddress('scaled_lepW_pt',lepW_pt)
		hadt_pt 	  = array('d',[-1.0]);  chain.SetBranchAddress('scaled_hadt_pt',hadt_pt)
		muon1_eta 	  = array('d',[100.0]); chain.SetBranchAddress('muon1_eta',muon1_eta)
		muon1_isLoose = array('I',[2]);  chain.SetBranchAddress('muon1_isLoose',muon1_isLoose)
		muon1_relPt   = array('d',[-1.0]);  chain.SetBranchAddress('muon1_relPt',muon1_relPt)
		muon1_dR 	  = array('d',[-1.0]);  chain.SetBranchAddress('muon1_dR',muon1_dR)
		ele1_eta 	  = array('d',[100.0]); chain.SetBranchAddress('ele1_eta',ele1_eta)
		ele1_isLoose  = array('I',[2]);  chain.SetBranchAddress('ele1_isLoose',ele1_isLoose)
		ele1_relPt 	  = array('d',[-1.0]);  chain.SetBranchAddress('ele1_relPt',ele1_relPt)
		ele1_dR 	  = array('d',[-1.0]);  chain.SetBranchAddress('ele1_dR',ele1_dR)
		lept_M 		  = array('d',[-1.0]);  chain.SetBranchAddress('scaled_lept_M',lept_M)
		hadt_tau21 	  = array('d',[-1.0]);  chain.SetBranchAddress('hadt_tau21',hadt_tau21)
		#find the distributions this tree will contribute to
		for i in range(len(self.dists)) :
			distname = self.dists[i].name
			#if events in this chain should be added to this distribution
			contribution = self.__contributesToDist__(file_ifd,i)
			if contribution != 0. :
				#set all of the branch addresses
				for branch in self.dists[i].all_branches :
					chain.SetBranchAddress(branch[0],branch[2])
				#copy over the relevant parts of the tree
				nEntries = chain.GetEntries()
				print 'Adding trees from '+ttree_dir_path+' to distribution '+distname+' ('+str(nEntries)+' entries)'
				added = 0
				for entry in range(nEntries) :
					chain.GetEntry(entry)
					keep = True
					muon_preselection = muon1_pt[0]>ele1_pt[0] and lepW_pt[0]>50. and hadt_pt[0]>300. and self.dists[i].hadt_M[0]>100.
					muon_kinematics   = muon1_pt[0]>40. and abs(muon1_eta[0])<2.4
					muon_ID = muon1_isLoose[0]==1
					muon_2D = muon1_relPt[0]>25. or muon1_dR[0]>0.5
					ele_preselection = ele1_pt[0]>muon1_pt[0] and lepW_pt[0]>50. and hadt_pt[0]>300. and self.dists[i].hadt_M[0]>100.
					ele_kinematics = ele1_pt[0]>40. and abs(ele1_eta[0])<2.4
					ele_ID = ele1_isLoose[0]==1
					ele_2D = ele1_relPt[0]>25. or ele1_dR[0]>0.5
					lep_top_mass = lept_M[0]>140. and lept_M[0]<250.
					muon_full_leptonic = muon_preselection and muon_kinematics and muon_ID and muon_2D and lep_top_mass
					muon_hadronic_pretag = muon_full_leptonic and hadt_tau21[0]>0.1
					ele_full_leptonic = ele_preselection and ele_kinematics and ele_ID and ele_2D and lep_top_mass
					ele_hadronic_pretag = ele_full_leptonic and hadt_tau21[0]>0.1
					signal_mass = self.dists[i].hadt_M[0]>140. and self.dists[i].hadt_M[0]<250.
					signal_tau32 = self.dists[i].hadt_tau32[0]<0.55 
					#make appropriate lepton cuts
					if distname.startswith('mu') :
						keep = muon_hadronic_pretag
					elif distname.startswith('el') :
						keep = ele_hadronic_pretag
					elif distname.startswith('allchannels') :
						keep = muon_hadronic_pretag or ele_hadronic_pretag
					if distname.find('plus')!=-1 :
						keep = keep and (self.dists[i].Q_l[0]>0 or self.dists[i].addTwice[0]==1)
					if distname.find('minus')!=-1 :
						keep = keep and (self.dists[i].Q_l[0]<0 or self.dists[i].addTwice[0]==1)
					if not keep :
						continue
					added+=1
					#change the event weight to add or subtract as appropriate
					self.dists[i].contrib_weight[0]=contribution
					self.dists[i].tree.Fill()
				print '	Added %d events'%(added)

	#build_templates builds the template histograms from the distributions
	def build_templates(self) :
		for dist in self.dists :
			if dist.name.find('fntmj')==-1 :
				print 'Building templates for distribution '+dist.name+'. . .'
				self.histos+=dist.build_templates(self.f_aux,self.sum_charge,self.parfile)
				print 'Done'

	#build_NTMJ_templates automatically generates templates for all of the data-driven NTMJ background
	#automatically accounting for systematics
	def build_NTMJ_templates(self) :
		ntmjdists = []
		for dist in self.dists :
			if dist.name.find('ntmj')!=-1 :
				ntmjdists.append(dist)
		for ntmjdist in ntmjdists :
			#first make room for all of the templates we'll eventually have
			#Nominal and due to fitting function parameters
			ntmjdist.__organizeFunctionTemplates__(self.parfile)
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
								if template_name.find(ntmjdist.systematics[k]+'__down')!=-1 :
									eventweight*=ntmjdist.systematics_arrays_down[k][0]; eventweight_opp*=ntmjdist.systematics_arrays_down[k][0]
								elif template_name.find(ntmjdist.systematics[k]+'__up')!=-1 : 
									eventweight*=ntmjdist.systematics_arrays_up[k][0]; eventweight_opp*=ntmjdist.systematics_arrays_up[k][0]
								else :
									eventweight*=ntmjdist.systematics_arrays[k][0]; eventweight_opp*=ntmjdist.systematics_arrays[k][0]
						#apply the function of the fitting parameters
						eventweight*=ntmjdist.factors[j]; eventweight_opp*=ntmjdist.factors[j]
						#if we want to add the event twice half the weight
						if ntmjdist.addTwice[0]==1 :
							eventweight*=0.5; eventweight_opp*=0.5
						#add to the running total
						if self.sum_charge or ntmjdist.name.startswith('allchannels') or (ntmjdist.Q_l[0]>0 and template_name.find('plus')!=-1) or (ntmjdist.Q_l[0]<0 and template_name.find('minus')!=-1) :
							sideband_counts[i][j]+=eventweight; sideband_xs[i][j]+=eventweight*ntmjdist.hadt_M[0]
						if ntmjdist.addTwice[0]==1 and (self.sum_charge or ntmjdist.name.startswith('allchannels') or (ntmjdist.Q_l[0]<0 and template_name.find('plus')!=-1) or (ntmjdist.Q_l[0]>0 and template_name.find('minus')!=-1)) :
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
							if template_name.find(ntmjdist.systematics[k]+'__down')!=-1 :
								eventweight*=ntmjdist.systematics_arrays_down[k][0]; eventweight_opp*=ntmjdist.systematics_arrays_down[k][0]
							elif template_name.find(ntmjdist.systematics[k]+'__up')!=-1 : 
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
					if self.sum_charge or ntmjdist.name.startswith('allchannels') or (ntmjdist.Q_l[0]>0 and template_name.find('plus')!=-1) or (ntmjdist.Q_l[0]<0 and template_name.find('minus')!=-1) :
						ntmjdist.__Fill__(j,ntmjdist.cstar[0],ntmjdist.x_F[0],ntmjdist.M[0],eventweight)
					if ntmjdist.addTwice[0]==1 and (self.sum_charge or ntmjdist.name.startswith('allchannels') or (ntmjdist.Q_l[0]<0 and template_name.find('plus')!=-1) or (ntmjdist.Q_l[0]>0 and template_name.find('minus')!=-1)) :
						ntmjdist.__Fill__(j,-1.0*ntmjdist.cstar[0],ntmjdist.x_F[0],ntmjdist.M[0],eventweight_opp)
			if not ntmjdist.name.startswith('allchannels') :
				#Save the templates to the auxiliary file
				self.f_aux.cd()
				returnhistos = []
				for histo in ntmjdist.all_histos :
					histo.Write()
				for i in range(len(ntmjdist.all_histos)/4) :
					returnhistos.append(convertTo1D(ntmjdist.all_histos[4*i]))
				#add to the list of new 1D templates
				self.histos+=returnhistos

	def make_plots(self) :
		#First make a list of all the channels in the file
		channel_names = []
		for dist in self.dists :
			if dist.name.endswith('DATA') :
				a = dist.name.split('__')
				channel_names.append(a[0])
		#Make lists of histogram stacks, residual plots, canvases, and pads
		x_stacks = []; y_stacks = []; z_stacks = []
		x_resids = []; y_resids = []; z_resids = []
		x_canvs  = []; y_canvs  = []; z_canvs  = []
		x_histo_pads = []; y_histo_pads = []; z_histo_pads = []
		x_resid_pads = []; y_resid_pads = []; z_resid_pads = []
		for channame in channel_names :
			#Projection histo stacks
			x_stacks.append(THStack(channame+'_x_stack',channame+' channel comparison plot, c* projection;;Events'))
			y_stacks.append(THStack(channame+'_y_stack',channame+' channel comparison plot, x_{F} projection;;Events'))
			z_stacks.append(THStack(channame+'_z_stack',channame+' channel comparison plot, M projection;;Events'))
			#Residual plots
			x_resids.append(TH1D(channame+'_x_residuals','; c*; data/MC',XBINS,XMIN,XMAX))
			y_resids.append(TH1D(channame+'_y_residuals','; |x_{F}|; data/MC',YBINS,YMIN,YMAX))
			z_resids.append(TH1D(channame+'_z_residuals','; M (GeV); data/MC',ZBINS,ZMIN,ZMAX))
			#Canvases
			x_canvs.append(TCanvas(channame+'_x_canvas',channame+'_x_canvas',900,900))
			y_canvs.append(TCanvas(channame+'_y_canvas',channame+'_y_canvas',900,900))
			z_canvs.append(TCanvas(channame+'_z_canvas',channame+'_z_canvas',900,900))
		#Set plot directories
		for i in range(len(channel_names)) :
			x_resids[i].SetDirectory(0); y_resids[i].SetDirectory(0); z_resids[i].SetDirectory(0)
		#build histogram stacks
		for i in range(len(channel_names)) :
			channame = channel_names[i]
			for dist in self.dists :
				if dist.name.startswith(channame) and not dist.name.endswith('DATA') :
					for j in range(len(dist.all_histo_names)) :
						endofname = dist.all_histo_names[j].split('__')[len(dist.all_histo_names[j].split('__'))-1].split('_')
						if 'up' not in endofname and 'down' not in endofname :
							if 'x' in endofname :
								x_histo = dist.all_histos[j]
								x_histo.SetMarkerStyle(21); x_histo.SetFillColor(dist.color); x_histo.SetMarkerColor(dist.color)
								x_stacks[i].Add(x_histo)
							if 'y' in endofname :
								y_histo = dist.all_histos[j]
								y_histo.SetMarkerStyle(21); y_histo.SetFillColor(dist.color); y_histo.SetMarkerColor(dist.color)
								y_stacks[i].Add(y_histo)
							if 'z' in endofname :
								z_histo = dist.all_histos[j]
								z_histo.SetMarkerStyle(21); z_histo.SetFillColor(dist.color); z_histo.SetMarkerColor(dist.color)
								z_stacks[i].Add(z_histo)
		#build residuals plots
		maxxdeviations = []; maxydeviations = []; maxzdeviations = []
		for i in range(len(channel_names)) :
			channame = channel_names[i]
			maxxdeviations.append(0.0); maxydeviations.append(0.0); maxzdeviations.append(0.0)
			for dist in self.dists :
				if dist.name.startswith(channame) and dist.name.endswith('DATA') :
					for j in range(1,XBINS+1) :
						if x_stacks[i].GetStack().Last().GetBinContent(j) != 0 :
							content = dist.all_histos[1].GetBinContent(j)/x_stacks[i].GetStack().Last().GetBinContent(j)
							x_resids[i].SetBinContent(j,content)
							if dist.all_histos[1].GetBinContent(j) != 0 :
								error = content*sqrt(1./dist.all_histos[1].GetBinContent(j) + 1./x_stacks[i].GetStack().Last().GetBinContent(j))
								x_resids[i].SetBinError(j,error)
								maxxdeviations[i] = max(maxxdeviations[i],max(abs(content+error-1.0),abs(content-error-1.0)))
					for j in range(1,YBINS+1) :
						if y_stacks[i].GetStack().Last().GetBinContent(j) != 0 :
							content = dist.all_histos[2].GetBinContent(j)/y_stacks[i].GetStack().Last().GetBinContent(j)
							y_resids[i].SetBinContent(j,content)
							if dist.all_histos[2].GetBinContent(j) != 0 :
								error = content*sqrt(1./dist.all_histos[2].GetBinContent(j) + 1./y_stacks[i].GetStack().Last().GetBinContent(j))
								y_resids[i].SetBinError(j,error)
								maxydeviations[i] = max(maxydeviations[i],max(abs(content+error-1.0),abs(content-error-1.0)))
					for j in range(1,ZBINS+1) :
						if z_stacks[i].GetStack().Last().GetBinContent(j) != 0 :
							content = dist.all_histos[3].GetBinContent(j)/z_stacks[i].GetStack().Last().GetBinContent(j)
							z_resids[i].SetBinContent(j,content)
							if dist.all_histos[3].GetBinContent(j) != 0 :
								error = content*sqrt(1./dist.all_histos[3].GetBinContent(j) + 1./z_stacks[i].GetStack().Last().GetBinContent(j))
								z_resids[i].SetBinError(j,error)
								maxzdeviations[i] = max(maxzdeviations[i],max(abs(content+error-1.0),abs(content-error-1.0)))
		#reset stack maxima
		for i in range(len(channel_names)) :
			for dist in self.dists :
				if dist.name.startswith(channel_names[i]) and dist.name.endswith('DATA') :
					xmaxdata = dist.all_histos[1].GetMaximum() 
					ymaxdata = dist.all_histos[2].GetMaximum() 
					zmaxdata = dist.all_histos[3].GetMaximum()
					x_stacks[i].SetMaximum(1.02*max(x_stacks[i].GetMaximum(),xmaxdata+sqrt(xmaxdata)))
					y_stacks[i].SetMaximum(1.02*max(y_stacks[i].GetMaximum(),ymaxdata+sqrt(ymaxdata)))
					z_stacks[i].SetMaximum(1.02*max(z_stacks[i].GetMaximum(),zmaxdata+sqrt(zmaxdata)))
		#Set histogram and residual plot properties
		for i in range(len(channel_names)) :
			x_resids[i].SetStats(0)
			x_resids[i].GetXaxis().SetLabelSize((0.05*0.72)/0.28); x_resids[i].GetXaxis().SetTitleOffset(0.8)
			x_resids[i].GetYaxis().SetLabelSize((0.05*0.72)/0.28); x_resids[i].GetYaxis().SetTitleOffset(0.4)
			x_resids[i].GetXaxis().SetTitleSize((0.72/0.28)*x_resids[i].GetXaxis().GetTitleSize())
			x_resids[i].GetYaxis().SetTitleSize((0.72/0.28)*x_resids[i].GetYaxis().GetTitleSize())
			maxx = 1.0+1.1*maxxdeviations[i]
			minx = 1.0-1.1*maxxdeviations[i]
			x_resids[i].GetYaxis().SetRangeUser(minx,maxx)
			x_resids[i].GetYaxis().SetNdivisions(503)
			y_resids[i].SetStats(0)
			y_resids[i].GetXaxis().SetLabelSize((0.05*0.72)/0.28); y_resids[i].GetXaxis().SetTitleOffset(0.8)
			y_resids[i].GetYaxis().SetLabelSize((0.05*0.72)/0.28); y_resids[i].GetYaxis().SetTitleOffset(0.4)
			y_resids[i].GetXaxis().SetTitleSize((0.72/0.28)*y_resids[i].GetXaxis().GetTitleSize())
			y_resids[i].GetYaxis().SetTitleSize((0.72/0.28)*y_resids[i].GetYaxis().GetTitleSize())
			maxy = 1.0+1.1*maxydeviations[i]
			miny = 1.0-1.1*maxydeviations[i]
			y_resids[i].GetYaxis().SetRangeUser(miny,maxy)
			y_resids[i].GetYaxis().SetNdivisions(503)
			z_resids[i].SetStats(0)
			z_resids[i].GetXaxis().SetLabelSize((0.05*0.72)/0.28); z_resids[i].GetXaxis().SetTitleOffset(0.8)
			z_resids[i].GetYaxis().SetLabelSize((0.05*0.72)/0.28); z_resids[i].GetYaxis().SetTitleOffset(0.4)
			z_resids[i].GetXaxis().SetTitleSize((0.72/0.28)*z_resids[i].GetXaxis().GetTitleSize())
			z_resids[i].GetYaxis().SetTitleSize((0.72/0.28)*z_resids[i].GetYaxis().GetTitleSize())
			maxz = 1.0+1.1*maxzdeviations[i]
			minz = 1.0-1.1*maxzdeviations[i]
			z_resids[i].GetYaxis().SetRangeUser(minz,maxz)
			z_resids[i].GetYaxis().SetNdivisions(503)
		#Build a legend
		leg = TLegend(0.62,0.67,0.9,0.9)
		for dist in self.dists :
			if dist.name.startswith(channel_names[0]) and not dist.name.endswith('DATA') :
				leg.AddEntry(dist.all_histos[1],dist.name.lstrip(channel_names[0]+'__'),'F')
			elif dist.name.startswith(channel_names[0]) and dist.name.endswith('DATA') :
				leg.AddEntry(dist.all_histos[1],dist.name.lstrip(channel_names[0]+'__'),'PE')
		#Build the lines that go at 1 on the residuals plots
		xline = TLine(XMIN,1.0,XMAX,1.0); xline.SetLineWidth(2); xline.SetLineStyle(2)
		yline = TLine(YMIN,1.0,YMAX,1.0); yline.SetLineWidth(2); yline.SetLineStyle(2)
		zline = TLine(ZMIN,1.0,ZMAX,1.0); zline.SetLineWidth(2); zline.SetLineStyle(2)
		#plot stacks with data overlaid and residuals
		for i in range(len(channel_names)) :
			channame = channel_names[i]
			for dist in self.dists :
				if dist.name.startswith(channame) and dist.name.endswith('DATA') :
					x_canvs[i].cd() 
					x_histo_pad=TPad(channame+'_x_histo_pad',channame+'_x_histo_pad',0,0.25,1,1)
					x_resid_pad=TPad(channame+'_x_residuals_pad',channame+'_x_residuals_pad',0,0,1.,0.25)
					x_histo_pad.SetCanvas(x_canvs[i]); x_resid_pad.SetCanvas(x_canvs[i])
					x_histo_pad.SetLeftMargin(0.16); x_histo_pad.SetRightMargin(0.05) 
					x_histo_pad.SetTopMargin(0.11);	 x_histo_pad.SetBottomMargin(0.02)
					x_histo_pad.SetBorderMode(0)
					x_resid_pad.SetLeftMargin(0.16); x_resid_pad.SetRightMargin(0.05)
					x_resid_pad.SetTopMargin(0.0);   x_resid_pad.SetBottomMargin(0.3)
					x_resid_pad.SetBorderMode(0)
					x_resid_pad.Draw(); x_histo_pad.Draw()
					x_histo_pad.cd(); 
					x_stacks[i].Draw(); dist.all_histos[1].Draw('SAME PE1X0'); x_stacks[i].GetXaxis().SetLabelOffset(999)
					leg.Draw()
					x_resid_pad.cd(); 
					x_resids[i].Draw('PE1X0'); xline.Draw()
					x_canvs[i].Update()
					self.f_aux.cd()
					x_canvs[i].Write()
					
					y_canvs[i].cd() 
					y_histo_pad=TPad(channame+'_y_histo_pad',channame+'_y_histo_pad',0,0.3,1,1)
					y_resid_pad=TPad(channame+'_y_residuals_pad',channame+'_y_residuals_pad',0,0,1.,0.3)
					y_histo_pad.SetCanvas(y_canvs[i]); y_resid_pad.SetCanvas(y_canvs[i])
					y_histo_pad.SetLeftMargin(0.16); y_histo_pad.SetRightMargin(0.05) 
					y_histo_pad.SetTopMargin(0.11);	 y_histo_pad.SetBottomMargin(0.02)
					y_histo_pad.SetBorderMode(0)
					y_resid_pad.SetLeftMargin(0.16); y_resid_pad.SetRightMargin(0.05)
					y_resid_pad.SetTopMargin(0.0);   y_resid_pad.SetBottomMargin(0.3)
					y_resid_pad.SetBorderMode(0)
					y_resid_pad.Draw(); y_histo_pad.Draw()
					y_histo_pad.cd(); 
					y_stacks[i].Draw(); dist.all_histos[2].Draw('SAME PE1X0'); y_stacks[i].GetXaxis().SetLabelOffset(999)
					leg.Draw()
					y_resid_pad.cd(); 
					y_resids[i].Draw('PE1X0'); yline.Draw()
					y_canvs[i].Update()
					self.f_aux.cd()
					y_canvs[i].Write()
					
					z_canvs[i].cd() 
					z_histo_pad=TPad(channame+'_z_histo_pad',channame+'_z_histo_pad',0,0.3,1,1)
					z_resid_pad=TPad(channame+'_z_residuals_pad',channame+'_z_residuals_pad',0,0,1.,0.3)
					z_histo_pad.SetCanvas(z_canvs[i]); z_resid_pad.SetCanvas(z_canvs[i])
					z_histo_pad.SetLeftMargin(0.16); z_histo_pad.SetRightMargin(0.05) 
					z_histo_pad.SetTopMargin(0.11);	 z_histo_pad.SetBottomMargin(0.02)
					z_histo_pad.SetBorderMode(0)
					z_resid_pad.SetLeftMargin(0.16); z_resid_pad.SetRightMargin(0.05)
					z_resid_pad.SetTopMargin(0.0);   z_resid_pad.SetBottomMargin(0.3)
					z_resid_pad.SetBorderMode(0)
					z_resid_pad.Draw(); z_histo_pad.Draw()
					z_histo_pad.cd(); 
					z_stacks[i].Draw(); dist.all_histos[3].Draw('SAME PE1X0'); z_stacks[i].GetXaxis().SetLabelOffset(999)
					leg.Draw()
					z_resid_pad.cd(); 
					z_resids[i].Draw('PE1X0'); zline.Draw()
					z_canvs[i].Update()
					self.f_aux.cd()
					z_canvs[i].Write()

	#__addAllDistributions__ sets up all of the final distributions depending on whether we want the charges summed
	def __addAllDistributions__(self) :
		lepprefixes = ['allchannels','mu','el']
		std_reweights = ['weight','sf_top_pT','sf_pileup']
		std_systematics = ['sf_lep_ID']
		for lepprefix in lepprefixes :
			chargeseps = ['']
			if not self.sum_charge and lepprefix != 'allchannels' :
				chargeseps = ['plus','minus']
			for chargesep in chargeseps :
				self.dists.append(distribution(lepprefix+chargesep+'__DATA',kBlack,'data distribution',None,None,None,None))
				self.dists.append(distribution(lepprefix+chargesep+'__fg0',kBlue,'0th gg (qg,q_{i}q_{j},etc.) distribution',None,std_reweights,std_systematics,'#scale#*(3.-#Rbck#-#Rntmj#)*(2.-#Rqqbar#)'))
				self.dists.append(distribution(lepprefix+chargesep+'__fqs0',kRed,'0th Symmetric q#bar{q} distribution',None,std_reweights,std_systematics,'#scale#*(3.-#Rbck#-#Rntmj#)*#Rqqbar#'))
				self.dists.append(distribution(lepprefix+chargesep+'__fqa0',kRed,'0th Antisymmetric q#bar{q} distribution','wqa0',std_reweights,std_systematics,'#scale#*(3.-#Rbck#-#Rntmj#)*#Rqqbar#*#Afb#'))
				self.dists.append(distribution(lepprefix+chargesep+'__fbck',kYellow,'background distribution',None,std_reweights,std_systematics,'#scale#*#Rbck#'))
				self.dists.append(distribution(lepprefix+chargesep+'__fntmj',kGreen,'NTMJ background distribution',None,std_reweights,std_systematics,'#scale#*#Rntmj#'))

	#__contributesToDist__ finds out whether the tree coming in should be added to the ith distribution
	def __contributesToDist__(self,ifd,i) :
		distname = self.dists[i].name
		if ((distname.startswith('mu') or distname.startswith('allchannels')) and ifd=='mudata' and (distname.endswith('DATA') or distname.endswith('fntmj'))) :
			return 1.0
		elif ((distname.startswith('el') or distname.startswith('allchannels')) and ifd=='eledata' and (distname.endswith('DATA') or distname.endswith('fntmj'))) :
			return 1.0
		elif (ifd == 'qq' and (distname.find('fqs')!=-1 or distname.find('fqa')!=-1)) :
			return LUMINOSITY
		elif (ifd == 'gg' and distname.find('fg')!=-1) or (ifd == 'bck' and distname.find('fbck')!=-1) :
			return LUMINOSITY
		elif (ifd == 'qq' or ifd == 'gg' or ifd == 'bck') and distname.find('fntmj')!=-1 :
			return -1.0*LUMINOSITY
		else :
			return 0.0

	#__del__ function
	def __del__(self) :
		#write to the file
		self.f_aux.cd()
		for dist in self.dists :
			if not dist.name.startswith('allchannels') :
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
	def __init__(self,name,color,formatted_name,sample_reweight,reweights,systematics,function) :
		self.name = name
		self.color = color
		self.formatted_name = formatted_name
		self.systematics = systematics
		self.function = function
		self.tree = TTree(self.name+'_tree',self.name+'_tree')
		self.tree.SetDirectory(0)
		self.all_branches = []
		self.cstar 				  = array('d',[100.]);  self.tree.Branch('cstar',self.cstar,'cstar/D'); 						self.all_branches.append(('cstar_scaled','cstar',self.cstar))
		self.x_F 				  = array('d',[100.]);  self.tree.Branch('x_F',self.x_F,'x_F/D'); 								self.all_branches.append(('x_F_scaled','x_F',self.x_F))
		self.M 					  = array('d',[-1.0]);  self.tree.Branch('M',self.M,'M/D'); 									self.all_branches.append(('M_scaled','M',self.M))
		self.hadt_M 			  = array('d',[-1.0]);  self.tree.Branch('hadt_M',self.hadt_M,'hadt_M/D'); 						self.all_branches.append(('scaled_hadt_M','hadt_M',self.hadt_M))
		self.hadt_tau32 		  = array('d',[-1.0]);  self.tree.Branch('hadt_tau32',self.hadt_tau32,'hadt_tau32/D'); 			self.all_branches.append(('hadt_tau32','hadt_tau32',self.hadt_tau32))
		self.Q_l 				  = array('i',[0]);     self.tree.Branch('Q_l',self.Q_l,'Q_l/I'); 								self.all_branches.append(('Q_l','Q_l',self.Q_l))
		self.addTwice 			  = array('I',[2]);     self.tree.Branch('addTwice',self.addTwice,'addTwice/i'); 				self.all_branches.append(('addTwice','addTwice',self.addTwice))
		self.contrib_weight 	  = array('d',[1.0]);   self.tree.Branch('contrib_weight',self.contrib_weight,'contrib_weight/D');
		self.sample_reweight 	  = array('d',[1.0]);   self.tree.Branch('sample_reweight',self.sample_reweight,'sample_reweight/D');
		self.sample_reweight_opp  = array('d',[1.0]);   self.tree.Branch('sample_reweight_opp',self.sample_reweight_opp,'sample_reweight_opp/D'); 
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
	def build_templates(self,aux_file,sumcharge,parfile) :
		#Nominal and due to fitting function parameters
		self.__organizeFunctionTemplates__(parfile)
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
						if template_name.find(self.systematics[j]+'__down')!=-1 :
							eventweight*=self.systematics_arrays_down[j][0]; eventweight_opp*=self.systematics_arrays_down[j][0]
						elif template_name.find(self.systematics[j]+'__up')!=-1 : 
							eventweight*=self.systematics_arrays_up[j][0]; eventweight_opp*=self.systematics_arrays_up[j][0]
						else :
							eventweight*=self.systematics_arrays[j][0]; eventweight_opp*=self.systematics_arrays[j][0]
				#apply the function of the fitting parameters
				eventweight*=self.factors[i]; eventweight_opp*=self.factors[i]
				#if we want to add the event twice half the weight
				if self.addTwice[0]==1 :
					eventweight*=0.5; eventweight_opp*=0.5
				#add to the template
				if sumcharge or self.name.startswith('allchannels') or (self.Q_l[0]>0 and template_name.find('plus')!=-1) or (self.Q_l[0]<0 and template_name.find('minus')!=-1) :
					self.__Fill__(i,self.cstar[0],self.x_F[0],self.M[0],eventweight)
				if self.addTwice[0]==1 and (sumcharge or self.name.startswith('allchannels') or (self.Q_l[0]<0 and template_name.find('plus')!=-1) or (self.Q_l[0]>0 and template_name.find('minus')!=-1)) :
					self.__Fill__(i,-1.0*self.cstar[0],self.x_F[0],self.M[0],eventweight_opp)
		
		#Save the templates to the auxiliary file
		aux_file.cd()
		returnhistos = []
		if not self.name.startswith('allchannels') :
			for histo in self.all_histos :
				histo.Write()
			for i in range(len(self.all_histos)/4) :
				returnhistos.append(convertTo1D(self.all_histos[4*i]))
		#return the list of new 1D templates
		return returnhistos

	#__organizeFunctionTemplates__ adds all of the necessary up/down templates for the fitting parameters
	#and returns the nominal factor for scaling the templates that are morphed for systematics
	def __organizeFunctionTemplates__(self,parfile_name) :
		#fit parameters and errors
		pars = []
		#Make appropriate templates from each line in the file
		parfile = open(parfile_name)
		for line in parfile :
			if line.startswith('#') :
				continue
			a = line.rstrip().split()
			if len(a) == 4 :
				[parname,initialvalue,downvalue,upvalue] = a
			elif len(a) == 5 :
				[parname,initialvalue,downvalue,upvalue,sigma] = a
			pars.append([parname,float(initialvalue),float(downvalue),float(upvalue)])
		#lists of new template prefactors and names
		self.factors = []; args = []; argispar = []
		if self.function != None :
			args = self.function.split('#')
		for j in range(len(args)) : #replace the parameters within the arguments
			argispar.append(False)
			for par in pars :
				if par[0] == args[j] :
					args[j] = par
					argispar[j] = True
					break
		#nominal distribution
		factorstring = ''
		for j in range(len(args)) :
			if not argispar[j] :
				factorstring+=args[j]
			else :
				factorstring+='('+str(args[j][1])+')'
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
			thisparname = ''
			for j in range(len(args)) :
				if not argispar[j] :
					factorstring_down+=args[j]; factorstring_up+=args[j]
				elif argispar[j] and countpars!=k :
					factorstring_down+='('+str(args[j][1])+')'; factorstring_up+='('+str(args[j][1])+')'
					countpars+=1
				elif argispar[j] and countpars==k :
					thisparname = args[j][0]
					factorstring_down+='('+str(args[j][2])+')'; factorstring_up+='('+str(args[j][3])+')'
					countpars+=1
			self.__addTemplate__(self.name+'__'+thisparname+'__'+'down',self.formatted_name+' '+thisparname+' down')
			self.__addTemplate__(self.name+'__'+thisparname+'__'+'up',self.formatted_name+' '+thisparname+' up')
			self.factors.append(eval(factorstring_down)); self.factors.append(eval(factorstring_up))
			print 'adding template with name %s and factor %s=%.4f'%(self.name+'__'+thisparname+'__'+'down',factorstring_down,self.factors[len(self.factors)-2])
			print 'adding template with name %s and factor %s=%.4f'%(self.name+'__'+thisparname+'__'+'up',factorstring_up,  self.factors[len(self.factors)-1])

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

	#__addTemplate__ function adds a 3D histogram and 1D projections to the file and to the lists
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

#convertTo1D takes a 3D distribution and makes it 1D for use with theta
def convertTo1D(original) :
	nBins = original.GetNbinsX()*original.GetNbinsY()*original.GetNbinsZ()
	newHisto = TH1F(original.GetName(),original.GetTitle(),nBins,0.,nBins-1.)
	newHisto.SetDirectory(0)
	for k in range(nBins) :
		newHisto.SetBinContent(k,original.GetBinContent(k))
	return newHisto