from ROOT import *
import glob
from array import array
from math import *
from distribution import distribution

#luminosity
LUMINOSITY = 19748.
LUMI_MIN_BIAS = 0.026

#histogram limits
XBINS = 20;		XMIN = -1.0;	XMAX = 1.0
YBINS = 6;		YMIN = 0.0;		YMAX = 0.6
ZBINS = 10;		ZMIN = 500.;	ZMAX = 2500.

#TDR Style
#gROOT.Macro('rootlogon.C')

##############################		   Template Group Class  		##############################

class template_group :
	#docstring
	"""template_group class"""
	
	#__init__function
	def __init__(self,outputname,parfilename) :
		#set the input variables
		self.step = parfilename.split('_')[0]
		self.f_nominal  = TFile(outputname+'nominal_'+self.step+'.root','recreate')
		self.f_simple_systematics = TFile(outputname+'simple_systematics_'+self.step+'.root','recreate')
		self.f_PDF_systematics = TFile(outputname+'PDF_systematics_'+self.step+'.root','recreate')
		self.f_JEC = TFile(outputname+'JEC_'+self.step+'.root','recreate')
		self.f_fit_parameters  = TFile(outputname+'fit_parameters_'+self.step+'.root','recreate')
		self.f_NTMJ = TFile(outputname+'NTMJ_'+self.step+'.root','recreate')
		self.f_aux = TFile(outputname+'aux_'+self.step+'.root','recreate')
		self.parfile = parfilename
		self.sum_charge = outputname.find('charge_summed')!=-1
		#final distributions
		print '	Initializing Distributions'
		self.dists = []
		self.__addAllDistributions__()
		print '	Done'

	#addTemplate function adds a new set of histograms for the given file, and sums the new template into the appropriate distribution
	def add_file_to_distributions(self,ttree_dir_path) :
		#chain up the ttree files
		jecwiggles = ['JES_up','JER_up','JES_down','JER_down']
		chains = [TChain('tree')]
		for i in range(len(jecwiggles)) :
			chains.append(TChain('tree'))
		filenamelist = glob.glob(ttree_dir_path+'/*_tree.root')
		for filename in filenamelist :
			iswiggle = False
			for i in range(len(jecwiggles)) :
				if ttree_dir_path.split('/')[len(ttree_dir_path.split('/'))-1].find('Run2012')!=-1 :
					chains[i+1].Add(filename)
				if filename.find(jecwiggles[i])!=-1 :
					chains[i+1].Add(filename)
					iswiggle = True
					break
			if not iswiggle :
				chains[0].Add(filename)
		for j in range(len(chains)) :
			chain = chains[j]
			#get all of the observables needed for making cuts
			muon1_pt 	  = array('d',[-1.0]);  chain.SetBranchAddress('muon1_pt',muon1_pt)
			ele1_pt 	  = array('d',[-1.0]);  chain.SetBranchAddress('ele1_pt',ele1_pt)
			lepW_pt 	  = array('d',[-1.0]);  chain.SetBranchAddress('scaled_lepW_pt',lepW_pt)
			hadt_pt 	  = array('d',[-1.0]);  chain.SetBranchAddress('scaled_hadt_pt',hadt_pt)
			muon1_eta 	  = array('d',[100.0]); chain.SetBranchAddress('muon1_eta',muon1_eta)
			muon1_isLoose = array('I',[2]);  	chain.SetBranchAddress('muon1_isLoose',muon1_isLoose)
			muon1_relPt   = array('d',[-1.0]);  chain.SetBranchAddress('muon1_relPt',muon1_relPt)
			muon1_dR 	  = array('d',[-1.0]);  chain.SetBranchAddress('muon1_dR',muon1_dR)
			mu_trigger 	  = array('I',[2]); 	chain.SetBranchAddress('mu_trigger',mu_trigger)
			ele1_eta 	  = array('d',[100.0]); chain.SetBranchAddress('ele1_eta',ele1_eta)
			ele1_isLoose  = array('I',[2]);  	chain.SetBranchAddress('ele1_isLoose',ele1_isLoose)
			ele1_relPt 	  = array('d',[-1.0]);  chain.SetBranchAddress('ele1_relPt',ele1_relPt)
			ele1_dR 	  = array('d',[-1.0]);  chain.SetBranchAddress('ele1_dR',ele1_dR)
			el_trigger 	  = array('I',[2]); 	chain.SetBranchAddress('el_trigger',el_trigger)
			lept_M 		  = array('d',[-1.0]);  chain.SetBranchAddress('scaled_lept_M',lept_M)
			hadt_M 		  = array('d',[-1.0]);  chain.SetBranchAddress('scaled_hadt_M',hadt_M)
			hadt_tau21 	  = array('d',[-1.0]);  chain.SetBranchAddress('hadt_tau21',hadt_tau21)
			hadt_tau32 	  = array('d',[-1.0]);  chain.SetBranchAddress('hadt_tau32',hadt_tau32)
			#build the list of distributions this tree will contribute to
			dists_to_add_to = []
			for i in range(len(self.dists)) :
				distname = self.dists[i].name
				#if events in this chain should be added to this distribution
				contribution = self.__contributesToDist__(ttree_dir_path,distname,j,jecwiggles)
				if contribution != 0. :
					dists_to_add_to.append(self.dists[i])
			nEntries = chain.GetEntries()
			print 'Adding trees from '+ttree_dir_path+' to distributions: '
			for i in range(len(dists_to_add_to)) :
				print '	'+dists_to_add_to[i].name
			#loop over events in tree
			for entry in range(nEntries) :
				chain.GetEntry(entry)
				cuts = []
				muon_preselection = muon1_pt[0]>ele1_pt[0] and lepW_pt[0]>50. and hadt_pt[0]>300. and hadt_M[0]>100.
				muon_kinematics   = muon1_pt[0]>40. and abs(muon1_eta[0])<2.4
				muon_ID = muon1_isLoose[0]==1
				muon_2D = muon1_relPt[0]>25. or muon1_dR[0]>0.5
				ele_preselection = ele1_pt[0]>muon1_pt[0] and lepW_pt[0]>50. and hadt_pt[0]>300. and hadt_M[0]>100.
				ele_kinematics = ele1_pt[0]>40. and abs(ele1_eta[0])<2.4
				ele_ID = ele1_isLoose[0]==1
				ele_2D = ele1_relPt[0]>25. or ele1_dR[0]>0.5
				lep_top_mass = lept_M[0]>140. and lept_M[0]<250.
				muon_full_leptonic = muon_preselection and muon_kinematics and muon_ID and muon_2D and lep_top_mass
				muon_hadronic_pretag = muon_full_leptonic and hadt_tau21[0]>0.1
				ele_full_leptonic = ele_preselection and ele_kinematics and ele_ID and ele_2D and lep_top_mass
				ele_hadronic_pretag = ele_full_leptonic and hadt_tau21[0]>0.1
				cuts.append(muon_hadronic_pretag or ele_hadronic_pretag)
				if cuts.count(False) > 0 :
					continue
				#copy over to the relevant distributions
				for i in range(len(dists_to_add_to)) :
					distname = dists_to_add_to[i].name
					#set all of the branch addresses to copy over the relevant parts of the tree
					for branch in dists_to_add_to[i].all_branches :
						if branch[0]!=None :
							chain.SetBranchAddress(branch[0],branch[2])
					chain.GetEntry(entry)
					#make appropriate lepton cuts
					if distname.startswith('mu') :
						cuts.append(muon_hadronic_pretag)
					elif distname.startswith('el') :
						cuts.append(ele_hadronic_pretag)
					if distname.find('plus')!=-1 :
						cuts.append(dists_to_add_to[i].Q_l[0]>0 or dists_to_add_to[i].addTwice[0]==1)
					elif distname.find('minus')!=-1 :
						cuts.append(dists_to_add_to[i].Q_l[0]<0 or dists_to_add_to[i].addTwice[0]==1)
					if cuts.count(False) > 0 :
						continue
					#change the event weight to add or subtract as appropriate
					contribution = self.__contributesToDist__(ttree_dir_path,distname,j,jecwiggles)
					dists_to_add_to[i].contrib_weight[0]=contribution
					dists_to_add_to[i].tree.Fill()

	#build_templates builds the template histograms from the distributions
	def build_templates(self) :
		for dist in self.dists :
			if dist.name.find('ntmj')!=-1 :
				continue
			print '	Building templates for distribution '+dist.name
			dist.build_templates()
			print '	Done'

	#build_NTMJ_templates automatically generates templates for all of the data-driven NTMJ background
	#automatically accounting for systematics
	def build_NTMJ_templates(self) :
		for dist in self.dists :
			if dist.name.find('ntmj')!=-1 :
				print '	Fixing NTMJ templates for distribution '+dist.name
				dist.fix_NTMJ_templates(self.f_aux)
				print '	Done'

	def write_to_files(self) :
		self.f_aux.cd()
		for dist in self.dists :
			if not dist.name.startswith('allchannels') :
				dist.tree.Write()
			for temp in dist.all_templates :
				temp.histo_3D.Write()
				temp.histo_x.Write()
				temp.histo_y.Write()
				temp.histo_z.Write()
		for dist in self.dists :
			if dist.name.find('ntmj') != -1 :
				self.f_NTMJ.cd() 
			elif dist.name.find('JER')!=-1 or dist.name.find('JES')!=-1 :
				self.f_JEC.cd()
			else :	
				self.f_simple_systematics.cd()
			for temp in dist.all_templates :
				if temp.name.find('up')==-1 and temp.name.find('down')==-1 :
					self.f_nominal.cd()
				elif temp.name.find('pdf_lambda')!=-1 :
					self.f_PDF_systematics.cd()	
				elif temp.name.find('par_')!=-1 :
					self.f_fit_parameters.cd()
				new1Dtemp = temp.convertTo1D()
				new1Dtemp.Write()
		self.f_nominal.Close()
		self.f_simple_systematics.Close()
		self.f_PDF_systematics.Close()
		self.f_JEC.Close()
		self.f_fit_parameters.Close()
		self.f_NTMJ.Close()
		self.f_aux.Close()

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
					for j in range(len(dist.all_templates)) :
						temp = dist.all_templates[j]
						endofname = temp.name.split('__')[len(temp.name.split('__'))-1]
						if endofname.find('up')==-1 and endofname.find('down')==-1 :
							temp = dist.all_templates[j]
							x_histo = temp.histo_x
							x_histo.SetMarkerStyle(21); x_histo.SetFillColor(dist.color); x_histo.SetMarkerColor(dist.color)
							x_stacks[i].Add(x_histo)
							y_histo = temp.histo_y
							y_histo.SetMarkerStyle(21); y_histo.SetFillColor(dist.color); y_histo.SetMarkerColor(dist.color)
							y_stacks[i].Add(y_histo)
							z_histo = temp.histo_z
							z_histo.SetMarkerStyle(21); z_histo.SetFillColor(dist.color); z_histo.SetMarkerColor(dist.color)
							z_stacks[i].Add(z_histo)
		#build residuals plots
		maxxdeviations = []; maxydeviations = []; maxzdeviations = []
		for i in range(len(channel_names)) :
			channame = channel_names[i]
			maxxdeviations.append(0.0); maxydeviations.append(0.0); maxzdeviations.append(0.0)
			for dist in self.dists :
				if dist.name.startswith(channame) and dist.name.endswith('DATA') :
					temp = dist.all_templates[0]
					for j in range(1,XBINS+1) :
						if x_stacks[i].GetStack().Last().GetBinContent(j) != 0 :
							content = temp.histo_x.GetBinContent(j)/x_stacks[i].GetStack().Last().GetBinContent(j)
							x_resids[i].SetBinContent(j,content)
							if temp.histo_x.GetBinContent(j) != 0 :
								error = content*sqrt(1./temp.histo_x.GetBinContent(j) + 1./x_stacks[i].GetStack().Last().GetBinContent(j))
								x_resids[i].SetBinError(j,error)
								maxxdeviations[i] = max(maxxdeviations[i],max(abs(content+error-1.0),abs(content-error-1.0)))
					for j in range(1,YBINS+1) :
						if y_stacks[i].GetStack().Last().GetBinContent(j) != 0 :
							content = temp.histo_y.GetBinContent(j)/y_stacks[i].GetStack().Last().GetBinContent(j)
							y_resids[i].SetBinContent(j,content)
							if temp.histo_y.GetBinContent(j) != 0 :
								error = content*sqrt(1./temp.histo_y.GetBinContent(j) + 1./y_stacks[i].GetStack().Last().GetBinContent(j))
								y_resids[i].SetBinError(j,error)
								maxydeviations[i] = max(maxydeviations[i],max(abs(content+error-1.0),abs(content-error-1.0)))
					for j in range(1,ZBINS+1) :
						if z_stacks[i].GetStack().Last().GetBinContent(j) != 0 :
							content = temp.histo_z.GetBinContent(j)/z_stacks[i].GetStack().Last().GetBinContent(j)
							z_resids[i].SetBinContent(j,content)
							if temp.histo_z.GetBinContent(j) != 0 and :
								error = content*sqrt(1./temp.histo_z.GetBinContent(j) + 1./z_stacks[i].GetStack().Last().GetBinContent(j))
								z_resids[i].SetBinError(j,error)
								maxzdeviations[i] = max(maxzdeviations[i],max(abs(content+error-1.0),abs(content-error-1.0)))
		#reset stack maxima
		for i in range(len(channel_names)) :
			for dist in self.dists :
				if dist.name.startswith(channel_names[i]) and dist.name.endswith('DATA') :
					temp = dist.all_templates[0]
					xmaxdata = temp.histo_x.GetMaximum() 
					ymaxdata = temp.histo_y.GetMaximum() 
					zmaxdata = temp.histo_z.GetMaximum()
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
				leg.AddEntry(dist.all_templates[0].histo_x,dist.name.lstrip(channel_names[0]+'__'),'F')
			elif dist.name.startswith(channel_names[0]) and dist.name.endswith('DATA') :
				leg.AddEntry(dist.all_templates[0].histo_x,dist.name.lstrip(channel_names[0]+'__'),'PE')
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
	#also handles the JEC corrections
	def __addAllDistributions__(self) :
#		lepprefixes = ['allchannels','mu','el']
		lepprefixes = ['el']
		for lepprefix in lepprefixes :
			chargeseps = ['']
			if not self.sum_charge and lepprefix != 'allchannels' :
				chargeseps = ['plus','minus']
			for chargesep in chargeseps :
				self.dists.append(distribution(self.parfile,lepprefix+chargesep+'__DATA',kBlack,'data distribution',None,None))
				self.dists.append(distribution(self.parfile,lepprefix+chargesep+'__fg0',kBlue,'0th gg (qg,q_{i}q_{j},etc.) distribution',None,'#scale#*(3.-#Rbck#-#Rntmj#)*(2.-#Rqqbar#)'))
				self.dists.append(distribution(self.parfile,lepprefix+chargesep+'__fqs0',kRed,'0th Symmetric q#bar{q} distribution',None,'#scale#*(3.-#Rbck#-#Rntmj#)*#Rqqbar#'))
				self.dists.append(distribution(self.parfile,lepprefix+chargesep+'__fqa0',kRed,'0th Antisymmetric q#bar{q} distribution','wqa0','#scale#*(3.-#Rbck#-#Rntmj#)*#Rqqbar#*#Afb#'))
				self.dists.append(distribution(self.parfile,lepprefix+chargesep+'__fbck',kYellow,'background distribution',None,'#scale#*#Rbck#'))
				self.dists.append(distribution(self.parfile,lepprefix+chargesep+'__fntmj',kGreen,'NTMJ background distribution',None,'#scale#*#Rntmj#'))
				if self.step == 'initial' or self.step == 'final' :
					JEC_wiggles = ['JES__up','JES__down','JER__up','JER__down']
					for JEC_wiggle in JEC_wiggles :
						self.dists.append(distribution(self.parfile,lepprefix+chargesep+'__fg0__'+JEC_wiggle,kBlue,'0th gg (qg,q_{i}q_{j},etc.) distribution',None,'#scale#*(3.-#Rbck#-#Rntmj#)*(2.-#Rqqbar#)'))
						self.dists.append(distribution(self.parfile,lepprefix+chargesep+'__fqs0__'+JEC_wiggle,kRed,'0th Symmetric q#bar{q} distribution',None,'#scale#*(3.-#Rbck#-#Rntmj#)*#Rqqbar#'))
						self.dists.append(distribution(self.parfile,lepprefix+chargesep+'__fqa0__'+JEC_wiggle,kRed,'0th Antisymmetric q#bar{q} distribution','wqa0','#scale#*(3.-#Rbck#-#Rntmj#)*#Rqqbar#*#Afb#'))
						self.dists.append(distribution(self.parfile,lepprefix+chargesep+'__fbck__'+JEC_wiggle,kYellow,'background distribution',None,'#scale#*#Rbck#'))
						self.dists.append(distribution(self.parfile,lepprefix+chargesep+'__fntmj__'+JEC_wiggle,kGreen,'NTMJ background distribution',None,'#scale#*#Rntmj#'))

	#__contributesToDist__ finds out whether the tree coming in should be added to the ith distribution
	def __contributesToDist__(self,ttree_path,distname,j,jecwiggles) :
		samplename = ttree_path.split('/')[len(ttree_path.split('/'))-1]
		bkg_names = ['dilep_TT','had_TT','T_s','T_t','T_tW','Tbar_s','Tbar_t','Tbar_tW']
		isbkg = False
		for bkg_name in bkg_names :
			if samplename.find(bkg_name)!=-1 :
				isbkg=True
				break
		iswig = distname.find('JES')!=-1 or distname.find('JER')!=-1
		wig = ''
		if j-1>=0 :
			wig = jecwiggles[j-1]
		if (iswig and (wig == '' or distname.find(wig.split('_')[0]+'__'+wig.split('_')[1])==-1)) or (not iswig and wig != '') :
			return 0.0
		if ((distname.startswith('mu') or distname.startswith('allchannels')) and samplename.find('SingleMu')!=-1 and (distname.endswith('DATA') or distname.find('fntmj')!=-1)) :
			return 1.0
		elif ((distname.startswith('el') or distname.startswith('allchannels')) and samplename.find('SingleEl')!=-1 and (distname.endswith('DATA') or distname.find('fntmj')!=-1)) :
			return 1.0
		elif (samplename.find('qq_semilep_TT')!=-1 and (distname.find('fqs')!=-1 or distname.find('fqa')!=-1)) :
			return LUMINOSITY
		elif samplename.find('gg_semilep_TT')!=-1 and distname.find('fg')!=-1 :
			return LUMINOSITY
		elif isbkg and distname.find('fbck')!=-1 :
			return LUMINOSITY
		elif samplename.find('SingleMu')==-1 and samplename.find('SingleEl')==-1 and distname.find('fntmj')!=-1 :
			return -1.0*LUMINOSITY
		else :
			return 0.0
