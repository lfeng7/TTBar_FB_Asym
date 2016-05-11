from ROOT import *
import glob
from array import array
from math import *
from distribution import distribution
from distribution import LM1_LOW, HM_HI
from template import template

#luminosity
LUMINOSITY = 19748.

#TDR Style
#gROOT.Macro('rootlogon.C')

##############################		   Template Group Class  		##############################

class template_group :
	#docstring
	"""template_group class"""
	
	#__init__function
	def __init__(self,outputname,parfilename,sumcharges,includepdf,includejec) :
		#set the input variables
		self.step = parfilename.split('_')[0]
		self.outputname = outputname
		self.f = TFile(outputname+'.root','recreate')
		self.f_aux = TFile(outputname+'_aux.root','recreate')
		self.parfile = parfilename
		self.postfit_histo_file = 'postfit_histos_'+outputname.split('_')[0]+'_'
		if self.step == 'refined' :
			self.postfit_histo_file+='initial.root'
		if self.step == 'final' :
			self.postfit_histo_file+='refined.root'
		self.sum_charge = sumcharges
		self.include_PDF = includepdf
		self.include_JEC = includejec
		#final distributions
		print '	Initializing Distributions'
		self.dists = []
		self.__addAllDistributions__()
		print '	Done'

	#addTemplate function adds a new set of histograms for the given file, and sums the new template into the appropriate distribution
	def add_file_to_distributions(self,ttree_file_path) :
		#get the tree from the ttree file
		jecwiggles = ['JES_up','JER_up','JES_down','JER_down']
		f = TFile(ttree_file_path)
		tree = f.Get('tree')
		#get all of the observables needed for making cuts
		original_branches = []
		muon1_pt 	  = array('d',[-1.0]);  original_branches.append(['muon1_pt',muon1_pt])
		ele1_pt 	  = array('d',[-1.0]);  original_branches.append(['ele1_pt',ele1_pt])
		lepW_pt 	  = array('d',[-1.0]);  original_branches.append(['scaled_lepW_pt',lepW_pt])
		hadt_pt 	  = array('d',[-1.0]);  original_branches.append(['scaled_hadt_pt',hadt_pt])
		muon1_eta 	  = array('d',[100.0]); original_branches.append(['muon1_eta',muon1_eta])
		muon1_isLoose = array('I',[2]);  	original_branches.append(['muon1_isLoose',muon1_isLoose])
		muon1_relPt   = array('d',[-1.0]);  original_branches.append(['muon1_relPt',muon1_relPt])
		muon1_dR 	  = array('d',[-1.0]);  original_branches.append(['muon1_dR',muon1_dR])
		mu_trigger 	  = array('I',[2]); 	original_branches.append(['mu_trigger',mu_trigger])
		ele1_eta 	  = array('d',[100.0]); original_branches.append(['ele1_eta',ele1_eta])
		ele1_isLoose  = array('I',[2]);  	original_branches.append(['ele1_isLoose',ele1_isLoose])
		ele1_relPt 	  = array('d',[-1.0]);  original_branches.append(['ele1_relPt',ele1_relPt])
		ele1_dR 	  = array('d',[-1.0]);  original_branches.append(['ele1_dR',ele1_dR])
		el_trigger 	  = array('I',[2]); 	original_branches.append(['el_trigger',el_trigger])
		lept_M 		  = array('d',[-1.0]);  original_branches.append(['scaled_lept_M',lept_M])
		hadt_M 		  = array('d',[-1.0]);  original_branches.append(['scaled_hadt_M',hadt_M])
		hadt_tau21 	  = array('d',[-1.0]);  original_branches.append(['hadt_tau21',hadt_tau21])
		hadt_tau32 	  = array('d',[-1.0]);  original_branches.append(['hadt_tau32',hadt_tau32])
		#build the list of distributions this tree will contribute to
		dists_to_add_to = []
		conttodist = self.__contributesToDist__
		for i in range(len(self.dists)) :
			distname = self.dists[i].name
			#if events in this tree should be added to this distribution
			contribution = conttodist(ttree_file_path,distname,jecwiggles)
			if contribution != 0. :
				dists_to_add_to.append(self.dists[i])
		if len(dists_to_add_to)==0 :
			return
		#loop over events in tree
		nEntries = tree.GetEntries()
		print 'Adding trees from '+ttree_file_path[ttree_file_path.find('total_ttree_files'):]+' to distributions: '
		for i in range(len(dists_to_add_to)) :
			print '	'+dists_to_add_to[i].name
		#Loop optimization
		tsetbranchaddress = tree.SetBranchAddress
		tgetentry = tree.GetEntry
		for entry in range(nEntries) :
			for branch in original_branches :
				tsetbranchaddress(branch[0],branch[1])
			tgetentry(entry)
			cuts = []
			muon_preselection = mu_trigger[0]==1 and muon1_pt[0]>ele1_pt[0] and hadt_pt[0]>300. and hadt_M[0]>50.
			muon_kinematics   = muon1_pt[0]>40. and abs(muon1_eta[0])<2.4
			muon_ID = muon1_isLoose[0]==1
			muon_2D = muon1_relPt[0]>25. or muon1_dR[0]>0.5
			ele_preselection = el_trigger[0]==1 and ele1_pt[0]>muon1_pt[0] and hadt_pt[0]>300. and hadt_M[0]>50.
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
						tsetbranchaddress(branch[0],branch[2])
				tgetentry(entry)
				#make appropriate lepton cuts
				xtra_cuts = []
				if distname.startswith('mu') :
					xtra_cuts.append(muon_hadronic_pretag)
				elif distname.startswith('el') :
					xtra_cuts.append(ele_hadronic_pretag)
				if distname.find('plus')!=-1 :
					xtra_cuts.append(dists_to_add_to[i].Q_l[0]>0 or dists_to_add_to[i].addTwice[0]==1)
				elif distname.find('minus')!=-1 :
					xtra_cuts.append(dists_to_add_to[i].Q_l[0]<0 or dists_to_add_to[i].addTwice[0]==1)
				if xtra_cuts.count(False) > 0 :
					continue
				#change the event weight to add or subtract as appropriate
				contribution = conttodist(ttree_file_path,distname,jecwiggles)
				dists_to_add_to[i].contrib_weight[0]=contribution
				#find the right tree to add to
				dists_to_add_to[i].filltree()

	#build_templates builds the template histograms from the distributions
	def build_templates(self) :
		for dist in self.dists :
			if dist.name.find('ntmj')!=-1 :
				continue
			print '	Building templates for distribution '+dist.name
			dist.build_templates(self.dists)
			print '	Done'

	#build_NTMJ_templates automatically generates templates for all of the data-driven NTMJ background
	#automatically accounting for systematics
	def build_NTMJ_templates(self) :
		#Build the event number and conversion function dictionaries
		all_event_counts = {}
		conv_functions = {}
		#For each type of leptons
		for leptype in self.lepprefixes :
			all_event_counts[leptype] = {'nominal':[0.,0.,0.,0.,0.,0.]}
			conv_functions[leptype] = {'nominal':[]}
		#And each type of template
		for ltype in all_event_counts :
			for dist in self.dists :
				if dist.name.split('__')[1].find('ntmj')!=-1 and dist.name.split('__')[0].find(ltype)!=-1 :
					for temp in dist.all_templates :
						tempnamesplit = temp.name.split('__')
						if len(tempnamesplit) > 2 :
							newtemptype = tempnamesplit[2]+'__'+tempnamesplit[3]
							if newtemptype!='fit__up' and newtemptype!='fit__down' and newtemptype not in all_event_counts[ltype].keys() :
								all_event_counts[ltype][newtemptype] = [0.,0.,0.,0.,0.,0.]
		#Set of cuts
		all_tree_cuts = ['lm1_ps','lm2_ps','hm_ps','lm1_fs','lm2_fs','hm_fs']
		#Fill the dictionary
		print '	Filling event numbers'
		for dist in self.dists :
			if dist.name.split('__')[1].find('ntmj')!=-1 :
				print '		Getting numbers for distribution '+dist.name
				dist.add_to_NTMJ_dictionary(self.dists,all_event_counts,all_tree_cuts)
		print '	Done'
		#Determine the conversion functions
		print '	Finding conversion functions'
		for ltype in all_event_counts :
			for ttype in all_event_counts[ltype] :
				print '		Doing conversion function for '+ltype+' leptons and template type '+ttype+''
				print '			NUMBERS OF EVENTS: ' #DEBUG
				for j in range(len(all_event_counts[ltype][ttype])) : #DEBUG
					if all_event_counts[ltype][ttype][j] <= 0. :
						all_event_counts[ltype][ttype][j] = 0.5
					print '				'+str(all_event_counts[ltype][ttype][j]) #DEBUG
				#Find the y values and errors for the points on the graph
				l1y = all_event_counts[ltype][ttype][0]/all_event_counts[ltype][ttype][3]
				l2y = all_event_counts[ltype][ttype][1]/all_event_counts[ltype][ttype][4]
				hy  = all_event_counts[ltype][ttype][2]/all_event_counts[ltype][ttype][5]
				l1ye = l1y*sqrt((1./all_event_counts[ltype][ttype][3])+(1./all_event_counts[ltype][ttype][0]))
				l2ye = l2y*sqrt((1./all_event_counts[ltype][ttype][4])+(1./all_event_counts[ltype][ttype][1]))
				hye  = hy*sqrt((1./all_event_counts[ltype][ttype][5])+(1./all_event_counts[ltype][ttype][2]))
				print '		l1y = '+str(l1y) #DEBUG
				print '		l2y = '+str(l2y) #DEBUG
				print '		hy = '+str(hy) #DEBUG
				print '		l1ye = '+str(l1ye) #DEBUG
				print '		l2ye = '+str(l2ye) #DEBUG
				print '		hye = '+str(hye) #DEBUG
				#Build the TGraph to fit with the conversion function
				n=3
				xs  = array('d',[80.,120.,275.])
				xes = array('d',[20.,20.,25.])
				ys  = array('d',[l1y,l2y,hy])
				yes = array('d',[l1ye,l2ye,hye])
				gr = TGraphErrors(n,xs,ys,xes,yes)
				#Define the conversion function
				conv_func = TF1('conv_func','[0]*x+[1]',LM1_LOW,HM_HI)
				#Fit the graph
				gr.Fit('conv_func')
				#get the fitter
				fitter = TVirtualFitter.GetFitter()
				#Build the nominal, slope up/down, and intercept up/down functions
				nom_func = TF1('nom_func','[0]*x+[1]',LM1_LOW,HM_HI) 
				nom_func.SetParameter(0,conv_func.GetParameter(0)) 
				nom_func.SetParameter(1,conv_func.GetParameter(1))
				fit_up_func = TF1('fit_up_func','[0]*x+[1]+sqrt([2]*[2]+2*x*[3]*[3]+x*x*[4]*[4])',LM1_LOW,HM_HI) 
				fit_up_func.SetParameter(0,conv_func.GetParameter(0)) 
				fit_up_func.SetParameter(1,conv_func.GetParameter(1))
				fit_up_func.SetParameter(2,conv_func.GetParError(1))
				fit_up_func.SetParameter(3,fitter.GetCovarianceMatrixElement(0,1))
				fit_up_func.SetParameter(4,conv_func.GetParError(0))
				fit_down_func = TF1('fit_down_func','[0]*x+[1]-sqrt([2]*[2]+2*x*[3]*[3]+x*x*[4]*[4])',LM1_LOW,HM_HI) 
				fit_down_func.SetParameter(0,conv_func.GetParameter(0)) 
				fit_down_func.SetParameter(1,conv_func.GetParameter(1))
				fit_down_func.SetParameter(2,conv_func.GetParError(1))
				fit_down_func.SetParameter(3,fitter.GetCovarianceMatrixElement(0,1))
				fit_down_func.SetParameter(4,conv_func.GetParError(0))
				conv_functions[ltype][ttype] = [nom_func,fit_up_func,fit_down_func]
				#Make the plot of the fit
				canv = TCanvas(ltype+'_'+ttype+'_conv_func_canv',ltype+' '+ttype+' conversion function canvas',900,900)
				gr.SetTitle(ltype+' '+ttype+' conversion function fit')
				gr.GetXaxis().SetTitle('hadronic top candidate mass (GeV)')
				gr.GetYaxis().SetTitle('conversion rate (N_{passed}/N_{failed})')
				gr.GetXaxis().SetRangeUser(LM1_LOW-20.,HM_HI+20.)
				gr.GetYaxis().SetRangeUser(0.,1.1*fit_up_func.Eval(300.))
				gr.SetMarkerStyle(21)
				gr.Draw('AP')
				nom_func.SetLineWidth(3)
				nom_func.SetLineColor(kBlue)
				nom_func.Draw('L SAME')
				fit_up_func.SetLineWidth(3)
				fit_up_func.SetLineStyle(2)
				fit_up_func.SetLineColor(kBlue)
				fit_up_func.Draw('L SAME')
				fit_down_func.SetLineWidth(3)
				fit_down_func.SetLineStyle(2)
				fit_down_func.SetLineColor(kBlue)
				fit_down_func.Draw('L SAME')
				leg = TLegend(0.62,0.67,0.9,0.9)
				leg.AddEntry(gr,'measured rates','PE')
				leg.AddEntry(nom_func,'linear fit','L')
				leg.AddEntry(fit_up_func,'fit error','L')
				leg.Draw()
				self.f_aux.cd()
				canv.Write()
		print '	Done'
		#Apply the conversion functions
		print '	Applying conversion functions'
		for dist in self.dists :
			if dist.name.split('__')[1].find('ntmj')!=-1 :
				dist.apply_conversion_functions(self.dists,conv_functions)
		print '	Done'

	def write_to_files(self) :
		#save all the new stuff in the auxiliary file
		self.f_aux.cd()
		for dist in self.dists :
			if not dist.name.startswith('allchannels') :
				dist.signal_tree.Write()
				dist.antitagged_tree.Write()
				dist.sidebands_tree.Write()
			for temp in dist.all_templates :
				temp.histo_3D.Write()
				temp.histo_x.Write()
				temp.histo_y.Write()
				temp.histo_z.Write()
		self.f_aux.Close()
		#if it's the first or last run just write everything we have to the output file
		if self.step == 'initial' or self.step == 'final' :
			self.f.cd()
			for dist in self.dists :
				if dist.name.startswith('allchannels') :
					continue
				for temp in dist.all_templates :
					print '	Adding newly created template '+temp.name
					new1Dtemp = temp.convertTo1D()
					new1Dtemp.Write()
			self.f.Close()
		#otherwise copy the templates from the original file, replacing with new ones where necessary
		else :
			self.f.cd()
			#make a list of all the new template names
			new_template_names = []
			for dist in self.dists :
				if dist.name.startswith('allchannels') :
					continue
				for temp in dist.all_templates :
					new_template_names.append(temp.name)
					#save the template to the new file
					print '	Adding newly created template '+temp.name
					new1Dtemp = temp.convertTo1D()
					new1Dtemp.Write()
			#open the file from the initial step and get its list of keys
			inifile = TFile(self.outputname.rstrip('_refined.root')+'_initial.root')
			inifilekeys = inifile.GetListOfKeys()
			self.f.cd()
			#for each key
			for key in inifilekeys :
				keyname = key.GetName()
				#if its name is reproduced in this run, skip it
				if keyname in new_template_names :
					continue
				#otherwise get its histogram and write it to this output file
				print '	Copying over previously created template '+keyname
				(inifile.Get(keyname)).Write()

		
	def make_plots(self) :
		#First make a list of all the channels in the file
		channel_names = ['allchannels']
		for dist in self.dists :
			if dist.name.endswith('DATA') :
				a = dist.name.split('__')
				if a[0] != 'allchannels' :
					channel_names.append(a[0])
		#Make lists of histogram stacks, MC uncertainty graphs, residual plots, canvases, and pads
		x_stacks = []; y_stacks = []; z_stacks = []
		x_err_hs = []; y_err_hs = []; z_err_hs = []
		x_resids = []; y_resids = []; z_resids = []
		x_canvs  = []; y_canvs  = []; z_canvs  = []
		x_histo_pads = []; y_histo_pads = []; z_histo_pads = []
		x_resid_pads = []; y_resid_pads = []; z_resid_pads = []
		for channame in channel_names :
			#Projection histo stacks
			x_stacks.append(THStack(channame+'_x_stack',channame+' channel comparison plot, c* projection;;Events'))
			y_stacks.append(THStack(channame+'_y_stack',channame+' channel comparison plot, x_{F} projection;;Events'))
			z_stacks.append(THStack(channame+'_z_stack',channame+' channel comparison plot, M projection;;Events'))
			#MC Uncertainty graphs
			x_err_hs.append(self.dists[0].all_templates[0].histo_x.Clone()); x_err_hs[len(x_err_hs)-1].Reset(); x_err_hs[len(x_err_hs)-1].SetFillColor(kBlack)
			y_err_hs.append(self.dists[0].all_templates[0].histo_y.Clone()); y_err_hs[len(y_err_hs)-1].Reset(); y_err_hs[len(y_err_hs)-1].SetFillColor(kBlack)
			z_err_hs.append(self.dists[0].all_templates[0].histo_z.Clone()); z_err_hs[len(z_err_hs)-1].Reset(); z_err_hs[len(z_err_hs)-1].SetFillColor(kBlack)
			#Residual plots
			x_resids.append(self.dists[0].all_templates[0].histo_x.Clone())
			x_resids[len(x_resids)-1].SetTitle(';c*;(Data-MC)/#sigma'); x_resids[len(x_resids)-1].SetFillColor(kMagenta+4)
			y_resids.append(self.dists[0].all_templates[0].histo_y.Clone())
			y_resids[len(x_resids)-1].SetTitle(';x_{F};(Data-MC)/#sigma'); y_resids[len(x_resids)-1].SetFillColor(kMagenta+4)
			z_resids.append(self.dists[0].all_templates[0].histo_z.Clone())
			z_resids[len(x_resids)-1].SetTitle(';M (GeV);(Data-MC)/#sigma'); z_resids[len(x_resids)-1].SetFillColor(kMagenta+4)
			#Canvases
			x_canvs.append(TCanvas(channame+'_x_canvas',channame+'_x_canvas',900,900))
			y_canvs.append(TCanvas(channame+'_y_canvas',channame+'_y_canvas',900,900))
			z_canvs.append(TCanvas(channame+'_z_canvas',channame+'_z_canvas',900,900))
		#Set plot directories
		for i in range(len(channel_names)) :
			x_err_hs[i].SetDirectory(0); y_err_hs[i].SetDirectory(0); z_err_hs[i].SetDirectory(0)
			x_resids[i].SetDirectory(0); y_resids[i].SetDirectory(0); z_resids[i].SetDirectory(0)
		nxbins = x_resids[0].GetNbinsX(); nybins = y_resids[0].GetNbinsX(); nzbins = z_resids[0].GetNbinsX()
		#build 3D histograms out of 1D post-fit histograms from theta; add to histogram stacks and MC uncertainty graphs
		postfithistofile = TFile(self.postfit_histo_file)
		disttypes = []
		disttypecolors = []
		for dist in self.dists :
			disttype = dist.name.split('__')[1]
			if not disttype.endswith('DATA') and not disttype in disttypes :
				disttypes.append(disttype)
				disttypecolors.append(dist.color)
		for k in range(len(disttypes)) :
			for i in range(len(channel_names)) :
				channame = channel_names[i]
				if channame == 'allchannels' :
					continue
				#Get histogram from post-fit file and make a new template
				newname = channame+'__'+disttypes[k]
				newtemp = template(newname+'__POSTFIT',newname+'__POSTFIT')
				new1Dhisto = postfithistofile.Get(newname)
				newtemp.make_from_1D_histo(new1Dhisto)
				#Set attributes and add to histogram stack
				x_histo = newtemp.histo_x
				x_histo.SetFillColor(disttypecolors[k]); x_histo.SetLineColor(disttypecolors[k]); x_histo.SetMarkerStyle(21); x_histo.SetMarkerColor(disttypecolors[k])
				x_stacks[i].Add(x_histo,'hist')
				x_stacks[0].Add(x_histo,'hist')
				#increment error values
				for j in range(x_histo.GetSize()) :
					if not x_histo.IsBinOverflow(j) and not x_histo.IsBinUnderflow(j) :
						x_err_hs[i].SetBinContent(j,x_err_hs[i].GetBinContent(j)+x_histo.GetBinContent(j))
						x_err_hs[i].SetBinError(j,x_err_hs[i].GetBinError(j)+x_histo.GetBinError(j)**2)
						x_err_hs[0].SetBinContent(j,x_err_hs[0].GetBinContent(j)+x_histo.GetBinContent(j))
						x_err_hs[0].SetBinError(j,x_err_hs[0].GetBinError(j)+x_histo.GetBinError(j)**2)
				#Repeat all of the above for y
				y_histo = newtemp.histo_y
				y_histo.SetFillColor(disttypecolors[k]); y_histo.SetLineColor(disttypecolors[k]); y_histo.SetMarkerStyle(21); y_histo.SetMarkerColor(disttypecolors[k])
				y_stacks[i].Add(y_histo,'hist')
				y_stacks[0].Add(y_histo,'hist')
				#increment error values
				for j in range(y_histo.GetSize()) :
					if not y_histo.IsBinOverflow(j) and not y_histo.IsBinUnderflow(j) :
						y_err_hs[i].SetBinContent(j,y_err_hs[i].GetBinContent(j)+y_histo.GetBinContent(j))
						y_err_hs[i].SetBinError(j,y_err_hs[i].GetBinError(j)+y_histo.GetBinError(j)**2)
						y_err_hs[0].SetBinContent(j,y_err_hs[0].GetBinContent(j)+y_histo.GetBinContent(j))
						y_err_hs[0].SetBinError(j,y_err_hs[0].GetBinError(j)+y_histo.GetBinError(j)**2)
				#Repeat all of the above for z
				z_histo = newtemp.histo_z
				z_histo.SetFillColor(disttypecolors[k]); z_histo.SetLineColor(disttypecolors[k]); z_histo.SetMarkerStyle(21); z_histo.SetMarkerColor(disttypecolors[k])
				z_stacks[i].Add(z_histo,'hist')
				z_stacks[0].Add(z_histo,'hist')
				#increment error values
				for j in range(z_histo.GetSize()) :
					if not z_histo.IsBinOverflow(j) and not z_histo.IsBinUnderflow(j) :
						z_err_hs[i].SetBinContent(j,z_err_hs[i].GetBinContent(j)+z_histo.GetBinContent(j))
						z_err_hs[i].SetBinError(j,z_err_hs[i].GetBinError(j)+z_histo.GetBinError(j)**2)
						z_err_hs[0].SetBinContent(j,z_err_hs[0].GetBinContent(j)+z_histo.GetBinContent(j))
						z_err_hs[0].SetBinError(j,z_err_hs[0].GetBinError(j)+z_histo.GetBinError(j)**2)
		#Take the root of all of the final MC uncertainties
		for i in range(len(channel_names)) :
			for j in range(nxbins) :
				x_err_hs[i].SetBinError(j,sqrt(x_err_hs[i].GetBinError(j)))
			for j in range(nybins) :
				y_err_hs[i].SetBinError(j,sqrt(y_err_hs[i].GetBinError(j)))
			for j in range(nzbins) :
				z_err_hs[i].SetBinError(j,sqrt(z_err_hs[i].GetBinError(j)))
		#build residuals plots
		maxxdeviations = []; maxydeviations = []; maxzdeviations = []
		for i in range(len(channel_names)) :
			channame = channel_names[i]
			maxxdeviations.append(0.0); maxydeviations.append(0.0); maxzdeviations.append(0.0)
			for dist in self.dists :
				if dist.name.startswith(channame) and dist.name.endswith('DATA') :
					for temp in dist.all_templates :
						if len(temp.name.split('__'))>2 :
							continue
						for j in range(temp.histo_x.GetSize()) :
							if temp.histo_x.IsBinUnderflow(j) or temp.histo_x.IsBinOverflow(j) :
								continue
							if x_err_hs[i].GetBinContent(j)!=0. and temp.histo_x.GetBinContent(j)!=0. :
								delta = temp.histo_x.GetBinContent(j) - x_err_hs[i].GetBinContent(j)
								sigma = sqrt(temp.histo_x.GetBinError(j)**2 + x_err_hs[i].GetBinError(j)**2)
								content = delta/sigma
#								print '		%s XHISTO BIN %d = (%.2f - %.2f)/sqrt(%.2f^2 + %.2f^2) = %.2f/%.2f = %.2f'%(channame,j,temp.histo_x.GetBinContent(j),x_err_hs[i].GetBinContent(j),temp.histo_x.GetBinError(j),x_err_hs[i].GetBinError(j),delta,sigma,content) #DEBUG
								x_resids[i].SetBinContent(j,content)
								maxxdeviations[i] = max(maxxdeviations[i],abs(content))
						for j in range(temp.histo_y.GetSize()) :
							if temp.histo_y.IsBinUnderflow(j) or temp.histo_y.IsBinOverflow(j) :
								continue
							if y_err_hs[i].GetBinContent(j)!=0. and temp.histo_y.GetBinContent(j)!=0. :
								delta = temp.histo_y.GetBinContent(j) - y_err_hs[i].GetBinContent(j)
								sigma = sqrt(temp.histo_y.GetBinError(j)**2 + y_err_hs[i].GetBinError(j)**2)
								content = delta/sigma
#								print '		%s YHISTO BIN %d = (%.2f - %.2f)/sqrt(%.2f^2 + %.2f^2) = %.2f/%.2f = %.2f'%(channame,j,temp.histo_y.GetBinContent(j),y_err_hs[i].GetBinContent(j),temp.histo_y.GetBinError(j),y_err_hs[i].GetBinError(j),delta,sigma,content) #DEBUG
								y_resids[i].SetBinContent(j,content)
								maxydeviations[i] = max(maxydeviations[i],abs(content))
						for j in range(temp.histo_z.GetSize()) :
							if temp.histo_z.IsBinUnderflow(j) or temp.histo_z.IsBinOverflow(j) :
								continue
							if z_err_hs[i].GetBinContent(j)!=0. and temp.histo_z.GetBinContent(j)!=0. :
								delta = temp.histo_z.GetBinContent(j) - z_err_hs[i].GetBinContent(j)
								sigma = sqrt(temp.histo_z.GetBinError(j)**2 + z_err_hs[i].GetBinError(j)**2)
								content = delta/sigma
#								print '		%s ZHISTO BIN %d = (%.2f - %.2f)/sqrt(%.2f^2 + %.2f^2) = %.2f/%.2f = %.2f'%(channame,j,temp.histo_z.GetBinContent(j),z_err_hs[i].GetBinContent(j),temp.histo_z.GetBinError(j),z_err_hs[i].GetBinError(j),delta,sigma,content) #DEBUG
								z_resids[i].SetBinContent(j,content)
								maxzdeviations[i] = max(maxzdeviations[i],abs(content))
		#reset stack maxima
		for i in range(len(channel_names)) :
			for dist in self.dists :
				if dist.name.startswith(channel_names[i]) and dist.name.endswith('DATA') :
					for temp in dist.all_templates :
						if len(temp.name.split('__'))>2 :
							continue
						xmaxdata = temp.histo_x.GetMaximum() 
						ymaxdata = temp.histo_y.GetMaximum() 
						zmaxdata = temp.histo_z.GetMaximum()
						x_stacks[i].SetMaximum(1.02*max(x_stacks[i].GetMaximum(),xmaxdata+sqrt(xmaxdata)))
						y_stacks[i].SetMaximum(1.02*max(y_stacks[i].GetMaximum(),ymaxdata+sqrt(ymaxdata)))
						z_stacks[i].SetMaximum(1.02*max(z_stacks[i].GetMaximum(),zmaxdata+sqrt(zmaxdata)))
		#Set histogram, MC error graph, and residual plot properties
		for i in range(len(channel_names)) :
			x_resids[i].SetStats(0)
			x_resids[i].GetXaxis().SetLabelSize((0.05*0.72)/0.28); x_resids[i].GetXaxis().SetTitleOffset(0.8)
			x_resids[i].GetYaxis().SetLabelSize((0.05*0.72)/0.28); x_resids[i].GetYaxis().SetTitleOffset(0.4)
			x_resids[i].GetXaxis().SetTitleSize((0.72/0.28)*x_resids[i].GetXaxis().GetTitleSize())
			x_resids[i].GetYaxis().SetTitleSize((0.72/0.28)*x_resids[i].GetYaxis().GetTitleSize())
			maxx = 0.1+ceil(maxxdeviations[i])
			minx = -0.1-ceil(maxxdeviations[i])
			x_resids[i].GetYaxis().SetRangeUser(minx,maxx)
			x_resids[i].GetYaxis().SetNdivisions(503)
			x_resids[i].SetMarkerStyle(20)
			x_err_hs[i].SetFillStyle(3005)
			y_resids[i].SetStats(0)
			y_resids[i].GetXaxis().SetLabelSize((0.05*0.72)/0.28); y_resids[i].GetXaxis().SetTitleOffset(0.8)
			y_resids[i].GetYaxis().SetLabelSize((0.05*0.72)/0.28); y_resids[i].GetYaxis().SetTitleOffset(0.4)
			y_resids[i].GetXaxis().SetTitleSize((0.72/0.28)*y_resids[i].GetXaxis().GetTitleSize())
			y_resids[i].GetYaxis().SetTitleSize((0.72/0.28)*y_resids[i].GetYaxis().GetTitleSize())
			maxy = 0.1+ceil(maxydeviations[i])
			miny = -0.1-ceil(maxydeviations[i])
			y_resids[i].GetYaxis().SetRangeUser(miny,maxy)
			y_resids[i].GetYaxis().SetNdivisions(503)
			y_resids[i].SetMarkerStyle(20)
			y_err_hs[i].SetFillStyle(3005)
			z_resids[i].SetStats(0)
			z_resids[i].GetXaxis().SetLabelSize((0.05*0.72)/0.28); z_resids[i].GetXaxis().SetTitleOffset(0.8)
			z_resids[i].GetYaxis().SetLabelSize((0.05*0.72)/0.28); z_resids[i].GetYaxis().SetTitleOffset(0.4)
			z_resids[i].GetXaxis().SetTitleSize((0.72/0.28)*z_resids[i].GetXaxis().GetTitleSize())
			z_resids[i].GetYaxis().SetTitleSize((0.72/0.28)*z_resids[i].GetYaxis().GetTitleSize())
			maxz = 0.1+ceil(maxzdeviations[i])
			minz = -0.1-ceil(maxzdeviations[i])
			z_resids[i].GetYaxis().SetRangeUser(minz,maxz)
			z_resids[i].GetYaxis().SetNdivisions(503)
			z_resids[i].SetMarkerStyle(20)
			z_err_hs[i].SetFillStyle(3005)
		#Build a legend
		leg = TLegend(0.62,0.67,0.9,0.9)
		for dist in self.dists :
			if dist.name.startswith(channel_names[0]) and dist.name.find('__up')==-1 and dist.name.find('__down')==-1 and not dist.name.endswith('DATA') :
				dist.all_templates[0].histo_x.SetFillColor(dist.color)
				leg.AddEntry(dist.all_templates[0].histo_x,dist.name.lstrip(channel_names[0]+'__'),'F')
			elif dist.name.startswith(channel_names[0]) and dist.name.endswith('DATA') :
				leg.AddEntry(dist.all_templates[0].histo_x,dist.name.lstrip(channel_names[0]+'__'),'PE')
		#plot stacks with data overlaid and residuals
		plotfile = TFile('comparison_plots_'+self.step+'_step.root','recreate')
		for i in range(len(channel_names)) :
			channame = channel_names[i]
			for dist in self.dists :
				if dist.name.startswith(channame) and dist.name.endswith('DATA') :
					for temp in dist.all_templates :
						if len(temp.name.split('__'))>2 :
							continue
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
						temp.histo_x.SetMarkerStyle(20)
						x_stacks[i].Draw(); temp.histo_x.Draw('SAME PE1X0'); x_err_hs[i].Draw('SAME E2'); x_stacks[i].GetXaxis().SetLabelOffset(999)
						leg.Draw()
						x_resid_pad.cd(); 
						x_resids[i].Draw('B')
						x_canvs[i].Update()
						plotfile.cd()
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
						temp.histo_y.SetMarkerStyle(20)
						y_stacks[i].Draw(); temp.histo_y.Draw('SAME PE1X0'); y_err_hs[i].Draw('SAME E2'); y_stacks[i].GetXaxis().SetLabelOffset(999)
						leg.Draw()
						y_resid_pad.cd(); 
						y_resids[i].Draw('B')
						y_canvs[i].Update()
						plotfile.cd()
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
						temp.histo_z.SetMarkerStyle(20)
						z_stacks[i].Draw(); temp.histo_z.Draw('SAME PE1X0'); z_err_hs[i].Draw('SAME E2'); z_stacks[i].GetXaxis().SetLabelOffset(999)
						leg.Draw()
						z_resid_pad.cd(); 
						z_resids[i].Draw('B')
						z_canvs[i].Update()
						plotfile.cd()
						z_canvs[i].Write()

	#__addAllDistributions__ sets up all of the final distributions depending on whether we want the charges summed
	#also handles the JEC corrections
	def __addAllDistributions__(self) :
		PREFAC_1 = '(#NTOT#-#NBCK#*#Rbck#-#NNTMJ#*#Rntmj#)*(1./#NTTBAR#)'
		PREFAC_2 = '(#NTTBAR#-#NQQBAR#*#Rqqbar#)*(1./(#NTTBAR#-#NQQBAR#))'
		FGG  = '(1. + #mu#*(1.-#mu#)*(#NG1#/(#NTTBAR#-#NQQBAR#))'
		FGG += '+ (#mu#*#mu#+#d#*#d#)*(1.+#mu#)*(#NG2#/(#NTTBAR#-#NQQBAR#))'
		FGG += '+ (#mu#*#mu#+#d#*#d#)*(1.-5.*#mu#)*(#NG3#/(#NTTBAR#-#NQQBAR#))'
		FGG += '+ (#mu#*#mu#+#d#*#d#)*(#mu#*#mu#+#d#*#d#)*(#NG4#/(#NTTBAR#-#NQQBAR#)))'
		FQQ  = '(1. + (2.*#mu#+#mu#*#mu#-#d#*#d#)*(#NQ1#/#NQQBAR#) + (#mu#*#mu#+#d#*#d#)*(#NQ2#/#NQQBAR#))'
		fbck_func  = '#scale#*#Rbck#'
		fntmj_func = '(1.)*#scale#*#Rntmj#'
		fgg_func   = '#scale#*'+PREFAC_1+'*'+PREFAC_2+'*(1./'+FGG+')*(1.+ #mu#*(1.-#mu#)*#wg1#'
		fgg_func  += '+ (#mu#*#mu#+#d#*#d#)*(1.+#mu#)*#wg2# + (#mu#*#mu#+#d#*#d#)*(1.-5.*#mu#)*#wg3#'
		fgg_func  += '+ (#mu#*#mu#+#d#*#d#)*(#mu#*#mu#+#d#*#d#)*#wg4#)'
		fqq_func   = '#scale#*'+PREFAC_1+'*(1./'+FQQ+')*#Rqqbar#*(1.+#Afb#*#wqa0# + (2.*#mu#+#mu#*#mu#-#d#*#d#)*(#wqs1#+#Afb#*#wqa1#)'
		fqq_func  += '+ (#mu#*#mu#+#d#*#d#)*(#wqs2#+#Afb#*#wqa2#))'
		fntmj_func = '#scale#*#Rntmj#*(1.+0.*(#Rbck#+#Rqqbar#+#Afb#))'
		self.lepprefixes = []
		if self.step == 'final' :
			self.lepprefixes.append('allchannels')
		self.lepprefixes.append('mu')
		self.lepprefixes.append('el')
		for lepprefix in self.lepprefixes :
			chargeseps = ['']
			if not self.sum_charge and lepprefix != 'allchannels' :
				chargeseps = ['plus','minus']
			for chargesep in chargeseps :
				self.dists.append(distribution(self.parfile,self.include_PDF,lepprefix+chargesep+'__DATA',kBlack,'data distribution',None))
				self.dists.append(distribution(self.parfile,self.include_PDF,lepprefix+chargesep+'__fgg',kBlue,'gg (qg,q_{i}q_{j},etc.) distribution',fgg_func))
				self.dists.append(distribution(self.parfile,self.include_PDF,lepprefix+chargesep+'__fqq',kRed+2,'0th Symmetric q#bar{q} distribution',fqq_func))
				self.dists.append(distribution(self.parfile,self.include_PDF,lepprefix+chargesep+'__fbck',kYellow,'background distribution',fbck_func))
				self.dists.append(distribution(self.parfile,self.include_PDF,lepprefix+chargesep+'__fntmj',kGreen,'NTMJ background distribution',fntmj_func))
				if self.include_JEC and self.step == 'initial' :
					JEC_wiggles = ['JES__up','JES__down','JER__up','JER__down']
					for JEC_wiggle in JEC_wiggles :
						self.dists.append(distribution(self.parfile,self.include_PDF,lepprefix+chargesep+'__fgg__'+JEC_wiggle,kBlue,JEC_wiggle+'gg (qg,q_{i}q_{j},etc.) distribution',fgg_func))
						self.dists.append(distribution(self.parfile,self.include_PDF,lepprefix+chargesep+'__fqq__'+JEC_wiggle,kRed+2,JEC_wiggle+'0th Symmetric q#bar{q} distribution',fqq_func))
						self.dists.append(distribution(self.parfile,self.include_PDF,lepprefix+chargesep+'__fbck__'+JEC_wiggle,kYellow,JEC_wiggle+'background distribution',fbck_func))
						self.dists.append(distribution(self.parfile,self.include_PDF,lepprefix+chargesep+'__fntmj__'+JEC_wiggle,kGreen,JEC_wiggle+'NTMJ background distribution',fntmj_func))

	#__contributesToDist__ finds out whether the tree coming in should be added to the ith distribution
	def __contributesToDist__(self,ttree_path,distname,jecwiggles) :
		samplename = ttree_path.split('/')[len(ttree_path.split('/'))-1][:ttree_path.split('/')[len(ttree_path.split('/'))-1].find('skim_all.root')]
		bkg_names = ['dilep_TT','had_TT','T_s','T_t','T_tW','Tbar_s','Tbar_t','Tbar_tW']
		isbkg = False
		for bkg_name in bkg_names :
			if samplename.find(bkg_name)!=-1 :
				isbkg=True
				break
		iswig = len(distname.split('__'))>2 and (distname.split('__')[2].find('JES')!=-1 or distname.split('__')[2].find('JER')!=-1)
		wig = ''
		for j in range(len(jecwiggles)) :
			if samplename.find(jecwiggles[j])!=-1 :
				wig = jecwiggles[j]
				break
		if wig != '' and not self.include_JEC :
			return 0.0
		disttype = distname.split('__')[1]
		if (iswig and (wig == '' or distname.find(wig.split('_')[0]+'__'+wig.split('_')[1])==-1) and samplename.find('Run2012')==-1) or (not iswig and wig != '') :
			return 0.0
		if ((distname.startswith('mu') or distname.startswith('allchannels')) and samplename.find('SingleMu')!=-1 and (distname.endswith('DATA') or distname.find('fntmj')!=-1)) :
			return 1.0
		elif ((distname.startswith('el') or distname.startswith('allchannels')) and samplename.find('SingleEl')!=-1 and (distname.endswith('DATA') or distname.find('fntmj')!=-1)) :
			return 1.0
		elif samplename.find('qq_semilep_TT')!=-1 and disttype.find('fqq')!=-1 :
			return LUMINOSITY
		elif samplename.find('gg_semilep_TT')!=-1 and disttype.find('fgg')!=-1 :
			return LUMINOSITY
		elif isbkg and disttype.find('fbck')!=-1 :
			return LUMINOSITY
		elif samplename.find('SingleMu')==-1 and samplename.find('SingleEl')==-1 and disttype.find('fntmj')!=-1 :
			return -1.0*LUMINOSITY
		else :
			return 0.0