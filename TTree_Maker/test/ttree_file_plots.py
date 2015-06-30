#Takes in output files from the ttree maker and puts together the plots and stacks/lines them, then reoutputs them.
#imports
from ROOT import *
from optparse import OptionParser
from array import array
from math import *
import glob

#TDR Style
gROOT.Macro('rootlogon.C')

def newPlot(name,draw,cuts,nbins,low,high,title,options,weights,data) :
	histonames.append(name)
	drawstrings.append(draw)
	cutstrings.append(cuts)
	for i in range(len(filenames)) :
		histograms[i].append(TH1D(histonames[len(histonames)-1]+'_'+str(i),title+'_'+str(i),nbins,low,high))
	datahistograms.append(TH1D(histonames[len(histonames)-1]+'_data',title,nbins,low,high))
	histostacks.append(THStack(histonames[len(histonames)-1]+'_stack',title))
	optionstrings.append(options)
	weightstrings.append(weights)
	includedata.append(data)
	binningstrings.append('('+str(nbins)+','+str(low)+','+str(high)+')')

#OptionsParser
parser = OptionParser()
parser.add_option('--n', metavar='F', type='string', action='store',dest='outname',default='full_selection_plots',help='')
parser.add_option('--lep', metavar='F', type='string', action='store',dest='lep',default='mu',help='')
(options, args) = parser.parse_args()
leptype = 'none'
if 'mu' in options.lep.lower() :
	leptype = 'mu'
elif 'el' in options.lep.lower() :
	leptype = 'el'

#lists of filenames, fill colors, and pointers
filenames = []
fillcolors = []
sample_weights = []
data_filenames = []
#Single top
filenames.append('T_s'); 	 fillcolors.append(kMagenta)
filenames.append('T_t'); 	 fillcolors.append(kMagenta)
filenames.append('T_tW'); 	 fillcolors.append(kMagenta)
filenames.append('Tbar_s');  fillcolors.append(kMagenta)
filenames.append('Tbar_t');  fillcolors.append(kMagenta)
filenames.append('Tbar_tW'); fillcolors.append(kMagenta)
#DYnJets
filenames.append('DY1Jets'); fillcolors.append(kAzure-2)
filenames.append('DY2Jets'); fillcolors.append(kAzure-2)
filenames.append('DY3Jets'); fillcolors.append(kAzure-2)
filenames.append('DY4Jets'); fillcolors.append(kAzure-2)
#WnJets samples
filenames.append('W1Jets'); fillcolors.append(kGreen-3)
filenames.append('W2Jets'); fillcolors.append(kGreen-3)
filenames.append('W3Jets'); fillcolors.append(kGreen-3)
filenames.append('W4Jets'); fillcolors.append(kGreen-3)
#POWHEG TT
#dileptonic 
filenames.append('Powheg_dilep_TT'); 				 fillcolors.append(kRed-7)
filenames.append('Powheg_dilep_TT_SC'); 			 fillcolors.append(kRed-7)
filenames.append('Powheg_dilep_TT_Mtt_700_to_1000'); fillcolors.append(kRed-7)
filenames.append('Powheg_dilep_TT_Mtt_1000_to_Inf'); fillcolors.append(kRed-7)
#hadronic
filenames.append('Powheg_had_TT'); 				   fillcolors.append(kRed-7)
filenames.append('Powheg_had_TT_SC'); 			   fillcolors.append(kRed-7)
filenames.append('Powheg_had_TT_Mtt_700_to_1000'); fillcolors.append(kRed-7)
filenames.append('Powheg_had_TT_Mtt_1000_to_Inf'); fillcolors.append(kRed-7)
#semileptonic qq
filenames.append('Powheg_qq_semilep_TT'); 				  fillcolors.append(kRed+1)
filenames.append('Powheg_qq_semilep_TT_SC'); 			  fillcolors.append(kRed+1)
filenames.append('Powheg_qq_semilep_TT_Mtt_700_to_1000'); fillcolors.append(kRed+1)
filenames.append('Powheg_qq_semilep_TT_Mtt_1000_to_Inf'); fillcolors.append(kRed+1)
#semileptonic gg
filenames.append('Powheg_gg_semilep_TT'); 				  fillcolors.append(kRed+1)
filenames.append('Powheg_gg_semilep_TT_SC'); 			  fillcolors.append(kRed+1)
filenames.append('Powheg_gg_semilep_TT_Mtt_700_to_1000'); fillcolors.append(kRed+1)
filenames.append('Powheg_gg_semilep_TT_Mtt_1000_to_Inf'); fillcolors.append(kRed+1)
#data
data_filenames.append('Powheg_semilep_TT')
if leptype == 'mu' :
	data_filenames.append('SingleMu_Run2012A')
	data_filenames.append('SingleMu_Run2012B')
	data_filenames.append('SingleMu_Run2012C')
	data_filenames.append('SingleMu_Run2012D')
elif leptype == 'el' :
	data_filenames.append('SingleEl_Run2012A')
	data_filenames.append('SingleEl_Run2012B')
	data_filenames.append('SingleEl_Run2012C')
	data_filenames.append('SingleEl_Run2012D')

#lists of files
filelists = []
data_filelists = []
for filename in filenames :
	path = '../'+filename+'/'+filename+'_*_skim_tree.root'
	filelists.append(glob.glob(path))
for filename in data_filenames :
	path = '../'+filename+'/'+filename+'_*_skim_tree.root'
	data_filelists.append(glob.glob(path))
#lists of histogram names, temporary histograms for MC and data, histogram stacks, strings of what to draw, cutstrings, 
#options strings, weight strings, and renormalization and include data options
histonames = []
histostacks = []
histograms = []
for i in range(len(filenames)) :
	histograms.append([])
binningstrings = []
datahistograms = []
drawstrings = []
cutstrings = []
optionstrings = []
weightstrings = []
includedata = []

#Default Cut Strings
MARC_PRESELECTION = 'hadt_pt>300. && hadt_M>100. && lepW_pt>50.'
MARC_CUTS_MUONS = 'lepW_pt>50. && muon1_pt>40. && muon1_pt>ele1_pt && abs(muon1_eta)<2.4 && muon1_isLoose==1'
MARC_CUTS_MUONS+= ' && (muon1_relPt>25. || muon1_dR>0.5)'
MARC_CUTS_MUONS+= ' && lepb_M<50. && lept_M>140. && lept_M<250. && hadt_pt>400. && hadt_M>140. && hadt_M<250.'
MARC_CUTS_MUONS+= ' && hadt_tau32<0.55 && hadt_tau21>0.1'
MORE_SELECTIVE_MUONS = MARC_CUTS_MUONS+' && (muon2_isLoose!=1 || muon2_pt<40. || abs(muon2_eta)>2.4)'
MORE_SELECTIVE_MUONS+=' && (ele1_isLoose!=1 || ele1_pt<40. || abs(ele1_eta)>2.4)'
MORE_SELECTIVE_MUONS+=' && (ele2_isLoose!=1 || ele2_pt<40. || abs(ele2_eta)>2.4)'
MARC_CUTS_ELECTRONS = 'lepW_pt>50. && ele1_pt>40. && ele1_pt>muon1_pt && abs(ele1_eta)<2.4 && ele1_isLoose==1'
MARC_CUTS_ELECTRONS+= ' && (ele1_relPt>25. || ele1_dR>0.5)'
MARC_CUTS_ELECTRONS+= ' && lepb_M<50. && lept_M>140. && lept_M<250. && hadt_pt>400. && hadt_M>140. && hadt_M<250.'
MARC_CUTS_ELECTRONS+= ' && hadt_tau32<0.55 && hadt_tau21>0.1'
MORE_SELECTIVE_ELECTRONS = MARC_CUTS_ELECTRONS+' && (ele2_isLoose!=1 || ele2_pt<40. || abs(ele2_eta)>2.4)'
MORE_SELECTIVE_ELECTRONS+=' && (muon1_isLoose!=1 || muon1_pt<40. || abs(muon1_eta)>2.4)'
MORE_SELECTIVE_ELECTRONS+=' && (muon2_isLoose!=1 || muon2_pt<40. || abs(muon2_eta)>2.4)'
MARC_PLUS_TRIANGLE = MARC_CUTS_ELECTRONS+' && ele1_tri_el_val<ele1_tri_cut_val && ele1_tri_jet_val<ele1_tri_cut_val'
MORE_SELECTIVE_PLUS_TRIANGLE = MARC_PLUS_TRIANGLE+' && (ele2_isLoose!=1 || ele2_pt<25. || abs(ele2_eta)>2.4)'
MORE_SELECTIVE_PLUS_TRIANGLE+=' && (muon1_isLoose!=1 || muon1_pt<40. || abs(muon1_eta)>2.4)'
MORE_SELECTIVE_PLUS_TRIANGLE+=' && (muon2_isLoose!=1 || muon2_pt<40. || abs(muon2_eta)>2.4)'
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
#Weight strings
STD_WEIGHTS = 'weight*sf_pileup'

##############################################################################################################
############################						  PLOTS 					  ############################
##############################################################################################################

full_selection_cuts = MARC_PRESELECTION
if leptype == 'mu' :
	full_selection_cuts = MARC_CUTS_MUONS
elif leptype == 'el' :
	full_selection_cuts = MARC_CUTS_ELECTRONS

#M
newPlot('M','M',full_selection_cuts,25,0.,2500.,'M in Simulation and Data; M (GeV)','',STD_WEIGHTS,True)
#c*
newPlot('cstar','cstar',full_selection_cuts,20,-1.0,1.0,'c* in Simulation and Data; c*','',STD_WEIGHTS,True)
#xF
newPlot('x_F','abs(x_F)',full_selection_cuts,30,0.,0.6,'x_{F} in Simulation and Data; x_{F}','',STD_WEIGHTS,True)
#M (scaled)
newPlot('M_scaled','M_scaled',full_selection_cuts,25,0.,2500.,'scaled M in Simulation and Data; M (GeV)','',STD_WEIGHTS,True)
#c* (scaled)
newPlot('cstar_scaled','cstar_scaled',full_selection_cuts,20,-1.0,1.0,'scaled c* in Simulation and Data; c*','',STD_WEIGHTS,True)
#xF (scaled)
newPlot('x_F_scaled','abs(x_F_scaled)',full_selection_cuts,30,0.,0.6,'scaled x_{F} in Simulation and Data; x_{F}','',STD_WEIGHTS,True)

##leading muon pT
#newPlot('muon1_pt','muon1_pt',MARC_PRESELECTION,20,0.,100.,'leading muon p_{T}; p_{T} (GeV)','',STD_WEIGHTS,True)
##leading muon eta
#newPlot('muon1_eta','muon1_eta',MARC_PRESELECTION,35,-3.5,3.5,'leading muon #eta; #eta','',STD_WEIGHTS,True)
##second leading muon pT
#newPlot('muon2_pt','muon2_pt',MARC_PRESELECTION,20,0.,100.,'second leading muon p_{T}; p_{T} (GeV)','',STD_WEIGHTS,True)
##second leading muon eta
#newPlot('muon2_eta','muon2_eta',MARC_PRESELECTION,35,-3.5,3.5,'second leading muon #eta; #eta','',STD_WEIGHTS,True)
##leading electron pT
#newPlot('ele1_pt','ele1_pt',MARC_PRESELECTION,20,0.,100.,'leading electron p_{T}; p_{T} (GeV)','',STD_WEIGHTS,True)
##leading electron eta
#newPlot('ele1_eta','ele1_eta',MARC_PRESELECTION,35,-3.5,3.5,'leading electron #eta; #eta','',STD_WEIGHTS,True)
##second leading electron pT
#newPlot('ele2_pt','ele2_pt',MARC_PRESELECTION,20,0.,100.,'second leading electron p_{T}; p_{T} (GeV)','',STD_WEIGHTS,True)
##second leading electron eta
#newPlot('ele2_eta','ele2_eta',MARC_PRESELECTION,35,-3.5,3.5,'second leading electron #eta; #eta','',STD_WEIGHTS,True)
#leptonic W pT
#newPlot('lepW_pt','lepW_pt',MARC_PRESELECTION,20,0.,200.,'leptonic W candidate p_{T}; p_{T} (GeV)','',STD_WEIGHTS,True)
#leptonic b mass
#newPlot('lepb_M','lepb_M',MARC_PRESELECTION,20,0.,100.,'leptonic b candidate mass; M (GeV)','',STD_WEIGHTS,True)
#leptonic t mass
#newPlot('lept_M','lept_M',MARC_PRESELECTION,35,0.,350.,'leptonic top candidate mass; M (GeV)','',STD_WEIGHTS,True)
##hadronic t pT
#newPlot('hadt_pt','hadt_pt',ele_hadronic_pretag,50,250.,750.,'hadronic top candidate p_{T}; p_{T} (GeV)','',STD_WEIGHTS,True)
##hadronic t mass
#newPlot('hadt_M','hadt_M',ele_hadronic_pretag,35,0.,350.,'hadronic top candidate mass; M (GeV)','',STD_WEIGHTS,True)
##hadronic t tau32
#newPlot('hadt_tau32','hadt_tau32',ele_hadronic_pretag,20,0.,1.,'hadronic top candidate #tau_{32}; #tau_{32}','',STD_WEIGHTS,True)
##hadronic t tau21
#newPlot('hadt_tau21','hadt_tau21',MARC_PRESELECTION,20,0.,1.,'hadronic top candidate #tau_{21}; #tau_{21}','',STD_WEIGHTS,True)

##############################################################################################################
############################						  PLOTS 					  ############################
##############################################################################################################

#Chain up all data events
print 'chaining up data files'
data_chain = TChain('tree')
for data_file_list in data_filelists :
	for data_file in data_file_list :
		data_chain.Add(data_file)
print 'done'

#define the pruning cuts
prunestring = MARC_PRESELECTION
if leptype == 'el' :
	prunestring = MARC_PRESELECTION+' && ele1_pt>muon1_pt'
elif leptype == 'mu' :
	prunestring = MARC_PRESELECTION+' && muon1_pt>ele1_pt'

#Put the data into histograms
print 'Drawing events from data files into histograms...'
for j in range(len(histonames)) :
	print '	Doing '+histonames[j]+' ('+str(j+1)+' out of '+str(len(histonames))+')'
	tmp = datahistograms[j].Clone('tmp')
	data_chain.CopyTree(prunestring).Draw(drawstrings[j]+'>>tmp',cutstrings[j],'')
	datahistograms[j].Add(tmp)
	datahistograms[j].SetDirectory(0)
print 'Done'

#Chain up all the MC events, separated by file still though
MC_chains = []
for i in range(len(filelists)) :
	print 'Chaining up files for '+filenames[i]+' ('+str(i+1)+' out of '+str(len(filenames))+')'
	MC_chains.append(TChain('tree'))
	for filename in filelists[i] :
		MC_chains[i].Add(filename)
	print 'done'

#Build renormalization values
renorms = []
for j in range(len(histonames)) :
	renorms.append(0.0)
for i in range(len(MC_chains)) :
	print ( 'Drawing events from file '+filenames[i]+' into histograms to build renormalization values('
			+str(i+1)+' out of '+str(len(filenames))+')' )
	weight_array = array('d',[0.])
	MC_chains[i].SetBranchAddress('weight',weight_array)
	MC_chains[i].GetEntry(0)
	sample_weights.append(weight_array[0])
	for j in range(len(histonames)) :
		#Draw the plot I want from the tree
		print '	Drawing '+histonames[j]+' ('+str(j+1)+' out of '+str(len(histonames))+')'
		tmp = histograms[i][j].Clone('tmp')
		MC_chains[i].CopyTree(prunestring).Draw(drawstrings[j]+'>>tmp','('+weightstrings[j]+')*('+cutstrings[j]+')',optionstrings[j]+'')
		histograms[i][j].Add(tmp)
		#Set Colors and Fills
		histograms[i][j].SetFillColor(fillcolors[i])
		histograms[i][j].SetLineColor(fillcolors[i])
		histograms[i][j].SetMarkerStyle(21)
		histograms[i][j].SetDirectory(0)
		renorms[j] = renorms[j]+histograms[i][j].Integral()

#rescale MC histograms, add to stack
leg = TLegend(0.62,0.67,0.9,0.9)
for i in range(len(filenames)) :
	for j in range(len(histonames)) :
		#Renormalize to data
		if renorms[j]!=0. :
			histograms[i][j].Scale(datahistograms[j].Integral()/renorms[j])
		#add to histogram stack
		histostacks[j].Add(histograms[i][j])
	#Add to legend if need be
	if i==0 :
		leg.AddEntry(histograms[i][0],"Single Top","F")
	elif i==6 :
		leg.AddEntry(histograms[i][0],"Z/#gamma+Jets","F")
	elif i==10 :
		leg.AddEntry(histograms[i][0],"W+Jets","F")
	elif i==14 :
		leg.AddEntry(histograms[i][0],"Dileptonic/Hadronic t#bar{t}","F")
	elif i==len(filenames)-1 :
		leg.AddEntry(histograms[i][0],"Semileptonic t#bar{t}","F")
leg.AddEntry(datahistograms[0],"Data","LPE")

#plot on canvases
print 'plotting stacked plots'
canvs = []
for i in range(len(histonames)) :
	tmp_canv_name = histonames[i]+'_canv'
	tmp_canv = TCanvas(tmp_canv_name,tmp_canv_name,1200,900)
	tmp_canv.cd()
	datamax = datahistograms[i].GetMaximum()+sqrt(datahistograms[i].GetMaximum())
	histostacks[i].SetMaximum(max(1.02*datamax,1.02*histostacks[i].GetMaximum()))
	histostacks[i].Draw()
	if includedata[i] :
		datahistograms[i].SetMarkerStyle(20)
		datahistograms[i].DrawCopy("SAME PE1")
	leg.Draw()
	canvs.append(tmp_canv)

#Make line plots
unique_histos = []
histo_colors = []
unique_histo_maxes = []
#look through all the MC filetypes
print 'making line plots'
for i in range(len(filenames)) :
	#if we haven't yet seen a fill of this color add a new list of unique histograms and copy over
	if fillcolors[i] not in histo_colors :
		histo_colors.append(fillcolors[i])
		unique_histos.append([])
		for j in range(len(histonames)) :
			unique_histos[len(unique_histos)-1].append(histograms[i][j].Clone())
			unique_histos[len(unique_histos)-1][j].SetDirectory(0)
	#otherwise find the point in the list of unique histograms with histograms of this color and add these
	else : 
		for j in range(len(histonames)) :
			unique_histos[histo_colors.index(fillcolors[i])][j].Add(histograms[i][j].Clone())

#renormalize all the unique histograms to one and set the marker and line styles
for i in range(len(unique_histos)) :
	for j in range(len(unique_histos[i])) :
		if unique_histos[i][j].Integral() != 0. :
			unique_histos[i][j].Scale(1.0/unique_histos[i][j].Integral())
		unique_histos[i][j].SetFillStyle(0)
		unique_histos[i][j].SetLineWidth(3)
		unique_histos[i][j].SetLineStyle(1)
#renormalize the data histograms, too
for datahisto in datahistograms :
	if datahisto.Integral() != 0. :
		datahisto.Scale(1.0/datahisto.Integral())
	datahisto.SetFillStyle(0)
	datahisto.SetLineWidth(3)
	datahisto.SetFillColor(kBlack)
	datahisto.SetLineColor(kBlack)
	datahisto.SetLineStyle(1)
#find the max value on each unique histogram
for i in range(len(histonames)) :
	unique_histo_maxes.append(0.)
	for j in range(len(unique_histos)) :
		if 1.02*unique_histos[j][i].GetMaximum() > unique_histo_maxes[i] :
			unique_histo_maxes[i] = 1.02*unique_histos[j][i].GetMaximum()
#rebuild the legend
leg = TLegend(0.62,0.67,0.9,0.9)
for i in range(len(filenames)) :
	if i==0 :
		leg.AddEntry(histograms[i][0],"Single Top","F")
	elif i==6 :
		leg.AddEntry(histograms[i][0],"Z/#gamma+Jets","F")
	elif i==10 :
		leg.AddEntry(histograms[i][0],"W+Jets","F")
	elif i==14 :
		leg.AddEntry(histograms[i][0],"Dileptonic/Hadronic t#bar{t}","F")
	elif i==len(filenames)-1 :
		leg.AddEntry(histograms[i][0],"Semileptonic t#bar{t}","F")
leg.AddEntry(datahistograms[0],"Data","F")

#plot on canvases
other_canvs = []
for i in range(len(histonames)) :
	tmp_canv_name = histonames[i]+'_line_canv'
	tmp_canv = TCanvas(tmp_canv_name,tmp_canv_name,1200,900)
	tmp_canv.cd()
	other_canvs.append(tmp_canv)
	datamax = datahistograms[i].GetMaximum()
	unique_histos[0][i].SetMaximum(max(1.02*datamax,unique_histo_maxes[i]))
	unique_histos[0][i].Draw('')
	for j in range(1,len(unique_histos)) :
		unique_histos[j][i].Draw('SAME')
	if includedata[i] :
		datahistograms[i].Draw('SAME')
	leg.Draw()

#save in file
print 'saving plots'
# Make a new root file and save the histogram stacks to it
name = options.outname
if leptype == 'mu' :
	name += '_muons'
elif leptype == 'el' :
	name += '_electrons'
f = TFile(name+'.root', 'Recreate' )
f.cd()
for canv in canvs :
	canv.Write()
f.cd()
for canv in other_canvs :
	canv.Write()
f.Write()
f.Close()
