#Takes in output files from the ttree maker and puts together the plots and stacks/lines them, then reoutputs them.
#imports
from ROOT import *
from optparse import OptionParser
from array import array
from math import *

#TDR Style
gROOT.Macro('rootlogon.C')

#draw plots on a log scale?
draw_log = False

#OptionsParser
parser = OptionParser()
parser.add_option('--n', metavar='F', type='string', action='store',dest='outname',default='NEW_PLOTS',help='')
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
linestyles = []
sample_weights = []
#Single top
filenames.append('T_s'); 	 fillcolors.append(kMagenta);	linestyles.append(1)
filenames.append('T_t'); 	 fillcolors.append(kMagenta);	linestyles.append(2)
filenames.append('T_tW'); 	 fillcolors.append(kMagenta);	linestyles.append(3)
filenames.append('Tbar_s');  fillcolors.append(kMagenta);	linestyles.append(4)
filenames.append('Tbar_t');  fillcolors.append(kMagenta);	linestyles.append(5)
filenames.append('Tbar_tW'); fillcolors.append(kMagenta);	linestyles.append(6)
#DYnJets
filenames.append('DY1Jets'); fillcolors.append(kAzure-2);	linestyles.append(1)
filenames.append('DY2Jets'); fillcolors.append(kAzure-2);	linestyles.append(2)
filenames.append('DY3Jets'); fillcolors.append(kAzure-2);	linestyles.append(3)
filenames.append('DY4Jets'); fillcolors.append(kAzure-2);	linestyles.append(4)
#WnJets samples
filenames.append('W1Jets'); fillcolors.append(kGreen-3);	linestyles.append(1)
filenames.append('W2Jets'); fillcolors.append(kGreen-3);	linestyles.append(2)
filenames.append('W3Jets'); fillcolors.append(kGreen-3);	linestyles.append(3)
filenames.append('W4Jets'); fillcolors.append(kGreen-3);	linestyles.append(4)
#POWHEG TT
#dileptonic 
filenames.append('Powheg_dilep_TT'); 				 fillcolors.append(kRed-7);	linestyles.append(1)
filenames.append('Powheg_dilep_TT_SC'); 			 fillcolors.append(kRed-7);	linestyles.append(2)
filenames.append('Powheg_dilep_TT_Mtt_700_to_1000'); fillcolors.append(kRed-7);	linestyles.append(3)
filenames.append('Powheg_dilep_TT_Mtt_1000_to_Inf'); fillcolors.append(kRed-7);	linestyles.append(4)
##hadronic
#filenames.append('Powheg_had_TT'); 				   fillcolors.append(kRed-7);	linestyles.append(1)
#filenames.append('Powheg_had_TT_SC'); 			   fillcolors.append(kRed-7);	linestyles.append(2)
#filenames.append('Powheg_had_TT_Mtt_700_to_1000'); fillcolors.append(kRed-7);	linestyles.append(3)
#filenames.append('Powheg_had_TT_Mtt_1000_to_Inf'); fillcolors.append(kRed-7);	linestyles.append(4)
#semileptonic qq
filenames.append('Powheg_qq_semilep_TT'); 				  fillcolors.append(kRed+1);	linestyles.append(1)
filenames.append('Powheg_qq_semilep_TT_SC'); 			  fillcolors.append(kRed+1);	linestyles.append(2)
filenames.append('Powheg_qq_semilep_TT_Mtt_700_to_1000'); fillcolors.append(kRed+1);	linestyles.append(3)
filenames.append('Powheg_qq_semilep_TT_Mtt_1000_to_Inf'); fillcolors.append(kRed+1);	linestyles.append(4)
#semileptonic gg
filenames.append('Powheg_gg_semilep_TT'); 				  fillcolors.append(kRed+1);	linestyles.append(1)
filenames.append('Powheg_gg_semilep_TT_SC'); 			  fillcolors.append(kRed+1);	linestyles.append(2)
filenames.append('Powheg_gg_semilep_TT_Mtt_700_to_1000'); fillcolors.append(kRed+1);	linestyles.append(3)
filenames.append('Powheg_gg_semilep_TT_Mtt_1000_to_Inf'); fillcolors.append(kRed+1);	linestyles.append(4)
#data
data_filename = 'Powheg_semilep_TT'
if leptype == 'mu' :
	data_filename = 'SingleMu_Run2012'
elif leptype == 'el' :
	data_filename = 'SingleEl_Run2012'

for i in range(len(filenames)) :
	filenames[i] = filenames[i] + '_all.root'
data_filename = data_filename+'_all.root'

#list of files
filelist = []
for filename in filenames :
	filelist.append(TFile(filename))
data_file = TFile(data_filename)
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
renormalize = []
includedata = []

#control plots
control_plot_names  = []
control_plot_lines  = []
control_plot_arrows = []
set_maxes = []
include_data = []
sum_same_colors = []

#Cut Strings
MARC_CUTS_MUONS = 'lepW_pt>50. && muon1_pt>10. && muon1_pt>ele1_pt && abs(muon1_eta)<2.4 && muon1_isLoose==1'
MARC_CUTS_MUONS+= ' && (muon1_relPt>25. || muon1_dR>0.5)'
MARC_CUTS_MUONS+= ' && lepb_M<50. && lept_M>140. && lept_M<250. && hadt_pt>400. && hadt_M>140. && hadt_M<250.'
MARC_CUTS_MUONS+= ' && hadt_tau32<0.55 && hadt_tau21>0.1 && cutflow==0'
MORE_SELECTIVE_MUONS = MARC_CUTS_MUONS+' && (muon2_isLoose!=1 || muon2_pt<40. || abs(muon2_eta)>2.4)'
MORE_SELECTIVE_MUONS+=' && (ele1_isLoose!=1 || ele1_pt<25. || abs(ele1_eta)>2.4)'
MORE_SELECTIVE_MUONS+=' && (ele2_isLoose!=1 || ele2_pt<25. || abs(ele2_eta)>2.4)'
MARC_CUTS_ELECTRONS = 'lepW_pt>50. && ele1_pt>25. && ele1_pt>muon1_pt && abs(ele1_eta)<2.4 && ele1_isLoose==1'
MARC_CUTS_ELECTRONS+= ' && (ele1_relPt>25. || ele1_dR>0.5)'
MARC_CUTS_ELECTRONS+= ' && lepb_M<50. && lept_M>140. && lept_M<250. && hadt_pt>400. && hadt_M>140. && hadt_M<250.'
MARC_CUTS_ELECTRONS+= ' && hadt_tau32<0.55 && hadt_tau21>0.1 && cutflow==0'
MORE_SELECTIVE_ELECTRONS = MARC_CUTS_ELECTRONS+' && (ele2_isLoose!=1 || ele2_pt<25. || abs(ele2_eta)>2.4)'
MORE_SELECTIVE_ELECTRONS+=' && (muon1_isLoose!=1 || muon1_pt<40. || abs(muon1_eta)>2.4)'
MORE_SELECTIVE_ELECTRONS+=' && (muon2_isLoose!=1 || muon2_pt<40. || abs(muon2_eta)>2.4)'
MARC_PLUS_TRIANGLE = MARC_CUTS_ELECTRONS+' && ele1_tri_el_val<ele1_tri_cut_val && ele1_tri_jet_val<ele1_tri_cut_val'
MORE_SELECTIVE_PLUS_TRIANGLE = MARC_PLUS_TRIANGLE+' && (ele2_isLoose!=1 || ele2_pt<25. || abs(ele2_eta)>2.4)'
MORE_SELECTIVE_PLUS_TRIANGLE+=' && (muon1_isLoose!=1 || muon1_pt<40. || abs(muon1_eta)>2.4)'
MORE_SELECTIVE_PLUS_TRIANGLE+=' && (muon2_isLoose!=1 || muon2_pt<40. || abs(muon2_eta)>2.4)'
#Weight strings
STD_WEIGHTS = 'weight*sf_pileup'

##############################################################################################################
############################						  PLOTS 					  ############################
##############################################################################################################

if leptype == 'mu' :

	#M
	histonames.append('M_muons')
	drawstrings.append('M')
	cutstrings.append(MARC_CUTS_MUONS)
	xbins = 25; xlow = 0.; xhigh = 2500.;
	title = 'M in Simulation and Data; M (GeV)'
	for i in range(len(filenames)) :
		histograms[i].append(TH1D(histonames[len(histonames)-1]+'_'+str(i),title+'_'+str(i),xbins,xlow,xhigh))
	datahistograms.append(TH1D(histonames[len(histonames)-1]+'_data',title,xbins,xlow,xhigh))
	histostacks.append(THStack(histonames[len(histonames)-1]+'_stack',title))
	optionstrings.append('')
	weightstrings.append(STD_WEIGHTS)
	includedata.append(True)
	renormalize.append(True)
	binningstrings.append('('+str(xbins)+','+str(xlow)+','+str(xhigh)+')')
	
	#M
	histonames.append('M_selective_muons')
	drawstrings.append('M')
	cutstrings.append(MORE_SELECTIVE_MUONS)
	xbins = 25; xlow = 0.; xhigh = 2500.;
	title = 'M in Simulation and Data; M (GeV)'
	for i in range(len(filenames)) :
		histograms[i].append(TH1D(histonames[len(histonames)-1]+'_'+str(i),title+'_'+str(i),xbins,xlow,xhigh))
	datahistograms.append(TH1D(histonames[len(histonames)-1]+'_data',title,xbins,xlow,xhigh))
	histostacks.append(THStack(histonames[len(histonames)-1]+'_stack',title))
	optionstrings.append('')
	weightstrings.append(STD_WEIGHTS)
	includedata.append(True)
	renormalize.append(True)
	binningstrings.append('('+str(xbins)+','+str(xlow)+','+str(xhigh)+')')

elif leptype == 'el' :

	#M
	histonames.append('M_electrons')
	drawstrings.append('M')
	cutstrings.append(MARC_CUTS_ELECTRONS)
	xbins = 25; xlow = 0.; xhigh = 2500.;
	title = 'M in Simulation and Data; M (GeV)'
	for i in range(len(filenames)) :
		histograms[i].append(TH1D(histonames[len(histonames)-1]+'_'+str(i),title+'_'+str(i),xbins,xlow,xhigh))
	datahistograms.append(TH1D(histonames[len(histonames)-1]+'_data',title,xbins,xlow,xhigh))
	histostacks.append(THStack(histonames[len(histonames)-1]+'_stack',title))
	optionstrings.append('')
	weightstrings.append(STD_WEIGHTS)
	includedata.append(True)
	renormalize.append(True)
	binningstrings.append('('+str(xbins)+','+str(xlow)+','+str(xhigh)+')')
	
	#M
	histonames.append('M_selective_electrons')
	drawstrings.append('M')
	cutstrings.append(MORE_SELECTIVE_ELECTRONS)
	xbins = 25; xlow = 0.; xhigh = 2500.;
	title = 'M in Simulation and Data; M (GeV)'
	for i in range(len(filenames)) :
		histograms[i].append(TH1D(histonames[len(histonames)-1]+'_'+str(i),title+'_'+str(i),xbins,xlow,xhigh))
	datahistograms.append(TH1D(histonames[len(histonames)-1]+'_data',title,xbins,xlow,xhigh))
	histostacks.append(THStack(histonames[len(histonames)-1]+'_stack',title))
	optionstrings.append('')
	weightstrings.append(STD_WEIGHTS)
	includedata.append(True)
	renormalize.append(True)
	binningstrings.append('('+str(xbins)+','+str(xlow)+','+str(xhigh)+')')
	
	#M
	histonames.append('M_electrons_plus_triangle')
	drawstrings.append('M')
	cutstrings.append(MARC_PLUS_TRIANGLE)
	xbins = 25; xlow = 0.; xhigh = 2500.;
	title = 'M in Simulation and Data; M (GeV)'
	for i in range(len(filenames)) :
		histograms[i].append(TH1D(histonames[len(histonames)-1]+'_'+str(i),title+'_'+str(i),xbins,xlow,xhigh))
	datahistograms.append(TH1D(histonames[len(histonames)-1]+'_data',title,xbins,xlow,xhigh))
	histostacks.append(THStack(histonames[len(histonames)-1]+'_stack',title))
	optionstrings.append('')
	weightstrings.append(STD_WEIGHTS)
	includedata.append(True)
	renormalize.append(True)
	binningstrings.append('('+str(xbins)+','+str(xlow)+','+str(xhigh)+')')
	
	#M
	histonames.append('M_selective_electrons_plus_triangle')
	drawstrings.append('M')
	cutstrings.append(MORE_SELECTIVE_PLUS_TRIANGLE)
	xbins = 25; xlow = 0.; xhigh = 2500.;
	title = 'M in Simulation and Data; M (GeV)'
	for i in range(len(filenames)) :
		histograms[i].append(TH1D(histonames[len(histonames)-1]+'_'+str(i),title+'_'+str(i),xbins,xlow,xhigh))
	datahistograms.append(TH1D(histonames[len(histonames)-1]+'_data',title,xbins,xlow,xhigh))
	histostacks.append(THStack(histonames[len(histonames)-1]+'_stack',title))
	optionstrings.append('')
	weightstrings.append(STD_WEIGHTS)
	includedata.append(True)
	renormalize.append(True)
	binningstrings.append('('+str(xbins)+','+str(xlow)+','+str(xhigh)+')')

##############################################################################################################
############################						  PLOTS 					  ############################
##############################################################################################################

#Draw data events into histograms
data_tree_unpruned = data_file.Get('tree')
prunestring = 'cutflow == 0'
if 'El' in data_filename :
	prunestring = MARC_CUTS_ELECTRONS
elif 'Mu' in data_filename :
	prunestring = MARC_CUTS_MUONS
data_tree = data_tree_unpruned.CopyTree(prunestring)
trees = []
print 'Drawing events from data file into histograms...'
for j in range(len(histonames)) :
	data_tree.Draw(drawstrings[j]+'>>tmp'+binningstrings[j]+'',cutstrings[j],' ')
	datahistograms[j]=(gDirectory.Get('tmp')).Clone(histonames[j]+'_data')
	datahistograms[j].SetDirectory(0)
print 'Done'
#Build renormalization values
leg = TLegend(0.62,0.67,0.9,0.9)
renorms = []
for j in range(len(histonames)) :
	renorms.append(0.0)
for i in range(len(filenames)) :
	unpruned = filelist[i].Get('tree')
	trees.append(unpruned.CopyTree(prunestring))
	print ( 'Drawing events from file '+filenames[i]+' into histograms to build renormalization values('
			+str(i+1)+' out of '+str(len(filenames))+')' )
	weight_array = array('d',[0.])
	trees[i].SetBranchAddress('weight',weight_array)
	trees[i].GetEntry(0)
	sample_weights.append(weight_array[0])
	for j in range(len(histonames)) :
		#Draw the plot I want from the tree
		print '	Drawing '+histonames[j]+' ('+str(j+1)+' out of '+str(len(histonames))+')'
		trees[i].Draw(drawstrings[j]+'>>tmp'+binningstrings[j]+'','('+weightstrings[j]+')*('+cutstrings[j]+')',optionstrings[j]+'')
		#Get the plot back
		print '		Cloning'
		histograms[i][j] = (gDirectory.Get('tmp')).Clone(histonames[j])
		histograms[i][j].SetDirectory(0)
		print '		Building renormalization value'
		renorms[j] = renorms[j]+histograms[i][j].Integral()
#Draw MC events into temporary histograms, add to stack
for i in range(len(filenames)) :
	for j in range(len(histonames)) :
		#Set Colors and Fills
		histograms[i][j].SetFillColor(fillcolors[i])
		histograms[i][j].SetLineColor(fillcolors[i])
		if renormalize[j] :
			histograms[i][j].SetMarkerStyle(21)
		else :
			histograms[i][j].SetMarkerStyle(25)
			histograms[i][j].SetLineWidth(3)
		#Renormalize to data if need be
		if renormalize[j] :
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
print 'plotting'
canvs = []
for i in range(len(histonames)) :
	tmp_canv_name = histonames[i]+'_canv'
	tmp_canv = TCanvas(tmp_canv_name,tmp_canv_name,1200,900)
	tmp_canv.cd()
	datamax = datahistograms[i].GetMaximum()+sqrt(datahistograms[i].GetMaximum())
	if renorms[i] :
		histostacks[i].SetMaximum(max(1.02*datamax,1.02*histostacks[i].GetMaximum()))
		histostacks[i].Draw()
	else :
		histostacks[i].SetMaximum(max(1.02*datamax,1.02*histostacks[i].GetMaximum("nostack")))
		histostacks[i].Draw("nostack")
	if includedata[i] :
		datahistograms[i].SetMarkerStyle(20)
		datahistograms[i].Draw("SAME PE1")
	leg.Draw()
	canvs.append(tmp_canv)

# Make a new root file and save the histogram stacks to it
name = options.outname
if draw_log :
	name += '_log'
else :
	name += '_linear'
f = TFile(name+'.root', 'Recreate' )
f.cd()
for canv in canvs :
	canv.Write()
f.Write()
f.Close()
for filep in filelist :
	filep.Close()
data_file.Close()
