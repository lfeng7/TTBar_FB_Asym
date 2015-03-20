#Takes in output files from the ttree maker and puts together the plots and stacks/lines them, then reoutputs them.
#imports
from ROOT import *
from optparse import OptionParser
from array import array
from math import *

#TDR Style
gROOT.Macro('rootlogon.C')

#OptionsParser
parser = OptionParser()
parser.add_option('--n', metavar='F', type='string', action='store',dest='outname',default='NEW_PLOTS_TYPE2_linear',help='')
(options, args) = parser.parse_args()

#lists of filenames, fill colors, and pointers
filenames = []
fillcolors = []
linestyles = []
sample_weights = []
#Single top
filenames.append('T_s_type2_all.root');		fillcolors.append(kMagenta);	linestyles.append(1)
filenames.append('T_t_type2_all.root');		fillcolors.append(kMagenta);	linestyles.append(2)
filenames.append('T_tW_type2_all.root');	fillcolors.append(kMagenta);	linestyles.append(3)
filenames.append('Tbar_s_type2_all.root');	fillcolors.append(kMagenta);	linestyles.append(4)
filenames.append('Tbar_t_type2_all.root');	fillcolors.append(kMagenta);	linestyles.append(5)
filenames.append('Tbar_tW_type2_all.root');	fillcolors.append(kMagenta);	linestyles.append(6)
#DYnJets
filenames.append('DY1Jets_type2_all.root');	fillcolors.append(kAzure-2);	linestyles.append(1)
filenames.append('DY2Jets_type2_all.root');	fillcolors.append(kAzure-2);	linestyles.append(2)
filenames.append('DY3Jets_type2_all.root');	fillcolors.append(kAzure-2);	linestyles.append(3)
filenames.append('DY4Jets_type2_all.root');	fillcolors.append(kAzure-2);	linestyles.append(4)
#WnJets samples
filenames.append('W1Jets_type2_all.root');	fillcolors.append(kGreen-3);	linestyles.append(1)
filenames.append('W2Jets_type2_all.root');	fillcolors.append(kGreen-3);	linestyles.append(2)
filenames.append('W3Jets_type2_all.root');	fillcolors.append(kGreen-3);	linestyles.append(3)
filenames.append('W4Jets_type2_all.root');	fillcolors.append(kGreen-3);	linestyles.append(4)
#POWHEG TT
#dileptonic 
filenames.append('Powheg_dilep_TT_type2_all.root');					fillcolors.append(kRed-7);	linestyles.append(1)
filenames.append('Powheg_dilep_TT_SC_type2_all.root');				fillcolors.append(kRed-7);	linestyles.append(2)
filenames.append('Powheg_dilep_TT_Mtt_700_to_1000_type2_all.root');	fillcolors.append(kRed-7);	linestyles.append(3)
filenames.append('Powheg_dilep_TT_Mtt_1000_to_Inf_type2_all.root');	fillcolors.append(kRed-7);	linestyles.append(4)
#hadronic
filenames.append('Powheg_had_TT_type2_all.root');					fillcolors.append(kRed-7);	linestyles.append(1)
filenames.append('Powheg_had_TT_SC_type2_all.root');				fillcolors.append(kRed-7);	linestyles.append(2)
filenames.append('Powheg_had_TT_Mtt_700_to_1000_type2_all.root');	fillcolors.append(kRed-7);	linestyles.append(3)
filenames.append('Powheg_had_TT_Mtt_1000_to_Inf_type2_all.root');	fillcolors.append(kRed-7);	linestyles.append(4)
#semileptonic qq
filenames.append('Powheg_qq_semilep_TT_type2_all.root');				 fillcolors.append(kRed+1);	linestyles.append(1)
filenames.append('Powheg_qq_semilep_TT_SC_type2_all.root');				 fillcolors.append(kRed+1);	linestyles.append(2)
filenames.append('Powheg_qq_semilep_TT_Mtt_700_to_1000_type2_all.root'); fillcolors.append(kRed+1);	linestyles.append(3)
filenames.append('Powheg_qq_semilep_TT_Mtt_1000_to_Inf_type2_all.root'); fillcolors.append(kRed+1);	linestyles.append(4)
#semileptonic gg
filenames.append('Powheg_gg_semilep_TT_type2_all.root');				 fillcolors.append(kRed+1);	linestyles.append(1)
filenames.append('Powheg_gg_semilep_TT_SC_type2_all.root');				 fillcolors.append(kRed+1);	linestyles.append(2)
filenames.append('Powheg_gg_semilep_TT_Mtt_700_to_1000_type2_all.root'); fillcolors.append(kRed+1);	linestyles.append(3)
filenames.append('Powheg_gg_semilep_TT_Mtt_1000_to_Inf_type2_all.root'); fillcolors.append(kRed+1);	linestyles.append(4)
#data
data_filename = 'SingleMu_Run2012_type2_all.root'

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
draw_log = []

##############################################################################################################
############################						  PLOTS 					  ############################
##############################################################################################################

#met pT
control_plot_names.append('met_pt')
control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
sum_same_colors.append(True); draw_log.append(False)
control_plot_lines[len(control_plot_lines)-1].append(TLine(10.,0.0,10.,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(10.,0.0,10.+40.,0.0,0.04,'|> SAME'))

#leading lepton pT
control_plot_names.append('lep1_pt')
control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
sum_same_colors.append(True); draw_log.append(False)
control_plot_lines[len(control_plot_lines)-1].append(TLine(45.,0.0,45.,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(45.,0.0,45.+40.,0.0,0.04,'|> SAME'))

#leading lepton eta
control_plot_names.append('lep1_eta')
control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
sum_same_colors.append(True); draw_log.append(False)
control_plot_lines[len(control_plot_lines)-1].append(TLine(2.1,0.0,2.1,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(2.1,0.0,2.1-0.4,0.0,0.04,'|> SAME'))
control_plot_lines[len(control_plot_lines)-1].append(TLine(-2.1,0.0,-2.1,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(-2.1,0.0,-2.1+0.4,0.0,0.04,'|> SAME'))

#second lepton pT
control_plot_names.append('lep2_pt')
control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
sum_same_colors.append(True); draw_log.append(False)
control_plot_lines[len(control_plot_lines)-1].append(TLine(45.,0.0,45.,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(45.,0.0,45.-40.,0.0,0.04,'|> SAME'))

#second lepton eta
control_plot_names.append('lep2_eta')
control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
sum_same_colors.append(True); draw_log.append(False)
control_plot_lines[len(control_plot_lines)-1].append(TLine(2.1,0.0,2.1,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(2.1,0.0,2.1+0.4,0.0,0.04,'|> SAME'))
control_plot_lines[len(control_plot_lines)-1].append(TLine(-2.1,0.0,-2.1,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(-2.1,0.0,-2.1-0.4,0.0,0.04,'|> SAME'))

#leading other lepton pT
control_plot_names.append('other_lep1_pt')
control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
sum_same_colors.append(True); draw_log.append(False)
control_plot_lines[len(control_plot_lines)-1].append(TLine(35.,0.0,35.,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(35.,0.0,35.-30.,0.0,0.04,'|> SAME'))

#leading other lepton eta
control_plot_names.append('other_lep1_eta')
control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
sum_same_colors.append(True); draw_log.append(False)
control_plot_lines[len(control_plot_lines)-1].append(TLine(2.5,0.0,2.5,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(2.5,0.0,2.5+0.4,0.0,0.04,'|> SAME'))
control_plot_lines[len(control_plot_lines)-1].append(TLine(-2.5,0.0,-2.5,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(-2.5,0.0,-2.5-0.4,0.0,0.04,'|> SAME'))

#leptonic bjet pT
control_plot_names.append('lep_bjet_pT')
control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
sum_same_colors.append(True); draw_log.append(False)
control_plot_lines[len(control_plot_lines)-1].append(TLine(25.,0.0,25.,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(25.,0.0,25.+30.,0.0,0.04,'|> SAME'))

#leptonic bjet dR
control_plot_names.append('lep_bjet_dPhi')
control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
sum_same_colors.append(True); draw_log.append(False)
control_plot_lines[len(control_plot_lines)-1].append(TLine(pi/2.,0.0,pi/2.,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(pi/2.,0.0,pi/2.-0.4,0.0,0.04,'|> SAME'))

#leptonic bjet combined mass
control_plot_names.append('lep_bjet_comb_mass')
control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
sum_same_colors.append(True); draw_log.append(False)
control_plot_lines[len(control_plot_lines)-1].append(TLine(140.,0.0,140.,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(140.,0.0,140.+40.,0.0,0.04,'|> SAME'))
control_plot_lines[len(control_plot_lines)-1].append(TLine(210.,0.0,210.,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(210.,0.0,210.-40.,0.0,0.04,'|> SAME'))

#leptonic bjet CSV
control_plot_names.append('lep_bjet_CSV')
control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
sum_same_colors.append(True); draw_log.append(False)
control_plot_lines[len(control_plot_lines)-1].append(TLine(0.244,0.0,0.244,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(0.244,0.0,0.244+0.15,0.0,0.04,'|> SAME'))

#type 1 control plots 
if 'type1' in filenames[0] :
	#top cand pT
	control_plot_names.append('t1_top_pT')
	control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
	sum_same_colors.append(True); draw_log.append(False)
	control_plot_lines[len(control_plot_lines)-1].append(TLine(400.,0.0,400.,0.0))
	control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(400.,0.0,400.+30.,0.0,0.04,'|> SAME'))

	#top cand mass
	control_plot_names.append('t1_top_mass')
	control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
	sum_same_colors.append(True); draw_log.append(False)
	control_plot_lines[len(control_plot_lines)-1].append(TLine(130.,0.0,130.,0.0))
	control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(130.,0.0,130.+30.,0.0,0.04,'|> SAME'))

	#top cand tau32
	control_plot_names.append('t1_top_tau32')
	control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
	sum_same_colors.append(True); draw_log.append(False)
	control_plot_lines[len(control_plot_lines)-1].append(TLine(0.7,0.0,0.7,0.0))
	control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(0.7,0.0,0.7-0.15,0.0,0.04,'|> SAME'))

	#top cand dR
	control_plot_names.append('t1_top_dPhi')
	control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
	sum_same_colors.append(True); draw_log.append(False)
	control_plot_lines[len(control_plot_lines)-1].append(TLine(pi/2.,0.0,pi/2.,0.0))
	control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(pi/2.,0.0,pi/2.+0.4,0.0,0.04,'|> SAME'))

	#top cand multiplicity
	control_plot_names.append('t1_top_mult')
	control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
	sum_same_colors.append(True); draw_log.append(False)
	control_plot_lines[len(control_plot_lines)-1].append(TLine(1.0,0.0,1.0,0.0))
	control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(1.0,0.0,1.0-0.2,0.0,0.04,'|> SAME'))

#type 2 control plots 
elif 'type2' in filenames[0] :
	#hadronic W pT
	control_plot_names.append('t2_top_W_pT')
	control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
	sum_same_colors.append(True); draw_log.append(False)
	control_plot_lines[len(control_plot_lines)-1].append(TLine(200.,0.0,200.,0.0))
	control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(200.,0.0,200.+30.,0.0,0.04,'|> SAME'))

	#hadronic W mass
	control_plot_names.append('t2_top_W_mass')
	control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
	sum_same_colors.append(True); draw_log.append(False)
	control_plot_lines[len(control_plot_lines)-1].append(TLine(65.,0.0,65.,0.0))
	control_plot_lines[len(control_plot_lines)-1].append(TLine(105.,0.0,105.,0.0))

	#hadronic W tau21
	control_plot_names.append('t2_top_W_tau21')
	control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
	sum_same_colors.append(True); draw_log.append(False)
	control_plot_lines[len(control_plot_lines)-1].append(TLine(0.75,0.0,0.75,0.0))
	control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(0.75,0.0,0.75-0.15,0.0,0.04,'|> SAME'))

	#hadronic W dR
	control_plot_names.append('t2_top_W_dPhi')
	control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
	sum_same_colors.append(True); draw_log.append(False)
	control_plot_lines[len(control_plot_lines)-1].append(TLine(pi/2.,0.0,pi/2.,0.0))
	control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(pi/2.,0.0,pi/2.+0.4,0.0,0.04,'|> SAME'))

	#hadronic W multiplicity
	control_plot_names.append('t2_top_W_mult')
	control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
	sum_same_colors.append(True); draw_log.append(False)
	control_plot_lines[len(control_plot_lines)-1].append(TLine(1.0,0.0,1.0,0.0))
	control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(1.0,0.0,1.0-0.2,0.0,0.04,'|> SAME'))

	#combined top mass
	control_plot_names.append('t2_top_comb_mass')
	control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
	sum_same_colors.append(True); draw_log.append(False)
	control_plot_lines[len(control_plot_lines)-1].append(TLine(140.,0.0,140.,0.0))
	control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(140.,0.0,140.+25.,0.0,0.04,'|> SAME'))
	control_plot_lines[len(control_plot_lines)-1].append(TLine(210.,0.0,210.,0.0))
	control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(210.,0.0,210.-25.,0.0,0.04,'|> SAME'))

	#hadronic b dR
	control_plot_names.append('t2_top_b_dPhi')
	control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
	sum_same_colors.append(True); draw_log.append(False)
	control_plot_lines[len(control_plot_lines)-1].append(TLine(pi/2.,0.0,pi/2.,0.0))
	control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(pi/2.,0.0,pi/2.+0.4,0.0,0.04,'|> SAME'))

	#hadronic b dR WRT to hadronic W
	control_plot_names.append('t2_top_b_W_dR')
	control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
	sum_same_colors.append(True); draw_log.append(False)
	control_plot_lines[len(control_plot_lines)-1].append(TLine(0.6,0.0,0.6,0.0))
	control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(0.6,0.0,0.6+0.4,0.0,0.04,'|> SAME'))

	#hadronic top best combined mass
	control_plot_names.append('t2_top_best_comb_mass')
	control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
	sum_same_colors.append(True); draw_log.append(False)


#Pileup Comparison
histonames.append('pileup_comp')
drawstrings.append('pileup')
cutstrings.append('cutflow==0')
xbins = 40; xlow = 0.0; xhigh = 40.0;
title = 'Pileup in Simulation and Data; npv'
for i in range(len(filenames)) :
	histograms[i].append(TH1D(histonames[len(histonames)-1]+'_'+str(i),title+'_'+str(i),xbins,xlow,xhigh))
datahistograms.append(TH1D(histonames[len(histonames)-1]+'_data',title,xbins,xlow,xhigh))
histostacks.append(THStack(histonames[len(histonames)-1]+'_stack',title))
optionstrings.append('')
weightstrings.append('weight*sf_top_pT*sf_btag_eff*sf_pileup*sf_lep_ID')
includedata.append(True)
renormalize.append(True)
binningstrings.append('('+str(xbins)+','+str(xlow)+','+str(xhigh)+')')

#costheta Comparison
histonames.append('costheta_comp')
drawstrings.append('cstar')
cutstrings.append('cutflow==0')
xbins = 20; xlow = -1.0; xhigh = 1.0;
title = 'c^{*} in Simulation and Data; c^{*}'
for i in range(len(filenames)) :
	histograms[i].append(TH1D(histonames[len(histonames)-1]+'_'+str(i),title+'_'+str(i),xbins,xlow,xhigh))
datahistograms.append(TH1D(histonames[len(histonames)-1]+'_data',title,xbins,xlow,xhigh))
histostacks.append(THStack(histonames[len(histonames)-1]+'_stack',title))
optionstrings.append('')
weightstrings.append('weight*sf_top_pT*sf_btag_eff*sf_pileup*sf_lep_ID')
includedata.append(True)
renormalize.append(True)
binningstrings.append('('+str(xbins)+','+str(xlow)+','+str(xhigh)+')')

#x_F Comparison
histonames.append('x_F_comp')
drawstrings.append('x_F')
cutstrings.append('cutflow==0')
xbins = 30; xlow = 0.0; xhigh = 0.6;
title = '|x_{F}| in Simulation and Data; |x_{F}|'
for i in range(len(filenames)) :
	histograms[i].append(TH1D(histonames[len(histonames)-1]+'_'+str(i),title+'_'+str(i),xbins,xlow,xhigh))
datahistograms.append(TH1D(histonames[len(histonames)-1]+'_data',title,xbins,xlow,xhigh))
histostacks.append(THStack(histonames[len(histonames)-1]+'_stack',title))
optionstrings.append('')
weightstrings.append('weight*sf_top_pT*sf_btag_eff*sf_pileup*sf_lep_ID')
includedata.append(True)
renormalize.append(True)
binningstrings.append('('+str(xbins)+','+str(xlow)+','+str(xhigh)+')')

#M Comparison
histonames.append('M_comp')
drawstrings.append('M')
cutstrings.append('cutflow==0')
xbins = 40; xlow = 350.; xhigh = 1750.;
title = 'M in Simulation and Data; M (GeV)'
for i in range(len(filenames)) :
	histograms[i].append(TH1D(histonames[len(histonames)-1]+'_'+str(i),title+'_'+str(i),xbins,xlow,xhigh))
datahistograms.append(TH1D(histonames[len(histonames)-1]+'_data',title,xbins,xlow,xhigh))
histostacks.append(THStack(histonames[len(histonames)-1]+'_stack',title))
optionstrings.append('')
weightstrings.append('weight*sf_top_pT*sf_btag_eff*sf_pileup*sf_lep_ID')
includedata.append(True)
renormalize.append(True)
binningstrings.append('('+str(xbins)+','+str(xlow)+','+str(xhigh)+')')


##############################################################################################################
############################						  PLOTS 					  ############################
##############################################################################################################

#Draw data events into histograms
data_tree = data_file.Get('tree')
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
	trees.append(filelist[i].Get('tree'))
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

#Do all the control plots stuff
control_plots = []
control_plot_max_ys = []
control_plot_canvs = []
#Get each plot from the files
for i in range(len(control_plot_names)) :
	control_plots.append([]); control_plot_max_ys.append(0.)
	#get this plot from all the MC files
	for j in range(len(filelist)) :
		newPlot = filelist[j].Get('control_plots_folder/'+control_plot_names[i]).Clone()
		newPlot.SetDirectory(0)
		newPlot.Scale(sample_weights[j])
		if sum_same_colors[i] :
			#add to this sample all of the other samples of this color
			for k in range(len(filelist)) :
				if fillcolors[k] == fillcolors[j] :
					addPlot = filelist[k].Get('control_plots_folder/'+control_plot_names[i]).Clone()
					addPlot.SetDirectory(0)
					addPlot.Scale(sample_weights[k])
					newPlot.Add(addPlot)
		#normalize
		newPlot.Scale(1.0/newPlot.Integral())
		#reset maximum if necessary
		if 1.02*newPlot.GetMaximum() > control_plot_max_ys[i] :
			control_plot_max_ys[i] = 1.02*newPlot.GetMaximum()
		#set line width and color accordingly
		newPlot.SetLineWidth(3); newPlot.SetLineColor(fillcolors[j]); newPlot.SetLineStyle(linestyles[j])
		control_plots[i].append(newPlot)
	#get this plot from the data file
	new_data_plot = data_file.Get('control_plots_folder/'+control_plot_names[i]).Clone()
	new_data_plot.SetDirectory(0)
	if new_data_plot.Integral() != 0. :
		new_data_plot.Scale(1.0/new_data_plot.Integral())
	if 1.02*new_data_plot.GetMaximum() > control_plot_max_ys[i] :
		control_plot_max_ys[i] = 1.02*new_data_plot.GetMaximum()
	#set line width and color accordingly
	new_data_plot.SetLineWidth(3); new_data_plot.SetLineColor(kBlack); new_data_plot.SetLineStyle(1)
	control_plots[i].append(new_data_plot)
	#If the line and arrow ys have to be reset, reset them.
	if set_maxes[i] :
		for line in control_plot_lines[i] :
			line.SetY2((0.98*control_plot_max_ys[i]))
		for arrow in control_plot_arrows[i] :
			arrow.SetY1((0.98*control_plot_max_ys[i])/2.)
			arrow.SetY2((0.98*control_plot_max_ys[i])/2.)
	#Set line and arrow widths and colors
	for line in control_plot_lines[i] :
		line.SetLineWidth(6); line.SetLineColor(kRed)
	for arrow in control_plot_arrows[i] :
		arrow.SetLineWidth(6); arrow.SetLineColor(kRed); arrow.SetAngle(40); arrow.SetFillColor(kRed)
	#make a canvas for this plot and draw all the plots on it
	tmp_canv_name = control_plot_names[i]+'_canv'
	tmp_canv = TCanvas(tmp_canv_name,tmp_canv_name,1200,900)
	if draw_log[i] :
		tmp_canv.SetLogy()
	tmp_canv.cd()
	drawn_colors = []
	control_plots[i][0].SetMaximum(control_plot_max_ys[i])
	control_plots[i][0].Draw()
	drawn_colors.append(fillcolors[0])
	for j in range(1,len(control_plots[i])-1) :
		if sum_same_colors[i] and fillcolors[j] in drawn_colors :
			continue
		control_plots[i][j].Draw("SAME")
		drawn_colors.append(fillcolors[j])
	if include_data[i] :
		control_plots[i][len(control_plots[i])-1].Draw("SAME")
	for line in control_plot_lines[i] :
		line.Draw("SAME")
	for arrow in control_plot_arrows[i] :
		arrow.Draw()
	leg.Draw()
	control_plot_canvs.append(tmp_canv)

# Make a new root file and save the histogram stacks to it
f = TFile(options.outname+'.root', 'Recreate' )
f.cd()
for canv in canvs :
	canv.Write()
for canv in control_plot_canvs :
	canv.Write()
f.Write()
f.Close()
for filep in filelist :
	filep.Close()
data_file.Close()
