#Takes in output files from the ttree maker and puts together the plots and stacks/lines them, then reoutputs them.
#imports
from ROOT import *
from optparse import OptionParser

#TDR Style
gROOT.Macro('rootlogon.C')

#OptionsParser
parser = OptionParser()
parser.add_option('--n', metavar='F', type='string', action='store',dest='outname',default='NEW_PLOTS',help='')
(options, args) = parser.parse_args()

#lists of filenames, fill colors, and pointers
filenames = []
fillcolors = []
linestyles = []
#Single top
filenames.append('T_s_type1_all.root');		fillcolors.append(kMagenta);	linestyles.append(1)
filenames.append('T_t_type1_all.root');		fillcolors.append(kMagenta);	linestyles.append(2)
filenames.append('T_tW_type1_all.root');	fillcolors.append(kMagenta);	linestyles.append(3)
filenames.append('Tbar_s_type1_all.root');	fillcolors.append(kMagenta);	linestyles.append(4)
filenames.append('Tbar_t_type1_all.root');	fillcolors.append(kMagenta);	linestyles.append(5)
filenames.append('Tbar_tW_type1_all.root');	fillcolors.append(kMagenta);	linestyles.append(6)
#DYnJets
filenames.append('DY1Jets_type1_all.root');	fillcolors.append(kAzure-2);	linestyles.append(1)
filenames.append('DY2Jets_type1_all.root');	fillcolors.append(kAzure-2);	linestyles.append(2)
filenames.append('DY3Jets_type1_all.root');	fillcolors.append(kAzure-2);	linestyles.append(3)
filenames.append('DY4Jets_type1_all.root');	fillcolors.append(kAzure-2);	linestyles.append(4)
#WnJets samples
filenames.append('W1Jets_type1_all.root');	fillcolors.append(kGreen-3);	linestyles.append(1)
filenames.append('W2Jets_type1_all.root');	fillcolors.append(kGreen-3);	linestyles.append(2)
filenames.append('W3Jets_type1_all.root');	fillcolors.append(kGreen-3);	linestyles.append(3)
filenames.append('W4Jets_type1_all.root');	fillcolors.append(kGreen-3);	linestyles.append(4)
#POWHEG TT
#dileptonic 
filenames.append('Powheg_dilep_TT_type1_all.root');					fillcolors.append(kRed-7);	linestyles.append(1)
filenames.append('Powheg_dilep_TT_SC_type1_all.root');				fillcolors.append(kRed-7);	linestyles.append(2)
filenames.append('Powheg_dilep_TT_Mtt_700_to_1000_type1_all.root');	fillcolors.append(kRed-7);	linestyles.append(3)
filenames.append('Powheg_dilep_TT_Mtt_1000_to_Inf_type1_all.root');	fillcolors.append(kRed-7);	linestyles.append(4)
#hadronic
filenames.append('Powheg_had_TT_type1_all.root');					fillcolors.append(kRed-7);	linestyles.append(1)
filenames.append('Powheg_had_TT_SC_type1_all.root');				fillcolors.append(kRed-7);	linestyles.append(2)
filenames.append('Powheg_had_TT_Mtt_700_to_1000_type1_all.root');	fillcolors.append(kRed-7);	linestyles.append(3)
filenames.append('Powheg_had_TT_Mtt_1000_to_Inf_type1_all.root');	fillcolors.append(kRed-7);	linestyles.append(4)
#semileptonic qq
filenames.append('Powheg_qq_semilep_TT_type1_all.root');				 fillcolors.append(kRed+1);	linestyles.append(1)
filenames.append('Powheg_qq_semilep_TT_SC_type1_all.root');				 fillcolors.append(kRed+1);	linestyles.append(2)
filenames.append('Powheg_qq_semilep_TT_Mtt_700_to_1000_type1_all.root'); fillcolors.append(kRed+1);	linestyles.append(3)
filenames.append('Powheg_qq_semilep_TT_Mtt_1000_to_Inf_type1_all.root'); fillcolors.append(kRed+1);	linestyles.append(4)
#semileptonic gg
filenames.append('Powheg_gg_semilep_TT_type1_all.root');				 fillcolors.append(kRed+1);	linestyles.append(1)
filenames.append('Powheg_gg_semilep_TT_SC_type1_all.root');				 fillcolors.append(kRed+1);	linestyles.append(2)
filenames.append('Powheg_gg_semilep_TT_Mtt_700_to_1000_type1_all.root'); fillcolors.append(kRed+1);	linestyles.append(3)
filenames.append('Powheg_gg_semilep_TT_Mtt_1000_to_Inf_type1_all.root'); fillcolors.append(kRed+1);	linestyles.append(4)
#data
data_filename = 'SingleMu_Run2012_type1_all.root'

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

##############################################################################################################
############################						  PLOTS 					  ############################
##############################################################################################################

#met pT
control_plot_names.append('met_pt')
control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
control_plot_lines[len(control_plot_lines)-1].append(TLine(10.,0.0,10.,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(10.,0.0,10.+40.,0.0))

#leading lepton pT
control_plot_names.append('lep1_pt')
control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
control_plot_lines[len(control_plot_lines)-1].append(TLine(45.,0.0,45.,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(45.,0.0,45.+40.,0.0))

#leading lepton eta
control_plot_names.append('lep1_eta')
control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
control_plot_lines[len(control_plot_lines)-1].append(TLine(2.1,0.0,2.1,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(2.1,0.0,2.1-0.4,0.0,'<'))
control_plot_lines[len(control_plot_lines)-1].append(TLine(-2.1,0.0,-2.1,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(-2.1,0.0,-2.1+0.4,0.0,'>'))

#second lepton pT
control_plot_names.append('lep2_pt')
control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
control_plot_lines[len(control_plot_lines)-1].append(TLine(45.,0.0,45.,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(45.,0.0,45.-40.,0.0,'<'))

#second lepton eta
control_plot_names.append('lep2_eta')
control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
control_plot_lines[len(control_plot_lines)-1].append(TLine(2.1,0.0,2.1,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(2.1,0.0,2.1+0.4,0.0,'>'))
control_plot_lines[len(control_plot_lines)-1].append(TLine(-2.1,0.0,-2.1,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(-2.1,0.0,-2.1-0.4,0.0,'<'))

#leading other lepton pT
control_plot_names.append('other_lep1_pt')
control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
control_plot_lines[len(control_plot_lines)-1].append(TLine(35.,0.0,35.,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(35.,0.0,35.-30.,0.0,'<'))

#leading other lepton eta
control_plot_names.append('other_lep1_eta')
control_plot_lines.append([]); control_plot_arrows.append([]); set_maxes.append(True); include_data.append(True)
control_plot_lines[len(control_plot_lines)-1].append(TLine(2.5,0.0,2.5,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(2.5,0.0,2.5+0.4,0.0,'>'))
control_plot_lines[len(control_plot_lines)-1].append(TLine(-2.5,0.0,-2.5,0.0))
control_plot_arrows[len(control_plot_arrows)-1].append(TArrow(-2.5,0.0,-2.5-0.4,0.0,'<'))

#Pileup Comparison
histonames.append('pileup_comp')
drawstrings.append('pileup_events')
cutstrings.append('nbTags>1')
xbins = 40; xlow = 0.0; xhigh = 40.0;
histograms.append(TH1D(histonames[len(histonames)-1],'Pileup in Simulation and Data; npv',xbins,xlow,xhigh))
datahistograms.append(TH1D(histonames[len(histonames)-1]+'_data','Pileup in Simulation and Data; npv',xbins,xlow,xhigh))
histostacks.append(THStack(histonames[len(histonames)-1]+'_stack','Pileup in Simulation and Data; npv'))
optionstrings.append('')
weightstrings.append('weight_top_pT*weight_btag_eff*weight_pileup*weight_tracking*weight_lep_ID*weight_lep_iso*weight_trig_eff*weight_GJR_scale')
includedata.append(True)
renormalize.append(True)
binningstrings.append('('+str(xbins)+','+str(xlow)+','+str(xhigh)+')')


##############################################################################################################
############################						  PLOTS 					  ############################
##############################################################################################################

#Draw data events into histograms
data_tree = data_file.Get('output')
for j in range(len(histonames)) :
	data_tree.Draw(drawstrings[j]+'>>tmp'+binningstrings[j]+'','('+weightstrings[j]+')*('+cutstrings[j]+')',' ')
	datahistograms[j]=(gDirectory.Get('tmp')).Clone(histonames[j]+'_data')
#Build renormalization values
renorms = []
for j in range(len(histonames)) :
	renorms.append(1.0)
for i in range(len(filenames)) :
	tree = filelist[i].Get('output')
	for j in range(len(histonames)) :
		#Draw the plot I want from the tree
		tree.Draw(drawstrings[j]+'>>tmp'+binningstrings[j]+'','('+weightstrings[j]+')*('+cutstrings[j]+')',optionstrings[j]+'')
		#Get the plot back
		histograms[j] = (gDirectory.Get('tmp')).Clone(histonames[j])
		#Rescale by cross section/number of events, add to total integral of stack
		histograms[j].Scale((25523595./245.8)*cross_section[i]/eventsGenerated[i])
		renorms[j] = renorms[j]+histograms[j].Integral()
#Draw MC events into temporary histograms, add to stack
leg = TLegend(0.62,0.67,0.9,0.9)
for i in range(len(filenames)) :
	tree = filelist[i].Get('output')
	print 'File: '+filenames[i]
	for j in range(len(histonames)) :
		print '	output->Draw('+drawstrings[j]+'>>tmp'+binningstrings[j]+',('+weightstrings[j]+')*('+cutstrings[j]+'),'+optionstrings[j]+')'
		tree.Draw(drawstrings[j]+'>>tmp'+binningstrings[j]+'','('+weightstrings[j]+')*('+cutstrings[j]+')',optionstrings[j]+'')
		histograms[j] = (gDirectory.Get('tmp')).Clone(histonames[j])
		#Set Colors and Fills
		histograms[j].SetFillColor(fillcolors[i])
		histograms[j].SetLineColor(fillcolors[i])
		if renormalize[j] :
			histograms[j].SetMarkerStyle(21)
		else :
			histograms[j].SetMarkerStyle(25)
			histograms[j].SetLineWidth(3)
		#Rescale by cross section/number of events, add to total integral of stack
		histograms[j].Scale((25523595./245.8)*cross_section[i]/eventsGenerated[i])
		#Renormalize to data if need be
		if renormalize[j] :
			histograms[j].Scale(datahistograms[j].Integral()/renorms[j])
		#add to histogram stack
		histostacks[j].Add(histograms[j].Clone())
	#Add to legend if need be
	if i==0 :
		leg.AddEntry(histograms[0],"Single Top","F")
	elif i==6 :
		leg.AddEntry(histograms[0],"Z/#gamma+Jets","F")
	elif i==10 :
		leg.AddEntry(histograms[0],"W+Jets","F")
	elif i==14 :
		leg.AddEntry(histograms[0],"Dileptonic/Hadronic t#bar{t}","F")
	elif i==len(filenames)-1 :
		leg.AddEntry(histograms[0],"Semileptonic t#bar{t}","F")

leg.AddEntry(datahistograms[0],"Data","LPE")

#plot on canvases
print 'plotting'
canvs = []
for i in range(len(histonames)) :
	tmp_canv_name = histonames[i]+'_canv'
	tmp_canv = TCanvas(tmp_canv_name,tmp_canv_name,1200,900)
	tmp_canv.cd()
	if renorms[i] :
		histostacks[i].SetMaximum(max(1.02*datahistograms[i].GetMaximum(),1.02*histostacks[i].GetMaximum()))
		histostacks[i].Draw()
	else :
		histostacks[i].SetMaximum(max(1.02*datahistograms[i].GetMaximum(),1.02*histostacks[i].GetMaximum("nostack")))
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
		#normalize
		newPlot.Scale(1.0/newPlot.Integral())
		#reset maximum if necessary
		if 1.02*newPlot.GetMaximum() > control_plot_max_ys[i] :
			control_plot_max_ys[i] = 1.02*newPlot.GetMaximum()
		#set line width and color accordingly
		newPlot.SetLineWidth(2); newPlot.SetLineColor(fillcolors[j]); newPlot.SetLineStyle(linestyles[j])
		control_plots[i].append(newPlot)
	#get this plot from the data file
	new_data_plot = data_file.Get('control_plots_folder/'+control_plot_names[i]).Clone()
	new_data_plot.SetDirectory(0)
	new_data_plot.Scale(1.0/new_data_plot.Integral())
	if 1.02*new_data_plot.GetMaximum() > control_plot_max_ys[i] :
		control_plot_max_ys[i] = 1.02*new_data_plot.GetMaximum()
	#set line width and color accordingly
	new_data_plot.SetLineWidth(2); new_data_plot.SetLineColor(kBlack); new_data_plot.SetLineStyle(1)
	control_plots[i].append(new_data_plot)
	#If the line and arrow ys have to be reset, reset them.
	if set_maxes[i] :
		control_plot_lines[i].SetY1(control_plot_max_ys[i])
		control_plot_lines[i].SetY2(control_plot_max_ys[i])
		control_plot_arrows[i].SetY1(control_plot_max_ys[i]/2.)
		control_plot_arrows[i].SetY2(control_plot_max_ys[i]/2.)
	#make a canvas for this plot and draw all the plots on it
	tmp_canv_name = control_plot_names[i]+'_canv'
	tmp_canv = TCanvas(tmp_canv_name,tmp_canv_name,1200,900)
	tmp_canv.cd()
	control_plots[i][0].Draw("AL")
	for j in range(1,len(control_plots[i])-1) :
		control_plots[i][j].Draw("L SAME")
	if include_data[i] :
		control_plots[i][len(control_plots[i])-1].Draw("L SAME")
	for j in range(len(control_plot_lines[i])) :
		control_plot_lines[i][j].Draw("SAME")
	for j in range(len(control_plot_arrows[i])) :
		control_plot_arrows[i][j].Draw("SAME")
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
