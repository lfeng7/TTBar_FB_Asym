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
	resid_title = ';'
	for i in range(1,len(title.split(';'))) :
		resid_title+=title.split(';')[i]
		if i<len(title.split(';'))-1 :
			resid_title+=';'
	resids.append(TH1D(histonames[len(histonames)-1]+'_resid',resid_title,nbins,low,high))
	resid_lines.append(TLine(low,1.0,high,1.0))
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
filenames.append('T_s'); 	 fillcolors.append(kYellow)#kMagenta)
filenames.append('T_t'); 	 fillcolors.append(kYellow)#kMagenta)
filenames.append('T_tW'); 	 fillcolors.append(kYellow)#kMagenta)
filenames.append('Tbar_s');  fillcolors.append(kYellow)#kMagenta)
filenames.append('Tbar_t');  fillcolors.append(kYellow)#kMagenta)
filenames.append('Tbar_tW'); fillcolors.append(kYellow)#kMagenta)
#DYnJets
filenames.append('DY1Jets'); fillcolors.append(kYellow)#kAzure-2)
filenames.append('DY2Jets'); fillcolors.append(kYellow)#kAzure-2)
filenames.append('DY3Jets'); fillcolors.append(kYellow)#kAzure-2)
filenames.append('DY4Jets'); fillcolors.append(kYellow)#kAzure-2)
#WnJets samples
#filenames.append('W1Jets'); fillcolors.append(kGreen-3)
#filenames.append('W2Jets'); fillcolors.append(kGreen-3)
#filenames.append('W3Jets'); fillcolors.append(kGreen-3)
#filenames.append('W4Jets'); fillcolors.append(kGreen-3)
#POWHEG TT
#dileptonic 
filenames.append('Powheg_dilep_TT'); 				 fillcolors.append(kYellow)#kRed-7)
#filenames.append('Powheg_dilep_TT_SC'); 			 fillcolors.append(kRed-7)
#filenames.append('Powheg_dilep_TT_Mtt_700_to_1000'); fillcolors.append(kRed-7)
#filenames.append('Powheg_dilep_TT_Mtt_1000_to_Inf'); fillcolors.append(kRed-7)
#hadronic
filenames.append('Powheg_had_TT'); 				   fillcolors.append(kYellow)#kRed-5)
#filenames.append('Powheg_had_TT_SC'); 			   fillcolors.append(kRed-7)
#filenames.append('Powheg_had_TT_Mtt_700_to_1000'); fillcolors.append(kRed-7)
#filenames.append('Powheg_had_TT_Mtt_1000_to_Inf'); fillcolors.append(kRed-7)
#semileptonic gg
filenames.append('Powheg_gg_semilep_TT'); 				  fillcolors.append(kBlue)#kRed+1)
#filenames.append('Powheg_gg_semilep_TT_SC'); 			  fillcolors.append(kRed+1)
#filenames.append('Powheg_gg_semilep_TT_Mtt_700_to_1000'); fillcolors.append(kRed+1)
#filenames.append('Powheg_gg_semilep_TT_Mtt_1000_to_Inf'); fillcolors.append(kRed+1)
#semileptonic qq
filenames.append('Powheg_qq_semilep_TT'); 				  fillcolors.append(kRed+1)
#filenames.append('Powheg_qq_semilep_TT_SC'); 			  fillcolors.append(kRed+1)
#filenames.append('Powheg_qq_semilep_TT_Mtt_700_to_1000'); fillcolors.append(kRed+1)
#filenames.append('Powheg_qq_semilep_TT_Mtt_1000_to_Inf'); fillcolors.append(kRed+1)
#data
#data_filenames.append('Powheg_semilep_TT')
if leptype == 'mu' :
	data_filenames.append('SingleMu_Run2012A')
	data_filenames.append('SingleMu_Run2012B')
	data_filenames.append('SingleMu_Run2012C')
	data_filenames.append('SingleMu_Run2012D')
#elif leptype == 'el' :
	data_filenames.append('SingleEl_Run2012A')
	data_filenames.append('SingleEl_Run2012B')
	data_filenames.append('SingleEl_Run2012C')
	data_filenames.append('SingleEl_Run2012D')

#lists of files
filelists = []
data_filelists = []
for filename in filenames :
	path = '../total_ttree_files/'+filename+'_skim_all.root'
	filelists.append(glob.glob(path))
for filename in data_filenames :
	path = '../total_ttree_files/'+filename+'_skim_all.root'
	data_filelists.append(glob.glob(path))
#lists of histogram names, temporary histograms for MC and data, histogram stacks, strings of what to draw, cutstrings, 
#options strings, weight strings, and renormalization and include data options
histonames = []
histostacks = []
resids = []
resid_lines = []
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
skim = 'hadt_M>100. && lepb_pt>0.'
mu_trigger = 'mu_trigger==1'
el_trigger = 'el_trigger==1'
lep_trigger_plots = skim+' && (('+mu_trigger+') || ('+el_trigger+'))' #just for presentation plots
mu_kinematics = 'muon1_pt>40. && abs(muon1_eta)<2.4'
el_kinematics = 'ele1_pt>40. && abs(ele1_eta)<2.4'
lep_kinematics_plots = skim+' && (('+mu_trigger+' && '+mu_kinematics+') || ('+el_trigger+' && '+el_kinematics+'))' #just for presentation plots
mu_ID_iso = 'muon1_isLoose==1 && (muon1_relPt>25. || muon1_dR>0.5)'
el_ID_iso = 'ele1_isLoose==1  && (ele1_relPt>25.  || ele1_dR>0.5)'
lep_ID_iso_plots = skim+' && (('+mu_trigger+' && '+mu_kinematics+' && '+mu_ID_iso+') || ('+el_trigger+' && '+el_kinematics+' && '+el_ID_iso+'))'
lept_mass = 'lept_M>140. && lept_M<250.'
lept_mass_plots = lep_ID_iso_plots+' && '+lept_mass
hadt_pt = 'hadt_pt>300.'
hadt_pt_plots = lept_mass_plots+' && '+hadt_pt
hadt_mass = 'hadt_M>140. && hadt_M<250.'
hadt_mass_plots = hadt_pt_plots+' && '+hadt_mass
hadt_structure = 'hadt_tau21>0.1 && hadt_tau32<0.55'
hadt_structure_plots = hadt_mass_plots+' && '+hadt_structure
#electron trigger selection cuts
ELECTRON_TRIGGER_CUTS = 'mu_trigger==1 && muon1_pt>40. && abs(muon1_eta)<2.4 && muon1_isLoose==1 && (muon1_relPt>25. || muon1_dR>0.5)'
ELECTRON_TRIGGER_CUTS += ' && ele1_pt>40. && abs(ele1_eta)<2.4 && ele1_isLoose==1 && (ele1_relPt>25. || ele1_dR>0.5)'
ELECTRON_TRIGGER_CUTS += ' && (ele2_isLoose!=1 || ele2_pt<40. || abs(ele2_eta)>2.4 || (ele2_relPt<25. && ele2_dR<0.5))'
ELECTRON_TRIGGER_CUTS += ' && (muon2_isLoose!=1 || muon2_pt<40. || abs(muon2_eta)>2.4 || (muon2_relPt<25. && muon2_dR<0.5))'
ELECTRON_TRIGGER_CUTS += ' && muon1_Q+ele1_Q==0 && lept_M>140. && lept_M<250.'
ELECTRON_TRIGGER_CUTS += ' && lepb_pt>50. && hadt_pt>50. && max(lepb_pt,hadt_pt)>150.'
#electron ID selection cuts
#electron 1 and/or 2 selection
ele1_tag_selec = '(ele1_pt>40. && abs(ele1_eta)<2.4 && ele1_isLoose==1 && (ele1_relPt>25. || ele1_dR>0.5))'
ele1_probe_selec = '(ele1_pt>40. && abs(ele1_eta)<2.4 && (ele1_relPt>25. || ele1_dR>0.5))'
ele2_tag_selec = '(ele2_pt>40. && abs(ele2_eta)<2.4 && ele2_isLoose==1 && (ele2_relPt>25. || ele2_dR>0.5))'
ele2_probe_selec = '(ele2_pt>40. && abs(ele2_eta)<2.4 && (ele2_relPt>25. || ele2_dR>0.5))'
ELECTRON_ID_CUTS = '(('+ele1_tag_selec+' && '+ele2_probe_selec+') || ('+ele2_tag_selec+' && '+ele1_probe_selec+'))'
#muon 1 and 2 rejection
ELECTRON_ID_CUTS += ' && (muon1_pt<40. || abs(muon1_eta)>2.4 || muon1_isLoose!=1 || (muon1_relPt<25. && muon1_dR<0.5))'
ELECTRON_ID_CUTS += ' && (muon2_pt<40. || abs(muon2_eta)>2.4 || muon2_isLoose!=1 || (muon2_relPt<25. && muon2_dR<0.5))'
#require two electrons to have opposite charges
ELECTRON_ID_CUTS += ' && (ele1_Q+ele2_Q==0)'
#make the dielectron fourvector and cut on its mass
eem = 'sqrt((sqrt(ele1_pt*ele1_pt*cosh(ele1_eta)*cosh(ele1_eta)+ele1_M*ele1_M)+sqrt(ele2_pt*ele2_pt*cosh(ele2_eta)*cosh(ele2_eta)+ele2_M*ele2_M))'
eem+='*(sqrt(ele1_pt*ele1_pt*cosh(ele1_eta)*cosh(ele1_eta)+ele1_M*ele1_M)+sqrt(ele2_pt*ele2_pt*cosh(ele2_eta)*cosh(ele2_eta)+ele2_M*ele2_M))'
eem+='-(ele1_pt*cos(ele1_phi)+ele2_pt*cos(ele2_phi))*(ele1_pt*cos(ele1_phi)+ele2_pt*cos(ele2_phi))'
eem+='-(ele1_pt*sin(ele1_phi)+ele2_pt*sin(ele2_phi))*(ele1_pt*sin(ele1_phi)+ele2_pt*sin(ele2_phi))'
eem+='-(ele1_pt*sinh(ele1_eta)+ele2_pt*sinh(ele2_eta))*(ele1_pt*sinh(ele1_eta)+ele2_pt*sinh(ele2_eta)))'
ELECTRON_ID_CUTS+=' && '+eem+'>12.'
ELECTRON_ID_CUTS+=' && ('+eem+'<76. || '+eem+'>106.)'
#other leptonic side cuts
ELECTRON_ID_CUTS+=' && lepW_pt[0]>50.'
#other jet requirements
ELECTRON_ID_CUTS+=' && lepb_pt>50. && hadt_pt>50. && max(lepb_pt,hadt_pt)>150.'
ELECTRON_ID_CUTS+=' && (lepb_csv>0.244 || hadt_csv>0.244)'
#Weight strings
STD_WEIGHTS = '19748.*weight*sf_pileup*sf_top_pT*sf_lep_ID*sf_trig_eff*1.2'

##############################################################################################################
############################						  PLOTS 					  ############################
##############################################################################################################

full_selection_cuts = hadt_structure_plots
#if leptype == 'mu' :
#	full_selection_cuts = skim+' && '+mu_trigger+' && '+mu_kinematics+' && '+mu_ID_iso+' && '+lept_mass+' && '+hadt_pt+' && '+hadt_mass+' && '+hadt_structure
#elif leptype == 'el' :
#	full_selection_cuts = skim+' && '+el_trigger+' && '+el_kinematics+' && '+el_ID_iso+' && '+lept_mass+' && '+hadt_pt+' && '+hadt_mass+' && '+hadt_structure

##M
#newPlot('M','M',full_selection_cuts,25,0.,2500.,'M in Simulation and Data; M (GeV)','',STD_WEIGHTS,True)
##c*
#newPlot('cstar','cstar',full_selection_cuts,20,-1.0,1.0,'c* in Simulation and Data; c*','',STD_WEIGHTS,True)
##xF
#newPlot('x_F','abs(x_F)',full_selection_cuts,30,0.,0.6,'x_{F} in Simulation and Data; x_{F}','',STD_WEIGHTS,True)
#M (scaled)
newPlot('M_scaled_skim','M_scaled',skim,25,0.,2500.,'scaled M in Simulation and Data, step 0; M (GeV); Events/100 GeV','',STD_WEIGHTS,True)
newPlot('M_scaled_trig','M_scaled',lep_trigger_plots,25,0.,2500.,'scaled M in Simulation and Data, step 1; M (GeV); Events/100 GeV','',STD_WEIGHTS,True)
newPlot('M_scaled_kine','M_scaled',lep_kinematics_plots,25,0.,2500.,'scaled M in Simulation and Data, step 2; M (GeV); Events/100 GeV','',STD_WEIGHTS,True)
newPlot('M_scaled_ID_iso','M_scaled',lep_ID_iso_plots,25,0.,2500.,'scaled M in Simulation and Data, step 3; M (GeV); Events/100 GeV','',STD_WEIGHTS,True)
#newPlot('M_scaled_lept_M','M_scaled',lept_mass_plots,25,0.,2500.,'scaled M in Simulation and Data, step 4; M (GeV); Events/100 GeV','',STD_WEIGHTS,True)
#newPlot('M_scaled_hadt_pt','M_scaled',hadt_pt_plots,25,0.,2500.,'scaled M in Simulation and Data, step 5; M (GeV); Events/100 GeV','',STD_WEIGHTS,True)
#newPlot('M_scaled_hadt_M','M_scaled',hadt_mass_plots,25,0.,2500.,'scaled M in Simulation and Data, step 6; M (GeV); Events/100 GeV','',STD_WEIGHTS,True)
#newPlot('M_scaled_hadt_s','M_scaled',hadt_structure_plots,25,0.,2500.,'scaled M in Simulation and Data, step 7; M (GeV); Events/100 GeV','',STD_WEIGHTS,True)
#c* (scaled)
#newPlot('cstar_scaled','cstar_scaled',full_selection_cuts,10,-1.0,1.0,'scaled c* in Simulation and Data; c*','',STD_WEIGHTS,True)
#xF (scaled)
#newPlot('x_F_scaled','abs(x_F_scaled)',full_selection_cuts,30,0.,0.6,'scaled x_{F} in Simulation and Data; x_{F}','',STD_WEIGHTS,True)

##leading muon pT
#newPlot('muon1_pt','muon1_pt',MARC_PRESELECTION,20,0.,100.,'leading muon p_{T}; p_{T} (GeV)','',STD_WEIGHTS,True)
###leading muon eta
##newPlot('muon1_eta','muon1_eta',MARC_PRESELECTION,35,-3.5,3.5,'leading muon #eta; #eta','',STD_WEIGHTS,True)
###second leading muon pT
##newPlot('muon2_pt','muon2_pt',MARC_PRESELECTION,20,0.,100.,'second leading muon p_{T}; p_{T} (GeV)','',STD_WEIGHTS,True)
###second leading muon eta
##newPlot('muon2_eta','muon2_eta',MARC_PRESELECTION,35,-3.5,3.5,'second leading muon #eta; #eta','',STD_WEIGHTS,True)
#leading electron pT
#newPlot('ele1_pt','ele1_pt',ELECTRON_ID_CUTS,15,0.,300.,'leading electron p_{T}; p_{T} (GeV)','',STD_WEIGHTS,True)
#leading electron eta
#newPlot('ele1_eta','ele1_eta',ELECTRON_ID_CUTS,35,-3.5,3.5,'leading electron #eta; #eta','',STD_WEIGHTS,True)
#second leading electron pT
#newPlot('ele2_pt','ele2_pt',ELECTRON_ID_CUTS,15,0.,300.,'second leading electron p_{T}; p_{T} (GeV)','',STD_WEIGHTS,True)
#second leading electron eta
#newPlot('ele2_eta','ele2_eta',ELECTRON_ID_CUTS,35,-3.5,3.5,'second leading electron #eta; #eta','',STD_WEIGHTS,True)
##leptonic W pT
#newPlot('lepW_pt','lepW_pt',MARC_PRESELECTION,20,0.,200.,'leptonic W candidate p_{T}; p_{T} (GeV)','',STD_WEIGHTS,True)
##leptonic b mass
#newPlot('lepb_M','lepb_M',MARC_PRESELECTION,20,0.,100.,'leptonic b candidate mass; M (GeV)','',STD_WEIGHTS,True)
##leptonic t mass
#newPlot('lept_M','lept_M',MARC_PRESELECTION,35,0.,350.,'leptonic top candidate mass; M (GeV)','',STD_WEIGHTS,True)
###hadronic t pT
##newPlot('hadt_pt','hadt_pt',ele_hadronic_pretag,50,250.,750.,'hadronic top candidate p_{T}; p_{T} (GeV)','',STD_WEIGHTS,True)
###hadronic t mass
#newPlot('hadt_M','hadt_M',ele_hadronic_pretag,35,0.,350.,'hadronic top candidate mass; M (GeV)','',STD_WEIGHTS,True)
#hadronic t pt
#newPlot('hadt_pt','hadt_pt',full_selection_cuts,25,0.,1000.,'hadronic top candidate p_{T}; p_{T} (GeV)','',STD_WEIGHTS,True)
###hadronic t tau32
#newPlot('hadt_tau32','hadt_tau32',ele_hadronic_pretag,20,0.,1.,'hadronic top candidate #tau_{32}; #tau_{32}','',STD_WEIGHTS,True)
###hadronic t tau21
##newPlot('hadt_tau21','hadt_tau21',MARC_PRESELECTION,20,0.,1.,'hadronic top candidate #tau_{21}; #tau_{21}','',STD_WEIGHTS,True)
#leading jet pT
#newPlot('jet1_pt','max(lepb_pt,hadt_pt)',ELECTRON_ID_CUTS,25,0.,1000.,'leading jet p_{T}; p_{T} (GeV)','',STD_WEIGHTS,True)
#sub-leading jet pT
#newPlot('jet2_pt','min(lepb_pt,hadt_pt)',ELECTRON_ID_CUTS,25,0.,500.,'sub-leading jet p_{T}; p_{T} (GeV)','',STD_WEIGHTS,True)

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
prunestring = ''
#prunestring = 'muon1_isLoose || ele1_isLoose'
#prunestring = MARC_PRESELECTION
#if leptype == 'el' :
#	prunestring = MARC_PRESELECTION+' && ele1_pt>muon1_pt'
#elif leptype == 'mu' :
#	prunestring = MARC_PRESELECTION+' && muon1_pt>ele1_pt'

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
		xtra = ''
		if prunestring!='' :
			xtra += ' && '
		xtra+='weight!=1.0'
		MC_chains[i].CopyTree(prunestring+xtra).Draw(drawstrings[j]+'>>tmp','('+weightstrings[j]+')*('+cutstrings[j]+')',optionstrings[j]+'')
		histograms[i][j].Add(tmp)
		#Set Colors and Fills
		histograms[i][j].SetFillColor(fillcolors[i])
		histograms[i][j].SetLineColor(fillcolors[i])
		histograms[i][j].SetMarkerStyle(21)
		histograms[i][j].SetDirectory(0)
#		renorms[j] = renorms[j]+histograms[i][j].Integral()

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
#	if i==0 :
#		leg.AddEntry(histograms[i][0],"Single Top","F")
#	elif i==-1:#6 :
#		leg.AddEntry(histograms[i][0],"Z/#gamma+Jets","F")
#	elif i==-1 :
#		leg.AddEntry(histograms[i][0],"W+Jets","F")
#	elif i==6:#10 :
#		leg.AddEntry(histograms[i][0],"Dileptonic t#bar{t}","F")
#	elif i==7:#11 :
#		leg.AddEntry(histograms[i][0],"Hadronic t#bar{t}","F")
#	elif i==len(filenames)-1 :
#		leg.AddEntry(histograms[i][0],"Semileptonic t#bar{t}","F")
	if i==0 :
		leg.AddEntry(histograms[i][0],"Background","F")
	elif i==len(filenames)-2 :
		leg.AddEntry(histograms[i][0],"gg/qg #rightarrow t#bar{t}","F")
	elif i==len(filenames)-1 :
		leg.AddEntry(histograms[i][0],"q#bar{q} #rightarrow t#bar{t}","F")
leg.AddEntry(datahistograms[0],"2012 Data","LPE")

#Make residual plots
for j in range(len(histonames)) :
	nhistobins = histograms[0][j].GetNbinsX()
	maxdev = 0.
	for i in range(nhistobins) :
		xvalue = histograms[0][j].GetBinCenter(i)
		xerr   = histograms[0][j].GetBinWidth(i)/2.
		devents = datahistograms[j].GetBinContent(i) 
		mcevents = histostacks[j].GetStack().Last().GetBinContent(i)
		yvalue = 1.0; yerr = 1.0
		if mcevents==0. or devents==0. :
			continue
		yvalue = devents/mcevents
		yerr = yvalue*sqrt(1./devents+1./mcevents)
		if abs(yvalue+yerr-1.0) > maxdev :
			maxdev = abs(yvalue+yerr-1.0)
		if abs(yvalue-yerr-1.0) > maxdev :
			maxdev = abs(yvalue-yerr-1.0)
		resids[j].SetBinContent(i,yvalue)
		resids[j].SetBinError(i,yerr)
	#set residue plot attributes
	resids[j].GetXaxis().SetLabelSize((0.05*0.72)/0.28); resids[j].GetXaxis().SetTitleOffset(0.8)
	resids[j].GetYaxis().SetLabelSize((0.05*0.72)/0.28); resids[j].GetYaxis().SetTitleOffset(0.4)
	resids[j].GetXaxis().SetTitleSize((0.72/0.28)*resids[j].GetXaxis().GetTitleSize())
	resids[j].GetYaxis().SetTitleSize((0.72/0.28)*resids[j].GetYaxis().GetTitleSize())
	resids[j].GetYaxis().SetTitle('Data/MC')
	print 'maxdev = '+str(maxdev)
	resids[j].GetYaxis().SetRangeUser(1.0-(1.05*maxdev),1.0+(1.05*maxdev))
	resids[j].GetYaxis().SetNdivisions(503)
	resids[j].SetStats(0)
	resid_lines[j].SetLineWidth(2); resid_lines[j].SetLineStyle(2)

#make a new file
name = options.outname
if leptype == 'mu' :
	name += '_muons'
elif leptype == 'el' :
	name += '_electrons'
f = TFile(name+'.root', 'Recreate' )

#plot on canvases
print 'plotting stacked plots'
for i in range(len(histonames)) :
	tmp_canv_name = histonames[i]+'_canv'
	tmp_canv = TCanvas(tmp_canv_name,tmp_canv_name,1200,900)
	tmp_canv.cd()
	histo_pad=TPad(histonames[i]+'_histo_pad',histonames[i]+'_histo_pad',0,0.25,1,1)
	resid_pad=TPad(histonames[i]+'_resid_pad',histonames[i]+'_resid_pad',0,0,1.,0.25)
	histo_pad.SetCanvas(tmp_canv); resid_pad.SetCanvas(tmp_canv)
	histo_pad.SetLeftMargin(0.16); histo_pad.SetRightMargin(0.05) 
	histo_pad.SetTopMargin(0.11);	 histo_pad.SetBottomMargin(0.02)
	histo_pad.SetBorderMode(0)
	resid_pad.SetLeftMargin(0.16); resid_pad.SetRightMargin(0.05)
	resid_pad.SetTopMargin(0.0);   resid_pad.SetBottomMargin(0.3)
	resid_pad.SetBorderMode(0)
	resid_pad.Draw(); histo_pad.Draw()
	histo_pad.cd(); 
	datamax = datahistograms[i].GetMaximum()+sqrt(datahistograms[i].GetMaximum())
	histostacks[i].SetMaximum(max(1.02*datamax,1.02*histostacks[i].GetMaximum()))
	histostacks[i].Draw()
	histostacks[i].GetXaxis().SetLabelOffset(999)
	if includedata[i] :
		datahistograms[i].SetMarkerStyle(20)
		datahistograms[i].DrawCopy("SAME PE1")
	leg.Draw()
	resid_pad.cd(); 
	resids[i].Draw('PE1X0'); resid_lines[i].Draw('SAME')
	tmp_canv.Update()
	f.cd()
	tmp_canv.Write()

##Make line plots
#unique_histos = []
#histo_colors = []
#unique_histo_maxes = []
##look through all the MC filetypes
#print 'making line plots'
#for i in range(len(filenames)) :
#	#if we haven't yet seen a fill of this color add a new list of unique histograms and copy over
#	if fillcolors[i] not in histo_colors :
#		histo_colors.append(fillcolors[i])
#		unique_histos.append([])
#		for j in range(len(histonames)) :
#			unique_histos[len(unique_histos)-1].append(histograms[i][j].Clone())
#			unique_histos[len(unique_histos)-1][j].SetDirectory(0)
#	#otherwise find the point in the list of unique histograms with histograms of this color and add these
#	else : 
#		for j in range(len(histonames)) :
#			unique_histos[histo_colors.index(fillcolors[i])][j].Add(histograms[i][j].Clone())
#
##renormalize all the unique histograms to one and set the marker and line styles
#for i in range(len(unique_histos)) :
#	for j in range(len(unique_histos[i])) :
#		if unique_histos[i][j].Integral() != 0. :
#			unique_histos[i][j].Scale(1.0/unique_histos[i][j].Integral())
#		unique_histos[i][j].SetFillStyle(0)
#		unique_histos[i][j].SetLineWidth(3)
#		unique_histos[i][j].SetLineStyle(1)
##renormalize the data histograms, too
#for datahisto in datahistograms :
#	if datahisto.Integral() != 0. :
#		datahisto.Scale(1.0/datahisto.Integral())
#	datahisto.SetFillStyle(0)
#	datahisto.SetLineWidth(3)
#	datahisto.SetFillColor(kBlack)
#	datahisto.SetLineColor(kBlack)
#	datahisto.SetLineStyle(1)
##find the max value on each unique histogram
#for i in range(len(histonames)) :
#	unique_histo_maxes.append(0.)
#	for j in range(len(unique_histos)) :
#		if 1.02*unique_histos[j][i].GetMaximum() > unique_histo_maxes[i] :
#			unique_histo_maxes[i] = 1.02*unique_histos[j][i].GetMaximum()
##rebuild the legend
#leg = TLegend(0.62,0.67,0.9,0.9)
#for i in range(len(filenames)) :
#	if i==0 :
#		leg.AddEntry(histograms[i][0],"Single Top","F")
#	elif i==6 :
#		leg.AddEntry(histograms[i][0],"Z/#gamma+Jets","F")
#	elif i==10 :
#		leg.AddEntry(histograms[i][0],"W+Jets","F")
#	elif i==14 :
#		leg.AddEntry(histograms[i][0],"Dileptonic/Hadronic t#bar{t}","F")
#	elif i==len(filenames)-1 :
#		leg.AddEntry(histograms[i][0],"Semileptonic t#bar{t}","F")
#leg.AddEntry(datahistograms[0],"Data","F")
#
##plot on canvases
#other_canvs = []
#for i in range(len(histonames)) :
#	tmp_canv_name = histonames[i]+'_line_canv'
#	tmp_canv = TCanvas(tmp_canv_name,tmp_canv_name,1200,900)
#	tmp_canv.cd()
#	other_canvs.append(tmp_canv)
#	datamax = datahistograms[i].GetMaximum()
#	unique_histos[0][i].SetMaximum(max(1.02*datamax,unique_histo_maxes[i]))
#	unique_histos[0][i].Draw('')
#	for j in range(1,len(unique_histos)) :
#		unique_histos[j][i].Draw('SAME')
#	if includedata[i] :
#		datahistograms[i].Draw('SAME')
#	leg.Draw()

#save in file
print 'saving plots'
#save stuff
f.cd()
#for canv in other_canvs :
#	canv.Write()
f.Write()
f.Close()
