from ROOT import *
from array import array
from math import *
import os

filenames = []
shortnames = []
weights = []
#Single top
filenames.append('T_s');		shortnames.append('Single Top')
filenames.append('T_t');		shortnames.append('Single Top')
filenames.append('T_tW');		shortnames.append('Single Top')
filenames.append('Tbar_s');		shortnames.append('Single Top')
filenames.append('Tbar_t');		shortnames.append('Single Top')
filenames.append('Tbar_tW');	shortnames.append('Single Top')
#DYnJets
filenames.append('DY1Jets');	shortnames.append('DYJets')
filenames.append('DY2Jets');	shortnames.append('DYJets')
filenames.append('DY3Jets');	shortnames.append('DYJets')
filenames.append('DY4Jets');	shortnames.append('DYJets')
#WnJets samples
filenames.append('W1Jets');		shortnames.append('WJets')
filenames.append('W2Jets');		shortnames.append('WJets')
filenames.append('W3Jets');		shortnames.append('WJets')
filenames.append('W4Jets');		shortnames.append('WJets')
#POWHEG TT
#dileptonic 
filenames.append('Powheg_dilep_TT');						shortnames.append('Dileptonic TTBar')
filenames.append('Powheg_dilep_TT_SC');						shortnames.append('Dileptonic TTBar')
filenames.append('Powheg_dilep_TT_Mtt_700_to_1000');		shortnames.append('Dileptonic TTBar')
filenames.append('Powheg_dilep_TT_Mtt_1000_to_Inf');		shortnames.append('Dileptonic TTBar')
#hadronic
filenames.append('Powheg_had_TT');							shortnames.append('Hadronic TTBar')
filenames.append('Powheg_had_TT_SC');						shortnames.append('Hadronic TTBar')
filenames.append('Powheg_had_TT_Mtt_700_to_1000');			shortnames.append('Hadronic TTBar')
filenames.append('Powheg_had_TT_Mtt_1000_to_Inf');			shortnames.append('Hadronic TTBar')
#semileptonic qq
filenames.append('Powheg_qq_semilep_TT');					shortnames.append('Semileptonic TTBar')
filenames.append('Powheg_qq_semilep_TT_SC');				shortnames.append('Semileptonic TTBar')
filenames.append('Powheg_qq_semilep_TT_Mtt_700_to_1000');	shortnames.append('Semileptonic TTBar')
filenames.append('Powheg_qq_semilep_TT_Mtt_1000_to_Inf');	shortnames.append('Semileptonic TTBar')
#semileptonic gg
filenames.append('Powheg_gg_semilep_TT');					shortnames.append('Semileptonic TTBar')
filenames.append('Powheg_gg_semilep_TT_SC');				shortnames.append('Semileptonic TTBar')
filenames.append('Powheg_gg_semilep_TT_Mtt_700_to_1000');	shortnames.append('Semileptonic TTBar')
filenames.append('Powheg_gg_semilep_TT_Mtt_1000_to_Inf');	shortnames.append('Semileptonic TTBar')
#data
data_filename = 'SingleMu_Run2012'

#Get the files
files = []
for i in range(len(filenames)) :
	filenames[i] += '_type1_all.root'
	files.append(TFile(filenames[i]))
data_file = TFile(data_filename+'_type1_all.root')

N_CUTFLOWS = 10+1
events_at_cutflow = []
#For each file, cut on each cutflow value, then rescale the integral to get the number of events
for i in range(len(files)) :
	print 'Getting event numbers for file '+filenames[i]+' ('+str(i+1)+' out of '+str(len(files))+')'
	events_at_cutflow.append([])
	tree = files[i].Get('tree')
	#Get the weight of each sample
	weight_array = array('d',[0.])
	tree.SetBranchAddress('weight',weight_array)
	tree.GetEntry(0)
	weights.append(weight_array[0])
	#Get the total number of events of this type	
	tree.Draw('cutflow>>tmp('+str(N_CUTFLOWS)+',0,'+str(N_CUTFLOWS)+')','(weight)')
	events_at_cutflow[i].append(1.0*gDirectory.Get('tmp').Integral())
	#Get the number of events at each cutflow except the last
	for cutflow in range(1,N_CUTFLOWS) :
		tree.Draw('cutflow>>tmp('+str(N_CUTFLOWS)+',0,'+str(N_CUTFLOWS)+')','(weight)*(cutflow==0 || cutflow>'+str(cutflow)+')')
		events_at_cutflow[i].append(1.0*gDirectory.Get('tmp').Integral())
	#Get the final number of selected events
	tree.Draw('cutflow>>tmp('+str(N_CUTFLOWS)+',0,'+str(N_CUTFLOWS)+')','(weight)*(cutflow==0)')
	events_at_cutflow[i].append(1.0*gDirectory.Get('tmp').Integral())
	files[i].Close()

print 'Getting event numbers from data file (get comfortable)'
data_events_at_cutflow = []
data_tree = data_file.Get('tree')
#Get the total number of events
data_tree.Draw('cutflow>>tmp('+str(N_CUTFLOWS)+',0,'+str(N_CUTFLOWS)+')')
data_events_at_cutflow.append(1.0*gDirectory.Get('tmp').Integral())
for cutflow in range(1,N_CUTFLOWS) :
	#Get the number of events at each cutflow
	data_tree.Draw('cutflow>>tmp('+str(N_CUTFLOWS)+',0,'+str(N_CUTFLOWS)+')','(cutflow==0 || cutflow>'+str(cutflow)+')')
	data_events_at_cutflow.append(1.0*gDirectory.Get('tmp').Integral())
#Get the total number of reconstructed events
data_tree.Draw('cutflow>>tmp('+str(N_CUTFLOWS)+',0,'+str(N_CUTFLOWS)+')','(cutflow==0)')
data_events_at_cutflow.append(1.0*gDirectory.Get('tmp').Integral())
data_file.Close()

#print out the number of events in data and the efficiencies for the MC samples
#first line is just table headings for each cutflow and each sample type
cutflow_filename = 'cutflow_table'
if 'type1' in filenames[0] :
	cutflow_filename+='_type1'
elif 'type2' in filenames[0] :
	cutflow_filename+='_type2'
cutflow_filename+='.txt'
first_line = 'Cutflow 	Data Events 		'
added_shortnames = []
for shortname in shortnames :
	if not shortname in added_shortnames :
		added_shortnames.append(shortname)
		first_line += shortname+' eff		'
print first_line
os.system('echo "'+first_line+'" > '+cutflow_filename+'')
#second line is just the total number of events in data
second_line = '	 	'+str(data_events_at_cutflow[0])+' 			'
for shortname in added_shortnames :
	second_line += '					'
print second_line
os.system('echo "'+second_line+'" >> '+cutflow_filename+'')
#Each line after that is the cutflow number, then the number of events in data, then the eff. with uncertainty for each sample type
for cutflow in range(1,N_CUTFLOWS+1) :
	#Begin with the cutflow number and the number of events in data
	next_line = ''
	if cutflow!=N_CUTFLOWS :
		next_line += str(cutflow)
	else :
		next_line += 'tot'
	next_line+=' 		%.2f'%(data_events_at_cutflow[cutflow])
	for i in range(len(added_shortnames)) :
		#calculate the efficiency for this sample type
		prev_number_of_events = 0.
		prev_number_of_events_variance = 0.
		number_of_events = 0.
		number_of_events_variance = 0.
		for j in range(len(shortnames)) :
			if shortnames[j] == added_shortnames[i] :
				prev_number_of_events += events_at_cutflow[j][cutflow-1]
				prev_number_of_events_variance += (weights[j])*events_at_cutflow[j][cutflow-1]
				number_of_events += events_at_cutflow[j][cutflow]
				number_of_events_variance += (weights[j])*events_at_cutflow[j][cutflow]
		eff = 1.0
		if prev_number_of_events != 0. :
			eff = number_of_events/prev_number_of_events
		if number_of_events == 0. :
			eff = 0.0
		eff_err = 1.0
		if number_of_events!=0. and prev_number_of_events!=0. :
			frac_var1 = number_of_events_variance/(number_of_events*number_of_events)
			frac_var2 = prev_number_of_events_variance/(prev_number_of_events*prev_number_of_events)
			eff_err = eff*sqrt(frac_var1+frac_var2)
		#add to the line
		next_line+='		%.6f (%.6f)'%(eff,eff_err)
	print next_line
	os.system('echo "'+next_line+'" >> '+cutflow_filename+'')