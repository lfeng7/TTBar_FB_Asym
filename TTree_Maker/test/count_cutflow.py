from ROOT import *
from array import array
from math import *
import os, glob

cutflow_filename = 'cutflow_table_muons.txt'

filenames = []
shortnames = []
weights = []
#POWHEG TT
#semileptonic qq
filenames.append('Powheg_qq_semilep_TT');					shortnames.append('Semileptonic TTBar')
#filenames.append('Powheg_qq_semilep_TT_SC');				shortnames.append('Semileptonic TTBar')
#filenames.append('Powheg_qq_semilep_TT_Mtt_700_to_1000');	shortnames.append('Semileptonic TTBar')
#filenames.append('Powheg_qq_semilep_TT_Mtt_1000_to_Inf');	shortnames.append('Semileptonic TTBar')
#semileptonic gg
filenames.append('Powheg_gg_semilep_TT');					shortnames.append('Semileptonic TTBar')
#filenames.append('Powheg_gg_semilep_TT_SC');				shortnames.append('Semileptonic TTBar')
#filenames.append('Powheg_gg_semilep_TT_Mtt_700_to_1000');	shortnames.append('Semileptonic TTBar')
#filenames.append('Powheg_gg_semilep_TT_Mtt_1000_to_Inf');	shortnames.append('Semileptonic TTBar')
#dileptonic 
filenames.append('Powheg_dilep_TT');						shortnames.append('Dileptonic TTBar')
#filenames.append('Powheg_dilep_TT_SC');						shortnames.append('Dileptonic TTBar')
#filenames.append('Powheg_dilep_TT_Mtt_700_to_1000');		shortnames.append('Dileptonic TTBar')
#filenames.append('Powheg_dilep_TT_Mtt_1000_to_Inf');		shortnames.append('Dileptonic TTBar')
#hadronic
filenames.append('Powheg_had_TT');							shortnames.append('Hadronic TTBar')
#filenames.append('Powheg_had_TT_SC');						shortnames.append('Hadronic TTBar')
#filenames.append('Powheg_had_TT_Mtt_700_to_1000');			shortnames.append('Hadronic TTBar')
#filenames.append('Powheg_had_TT_Mtt_1000_to_Inf');			shortnames.append('Hadronic TTBar')
#WnJets samples
filenames.append('W1Jets');		shortnames.append('WJets')
filenames.append('W2Jets');		shortnames.append('WJets')
filenames.append('W3Jets');		shortnames.append('WJets')
filenames.append('W4Jets');		shortnames.append('WJets')
#DYnJets
filenames.append('DY1Jets');	shortnames.append('DYJets')
filenames.append('DY2Jets');	shortnames.append('DYJets')
filenames.append('DY3Jets');	shortnames.append('DYJets')
filenames.append('DY4Jets');	shortnames.append('DYJets')
#Single top
filenames.append('T_s');		shortnames.append('Single Top')
filenames.append('T_t');		shortnames.append('Single Top')
filenames.append('T_tW');		shortnames.append('Single Top')
filenames.append('Tbar_s');		shortnames.append('Single Top')
filenames.append('Tbar_t');		shortnames.append('Single Top')
filenames.append('Tbar_tW');	shortnames.append('Single Top')

muon_data_filenames = []; ele_data_filenames = []
#data
muon_data_filenames.append('SingleMu_Run2012A')
muon_data_filenames.append('SingleMu_Run2012B')
muon_data_filenames.append('SingleMu_Run2012C')
muon_data_filenames.append('SingleMu_Run2012D')
ele_data_filenames.append('SingleEl_Run2012A')
ele_data_filenames.append('SingleEl_Run2012B')
ele_data_filenames.append('SingleEl_Run2012C')
ele_data_filenames.append('SingleEl_Run2012D')

#Chain up the files
MC_chains = []
muon_data_chain = TChain('tree')
ele_data_chain = TChain('tree')
shortnames_done = []
for i in range(len(filenames)) :
	index = 0
	if shortnames[i] not in shortnames_done :
		MC_chains.append(TChain('tree'))
		index = len(MC_chains)-1
		shortnames_done.append(shortnames[i])
	else :
		index = shortnames_done.index(shortnames[i])
	filenamelist = glob.glob('../'+filenames[i]+'/'+filenames[i]+'*_skim_tree.root')
	for filename in filenamelist :
		MC_chains[index].Add(filename)
for muon_data_filename in muon_data_filenames :
	filenamelist = glob.glob('../'+muon_data_filename+'/'+muon_data_filename+'*_skim_tree.root')
	for filename in filenamelist :
		muon_data_chain.Add(filename)
for ele_data_filename in ele_data_filenames :
	filenamelist = glob.glob('../'+ele_data_filename+'/'+ele_data_filename+'*_skim_tree.root')
	for filename in filenamelist :
		ele_data_chain.Add(filename)


#Cut details
cutnames = []; cutstrings = []; prior_cutstrings = []
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
cutnames.append('muon skim'); 			   cutstrings.append('muon1_pt>ele1_pt'); 										   prior_cutstrings.append('muon1_pt>ele1_pt')
cutnames.append('muon preselection'); 	   cutstrings.append(muon_preselection); 										   prior_cutstrings.append('muon1_pt>ele1_pt')
cutnames.append('muon kinematics'); 	   cutstrings.append(muon_preselection+' && '+muon_kinematics); 				   prior_cutstrings.append(muon_preselection)
cutnames.append('muon ID'); 			   cutstrings.append(muon_preselection+' && '+muon_ID); 						   prior_cutstrings.append(muon_preselection)
cutnames.append('muon 2D cut'); 		   cutstrings.append(muon_preselection+' && '+muon_2D); 						   prior_cutstrings.append(muon_preselection)
cutnames.append('muon leptonic top mass'); cutstrings.append(muon_preselection+' && '+lep_top_mass); 					   prior_cutstrings.append(muon_preselection)
cutnames.append('muon full leptonic'); 	   cutstrings.append(muon_full_leptonic); 										   prior_cutstrings.append(muon_preselection)
cutnames.append('muon hadronic pretag');   cutstrings.append(muon_hadronic_pretag); 									   prior_cutstrings.append(muon_preselection)
cutnames.append('muon signal mass'); 	   cutstrings.append(muon_hadronic_pretag+' && '+signal_mass); 					   prior_cutstrings.append(muon_hadronic_pretag)
cutnames.append('muon signal tau32'); 	   cutstrings.append(muon_hadronic_pretag+' && '+signal_tau32); 				   prior_cutstrings.append(muon_hadronic_pretag)
cutnames.append('muon full selection');    cutstrings.append(muon_hadronic_pretag+' && '+signal_mass+' && '+signal_tau32); prior_cutstrings.append(muon_preselection)
#cutnames.append('ele skim'); 			  cutstrings.append('ele1_pt>muon1_pt'); 										 prior_cutstrings.append('ele1_pt>muon1_pt')
#cutnames.append('ele preselection'); 	  cutstrings.append(ele_preselection); 											 prior_cutstrings.append('ele1_pt>muon1_pt')
#cutnames.append('ele kinematics'); 		  cutstrings.append(ele_preselection+' && '+ele_kinematics); 					 prior_cutstrings.append(ele_preselection)
#cutnames.append('ele ID'); 				  cutstrings.append(ele_preselection+' && '+ele_ID); 							 prior_cutstrings.append(ele_preselection)
#cutnames.append('ele 2D cut'); 			  cutstrings.append(ele_preselection+' && '+ele_2D); 							 prior_cutstrings.append(ele_preselection)
#cutnames.append('ele leptonic top mass'); cutstrings.append(ele_preselection+' && '+lep_top_mass); 						 prior_cutstrings.append(ele_preselection)
#cutnames.append('ele full leptonic'); 	  cutstrings.append(ele_full_leptonic); 										 prior_cutstrings.append(ele_preselection)
#cutnames.append('ele hadronic pretag');   cutstrings.append(ele_hadronic_pretag); 										 prior_cutstrings.append(ele_preselection)
#cutnames.append('ele signal mass'); 	  cutstrings.append(ele_hadronic_pretag+' && '+signal_mass); 					 prior_cutstrings.append(ele_hadronic_pretag)
#cutnames.append('ele signal tau32'); 	  cutstrings.append(ele_hadronic_pretag+' && '+signal_tau32); 					 prior_cutstrings.append(ele_hadronic_pretag)
#cutnames.append('ele full selection');    cutstrings.append(ele_hadronic_pretag+' && '+signal_mass+' && '+signal_tau32); prior_cutstrings.append(ele_preselection)

data_events_at_cut = []; data_events_at_prior_cut = []
events_at_cut = []; events_at_prior_cut = []
dist = TH1D('dist','distribution',20,-1.0,1.0)
for i in range(len(cutnames)) :
	print 'Getting numbers of data events for cut '+cutnames[i]+' ('+str(i+1)+' out of '+str(len(cutnames))+')'
	if 'muon' in cutnames[i] :
		tmp = dist.Clone('tmp')
		muon_data_chain.CopyTree(cutstrings[i]).Draw('cstar>>tmp','weight')
		data_events_at_cut.append(tmp.Integral())
		tmp = dist.Clone('tmp')
		muon_data_chain.CopyTree(prior_cutstrings[i]).Draw('cstar>>tmp','weight')
		data_events_at_prior_cut.append(tmp.Integral())
	elif 'ele' in cutnames[i] :
		tmp = dist.Clone('tmp')
		ele_data_chain.CopyTree(cutstrings[i]).Draw('cstar>>tmp','weight')
		data_events_at_cut.append(tmp.Integral())
		tmp = dist.Clone('tmp')
		ele_data_chain.CopyTree(prior_cutstrings[i]).Draw('cstar>>tmp','weight')
		data_events_at_prior_cut.append(tmp.Integral())
	events_at_cut.append([]); events_at_prior_cut.append([])
	print 'Getting numbers of MC events for cut '+cutnames[i]+' ('+str(i+1)+' out of '+str(len(cutnames))+')'
	for j in range(len(MC_chains)) :
		print '	Doing '+shortnames_done[j]+' ('+str(j+1)+' out of '+str(len(MC_chains))+')'
		tmp = dist.Clone('tmp')
		MC_chains[j].CopyTree(cutstrings[i]).Draw('cstar>>tmp','weight*(21560109./245.8)')
		events_at_cut[i].append(tmp.Integral())
		tmp = dist.Clone('tmp')
		MC_chains[j].CopyTree(prior_cutstrings[i]).Draw('cstar>>tmp','weight*(21560109./245.8)')
		events_at_prior_cut[i].append(tmp.Integral())

#print out the number of events in data and the efficiencies for the MC samples
#first line is just table headings for each cutflow and each sample type
first_line = 'Cut 		 	Data Events 		'
for shortname in shortnames_done :
	first_line += shortname+' eff		'
print first_line
os.system('echo "'+first_line+'" > '+cutflow_filename+'')
#second line is just the total number of events in data
second_line = '	 	'+str(data_events_at_cut[0])+' 			'
for shortname in shortnames_done :
	second_line += '		'
print second_line
os.system('echo "'+second_line+'" >> '+cutflow_filename+'')
#Each line after that is the cutflow number, then the number of events in data, then the eff. with uncertainty for each sample type
for i in range(1,len(cutnames)) :
	#Begin with the cutflow number and the number of events in data
	next_line = ''
	next_line += cutnames[i]
	next_line+=' 	%.2f'%(data_events_at_cut[i])
	for j in range(len(MC_chains)) :
		#calculate the efficiency for this sample type
		eff = events_at_cut[i][j]/events_at_prior_cut[i][j]
		eff_err = eff*sqrt(1./events_at_cut[i][j]+1./events_at_prior_cut[i][j])
		#add to the line
		next_line+='		%.6f (%.6f)'%(eff,eff_err)
	print next_line
	os.system('echo "'+next_line+'" >> '+cutflow_filename+'')