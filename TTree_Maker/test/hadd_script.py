import os, glob
from ROOT import *

sample_names = []
sample_names.append('Powheg_qq_semilep_TT')
sample_names.append('Powheg_qq_semilep_TT_SC')
sample_names.append('Powheg_qq_semilep_TT_Mtt_700_to_1000')
sample_names.append('Powheg_qq_semilep_TT_Mtt_1000_to_Inf')
sample_names.append('Powheg_gg_semilep_TT')
sample_names.append('Powheg_gg_semilep_TT_SC')
sample_names.append('Powheg_gg_semilep_TT_Mtt_700_to_1000')
sample_names.append('Powheg_gg_semilep_TT_Mtt_1000_to_Inf')
sample_names.append('Powheg_dilep_TT')
sample_names.append('Powheg_dilep_TT_SC')
sample_names.append('Powheg_dilep_TT_Mtt_700_to_1000')
sample_names.append('Powheg_dilep_TT_Mtt_1000_to_Inf')
sample_names.append('Powheg_had_TT')
sample_names.append('Powheg_had_TT_SC')
sample_names.append('Powheg_had_TT_Mtt_700_to_1000')
sample_names.append('Powheg_had_TT_Mtt_1000_to_Inf')
sample_names.append('W1Jets')
sample_names.append('W2Jets')
sample_names.append('W3Jets')
sample_names.append('W4Jets')
sample_names.append('DY1Jets')
sample_names.append('DY2Jets')
sample_names.append('DY3Jets')
sample_names.append('DY4Jets')
sample_names.append('T_s')
sample_names.append('T_t')
sample_names.append('T_tW')
sample_names.append('Tbar_s')
sample_names.append('Tbar_t')
sample_names.append('Tbar_tW')
sample_names.append('SingleMu_Run2012A')
sample_names.append('SingleMu_Run2012B')
#sample_names.append('SingleMu_Run2012C')
sample_names.append('SingleMu_Run2012D')
sample_names.append('SingleEl_Run2012A')
sample_names.append('SingleEl_Run2012B')
#sample_names.append('SingleEl_Run2012C')
sample_names.append('SingleEl_Run2012D')

for name in sample_names :
    print 'doing '+name
    
#    #make directories
#    os.system('mkdir '+name)
    
    os.chdir(name)
    
#    #copy scripts
#    os.system('cp ../grid_sub.csh .; cp ../cleanup.bash .')
    
#    #make input file
#    directory = raw_input('nTuple directory for '+name+': ')
#    os.system('python ../make_ttree_input_file.py --directory '+directory)
    
#    #make ana.listOfJobs
#    eventType = 'none'
#    if 'qq_semilep' in name :
#        eventType = 'qq_semilep'
#    elif 'gg_semilep' in name :
#        eventType = 'gg_semilep'
#    elif 'dilep' in name :
#        eventType = 'dilep'
#    elif 'had' in name :
#        eventType = 'had'
#    nJobs = raw_input('number of jobs for '+name+': ')
#    generator = raw_input('MC generator for '+name+': ')
#    crossSection = raw_input('cross section for '+name+': ')
#    nEvents = raw_input('number of events for '+name+': ')
#    cmd = 'python ../make_list_of_jobs.py --event_type '+eventType+' --n_jobs '+nJobs+' --name '+name
#    cmd+= ' --generator '+generator+' --cross_section '+crossSection+' --n_events '+nEvents
#    if 'Run2012' in name :
#        cmd+=' --data yes'
#    os.system(cmd)
    
#    #submit jobs
#    os.system('tcsh grid_sub.csh')
    
    #skim files
    filelist = glob.glob('*_tree.root')
    for i in range(len(filelist)) :
        print ' '+str(i)+''
        f = TFile(filelist[i]); t = f.Get('tree')
        newTree = t.CopyTree('hadt_M>100.')
        name = filelist[i].replace('_tree.root','')+'_skim_tree.root'
        newFile = TFile(name,'recreate')
        newTree.Write()

    #hadd files
    #os.system('hadd '+name+'_skim_all.root '+name+'*skim_tree.root')
    #os.system('mv *_all.root ../total_ttree_files')
    os.chdir('..')
#os.chdir('total_ttree_files')
#os.system('hadd SingleMu_Run2012_all.root SingleMu_Run2012*_all.root')
#cmd = 'hadd Powheg_semilep_TT_all.root Powheg_qq_semilep_TT_all.root'
#cmd += ' Powheg_qq_semilep_TT_SC_all.root Powheg_gg_semilep_TT_all.root Powheg_gg_semilep_TT_SC_all.root'
#os.system(cmd)
#os.chdir('..')
