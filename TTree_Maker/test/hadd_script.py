import os, glob
from ROOT import *

sample_names = []
#sample_names.append('Powheg_qq_semilep_TT')
#sample_names.append('Powheg_gg_semilep_TT')
#sample_names.append('Powheg_dilep_TT')
#sample_names.append('Powheg_had_TT')
#sample_names.append('DY1Jets')
#sample_names.append('DY2Jets')
#sample_names.append('DY3Jets')
#sample_names.append('DY4Jets')
#sample_names.append('T_s')
#sample_names.append('T_t')
#sample_names.append('T_tW')
#sample_names.append('Tbar_s')
#sample_names.append('Tbar_t')
#sample_names.append('Tbar_tW')
#sample_names.append('SingleMu_Run2012A')
#sample_names.append('SingleMu_Run2012B')
sample_names.append('SingleMu_Run2012C')
#sample_names.append('SingleMu_Run2012D')
#sample_names.append('SingleEl_Run2012A')
#sample_names.append('SingleEl_Run2012B')
#sample_names.append('SingleEl_Run2012C')
#sample_names.append('SingleEl_Run2012D')

for name in sample_names :
    print 'doing '+name
    
#    #make directories
#    os.system('mkdir '+name)
    
    os.chdir(name)

#    #get rid of old files
#    os.system('bash cleanup.bash')
#    os.system('rm -rf output *.root')
#    os.system('mv ana.listOfJobs_all ana.listOfJobs')
    
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
 
#    #make list of failed jobs
#    os.system('python ../make_failed_job_list.py')

    #skim files
    filelist = glob.glob('*_tree.root')
    for i in range(len(filelist)) :
        print ' '+str(i)+''
        f = TFile(filelist[i]); t = f.Get('tree')
        newname = filelist[i].replace('_tree.root','')+'_skim_tree.root'
        newFile = TFile(newname,'recreate')
        newTree = t.CopyTree('hadt_M>100.')
        newTree.Write()
        newFile.Close()
    i = 0
    while i<len(filelist) :
        if filelist[i].find('JES')!=-1 or filelist[i].find('JER')!=-1 or filelist[i].find('skim')!=-1 :
            filelist.pop(i)
        else :
            i+=1
    cmd = 'hadd -f '+name+'_skim_all.root '+name+'_?_skim_tree.root'
    if len(filelist) > 10 :
        cmd += ' '+name+'_??_skim_tree.root'
        if len(filelist) > 100 :
            cmd += ' '+name+'_???_skim_tree.root'
            if len(filelist) > 1000 :
                cmd += ' '+name+'_????_skim_tree.root'
    os.system(cmd)
    if name.find('Run2012')==-1 :
        cmd = 'hadd -f '+name+'_JES_up_skim_all.root '+name+'_JES_up_*_skim_tree.root'
        os.system(cmd)
        cmd = 'hadd -f '+name+'_JES_down_skim_all.root '+name+'_JES_down_*_skim_tree.root'
        os.system(cmd)
        cmd = 'hadd -f '+name+'_JER_up_skim_all.root '+name+'_JER_up_*_skim_tree.root'
        os.system(cmd)
        cmd = 'hadd -f '+name+'_JER_down_skim_all.root '+name+'_JER_down_*_skim_tree.root'
        os.system(cmd)
    os.system('mv *_all.root ../total_ttree_files')

    #skim files (e/mu selection)
    filelist = glob.glob('*_tree.root')
    i = 0
    while i<len(filelist) :
        if filelist[i].find('JES')!=-1 or filelist[i].find('JER')!=-1 or filelist[i].find('skim')!=-1 :
            filelist.pop(i)
        else :
            i+=1
    for i in range(len(filelist)) :
        print ' '+str(i)+''
        f = TFile(filelist[i]); t = f.Get('tree')
        cutstring = 'mu_trigger==1 && muon1_isLoose==1 && muon1_pt>40. && abs(muon1_eta)<2.4'
        cutstring+= ' && (muon1_relPt>25. || muon1_dR>0.5)'
        cutstring+= ' && ele1_pt>40. && abs(ele1_eta)<2.4'
        cutstring+= ' && (ele1_relPt>25. || ele1_dR>0.5)'
        newTree = t.CopyTree(cutstring)
        newname = filelist[i].replace('_tree.root','')+'_emu_skim_tree.root'
        newFile = TFile(newname,'recreate')
        newTree.Write()
        newFile.Close()
    cmd = 'hadd -f '+name+'_emu_skim_all.root '+name+'_?_emu_skim_tree.root'
    if len(filelist) > 10 :
        cmd += ' '+name+'_??_emu_skim_tree.root'
        if len(filelist) > 100 :
            cmd += ' '+name+'_???_emu_skim_tree.root'
            if len(filelist) > 1000 :
                cmd += ' '+name+'_????_emu_skim_tree.root'
    os.system(cmd)
    os.system('mv *_all.root ../total_ttree_files')
    os.system('rm -rf *emu_skim*.root')

    #skim files (e+/e- selection)
    filelist = glob.glob('*_tree.root')
    i = 0
    while i<len(filelist) :
        if filelist[i].find('JES')!=-1 or filelist[i].find('JER')!=-1 or filelist[i].find('skim')!=-1 :
            filelist.pop(i)
        else :
            i+=1
    for i in range(len(filelist)) :
        print ' '+str(i)+''
        f = TFile(filelist[i]); t = f.Get('tree')
        newname = filelist[i].replace('_tree.root','')+'_ee_skim_tree.root'
        newFile = TFile(newname,'recreate')
        ele1_tag_selec = '(ele1_pt>40. && abs(ele1_eta)<2.4 && ele1_isLoose==1 && (ele1_relPt>25. || ele1_dR>0.5))'
        ele1_probe_selec = '(ele1_pt>40. && abs(ele1_eta)<2.4 && (ele1_relPt>25. || ele1_dR>0.5))'
        ele2_tag_selec = '(ele2_pt>40. && abs(ele2_eta)<2.4 && ele2_isLoose==1 && (ele2_relPt>25. || ele2_dR>0.5))'
        ele2_probe_selec = '(ele2_pt>40. && abs(ele2_eta)<2.4 && (ele2_relPt>25. || ele2_dR>0.5))'
        cutstring = '(('+ele1_tag_selec+' && '+ele2_probe_selec+') || ('+ele2_tag_selec+' && '+ele1_probe_selec+'))'
        #muon 1 and 2 rejection
        cutstring += ' && (muon1_pt<40. || abs(muon1_eta)>2.4 || muon1_isLoose!=1 || (muon1_relPt<25. && muon1_dR<0.5))'
        cutstring += ' && (muon2_pt<40. || abs(muon2_eta)>2.4 || muon2_isLoose!=1 || (muon2_relPt<25. && muon2_dR<0.5))'
        #require two electrons to have opposite charges
        cutstring += ' && (ele1_Q+ele2_Q==0)'
        #other leptonic side cuts
        cutstring+=' && lepW_pt[0]>50.'
        #other jet requirements
        cutstring+=' && lepb_pt>50. && hadt_pt>50. && max(lepb_pt,hadt_pt)>150.'
        cutstring+=' && (lepb_csv>0.244 || hadt_csv>0.244)'
        newTree = t.CopyTree(cutstring)
        newTree.Write()
        newFile.Close()
    cmd = 'hadd -f '+name+'_ee_skim_all.root '+name+'_?_ee_skim_tree.root'
    if len(filelist) > 10 :
        cmd += ' '+name+'_??_ee_skim_tree.root'
        if len(filelist) > 100 :
            cmd += ' '+name+'_???_ee_skim_tree.root'
            if len(filelist) > 1000 :
                cmd += ' '+name+'_????_ee_skim_tree.root'
    os.system(cmd)
    os.system('mv *_all.root ../total_ttree_files')
    os.system('rm -rf *_ee_skim*.root')

    os.chdir('..')
#os.chdir('total_ttree_files')
#os.system('hadd -f SingleMu_Run2012_all.root SingleMu_Run2012*_all.root')
#cmd = 'hadd -f Powheg_semilep_TT_all.root Powheg_qq_semilep_TT_all.root'
#cmd += ' Powheg_qq_semilep_TT_SC_all.root Powheg_gg_semilep_TT_all.root Powheg_gg_semilep_TT_SC_all.root'
#os.system(cmd)
#os.chdir('..')
