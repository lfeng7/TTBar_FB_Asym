from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'TT_Tune4C_13TeV-pythia8-tauola'

config.section_('JobType')
config.JobType.psetName = '/uscms_data/d3/eminizer/ttbar_run2/CMSSW_7_2_0/src/Analysis/NTuplizer/test/b2gedmntuples_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.pyCfgParams = ['maxEvts=100', 'LHE=False']

config.section_('Data')
config.Data.inputDataset = '/TT_Tune4C_13TeV-pythia8-tauola/Phys14DR-PU40bx25_tsg_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.unitsPerJob = 1
config.Data.inputDBS = 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet'
config.Data.splitting = 'FileBased'
config.Data.publishDataName = 'TT_Tune4C_13TeV-pythia8-tauola'

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T3_US_FNALLPC'
