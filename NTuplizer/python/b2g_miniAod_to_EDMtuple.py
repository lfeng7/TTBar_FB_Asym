import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("B2G")

#######################################################
#######             Parameters             ############
#######################################################
# Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
options = VarParsing ('python')
options.register('ignoreTrigger',
                 1,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Ignore trigger in selection")
options.register('muOrEle',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Use muons (0) or electrons (1)")
options.register('useData',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Use data (1) or MC (0)")
options.parseArguments()
print options

#######################################################
#######             	I/O 	           ############
#######################################################

# Input file
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
	'/store/mc/Phys14DR/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/E873348E-BC70-E411-BFA8-0025907B4FD6.root'
))
# Max events
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(3000))
# Message Service
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
# Output file
process.out = cms.OutputModule("PoolOutputModule",
		   fileName = cms.untracked.string("ntuple.root"),
		   SelectEvents   = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
		   outputCommands = cms.untracked.vstring('drop *',
                     'keep *_*_*_B2G',
                     'keep *_generator_*_*',
                     'keep *_*prunedGenParticles*_*_*',
                     'keep *_pileup*_*_*',
                     'keep *_pdfWeights*_*_*'
))
process.outpath = cms.EndPath(process.out)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.out.dropMetaData = cms.untracked.string("DROPPED")

#######################################################
#######        Extra Jet Collections 	   ############
#######################################################

from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
from RecoJets.JetProducers.caTopTaggers_cff import cmsTopTagPFJetsCHS
from RecoJets.JetProducers.caTopTaggers_cff import caTopTagInfos
process.chs = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV()>1"))
#AK4 Jets, CHS
process.ak4PFJetsCHS = ak5PFJets.clone(
					src = 'chs',
					rParam = cms.double(0.4),
					jetPtMin = cms.double(10.0),
					jetAlgorithm = cms.string("AntiKt"))
#CA8 Jets, CHS
process.ca8PFJetsCHS = ak5PFJets.clone(
					src = 'chs',
					rParam = cms.double(0.8),
					jetPtMin = cms.double(3.0),
					jetAlgorithm = cms.string("CambridgeAachen"))
#CA8 Jets, CHS, Pruned
process.ca8PFJetsCHSPruned = ak5PFJets.clone(
					src = 'chs',
					rParam = cms.double(0.8),
					jetPtMin = cms.double(3.0),
					jetAlgorithm = cms.string("CambridgeAachen"),
					usePruning = cms.bool(True),
					useExplicitGhosts = cms.bool(True),
					jetCollInstanceName = cms.string("SubJets"),
					writeCompound = cms.bool(True),
					nFilt = cms.int32(2),
					zcut = cms.double(0.1),
					rcut_factor = cms.double(0.5))
#CMS Top Tag Jets, CHS
process.cmsTopTagPFJetsCHS = cmsTopTagPFJetsCHS.clone(src = cms.InputTag('chs'))
#CMS Top Tage Infos
process.caTopTagInfos = caTopTagInfos.clone(src = cms.InputTag("cmsTopTagPFJetsCHS"))
#CA8 NJettiness
process.ca8Njettiness = Njettiness.clone(src = 'ca8PFJetsCHS')
#CA8 Jets, CHS, with selection
process.selectedca8PFJetsCHS = cms.EDFilter('PFJetSelector',
					src = cms.InputTag('ca8PFJetsCHS'),
					cut = cms.string('pt > 25 && abs(eta) < 2.4'))
#CA8 Njettiness, with selection
process.selectedca8Njettiness = Njettiness.clone(src = 'selectedca8PFJetsCHS')
#final extra jet collections sequence
process.extraJetCols = cms.Sequence(
	process.ak4PFJetsCHS *
	process.ca8PFJetsCHS *
	process.ca8PFJetsCHSPruned *
	process.ca8Njettiness *
	process.cmsTopTagPFJetsCHS *
	process.caTopTagInfos
	#process.selectedca8PFJetsCHS *
	#process.selectedca8Njettiness
)

#######################################################
#######        	  Output Branches		   ############
#######################################################

#Trigger
process.trigger = cms.EDFilter('b2g_miniAodAnalyzer_trigger',
					printAll		= cms.bool(False), # Dump all trigger information
					isData			= cms.bool(True), # If false will not save HLT Objects
					bits			= cms.InputTag("TriggerResults", "", "HLT"),
	    			prescales 		= cms.InputTag("patTrigger"),
					objects			= cms.InputTag("selectedPatTrigger"),
					useTriggerList  = cms.bool(False), # Only save information for triggers in the list
					triggerList		= cms.vstring("HLT_Ele25WP60_SC4_Mass55_v1",
								      			  "HLT_Mu40_v1"
								     	  )
)
#General
process.general = cms.EDFilter('b2g_miniAodAnalyzer_general',
					vertecies		 = cms.InputTag("offlineSlimmedPrimaryVertices"),
	    			met 			 = cms.InputTag("slimmedMETs"),
					electrons		 = cms.InputTag("slimmedElectrons"),
		    		muons			 = cms.InputTag("slimmedMuons"),
		    		taus			 = cms.InputTag("slimmedTaus"),
					photons			 = cms.InputTag("slimmedPhotons"),
					ak4slimmed		 = cms.InputTag("slimmedJets"),
					ak8slimmed		 = cms.InputTag("slimmedJetsAK8"),
					ak8grommedMasses = cms.bool(True),
					pfcands			 = cms.InputTag("packedPFCandidates")
)
#Add extra jet collections running FastJet 'by hand'
#CA8 FastJet Jets
process.ca8jets = cms.EDFilter('b2g_miniAodAnalyzer_jets',
					vertecies	= cms.InputTag("offlineSlimmedPrimaryVertices"),
					pfcands		= cms.InputTag("packedPFCandidates"),
					jetAlgo  	= cms.string("CA"), # CA, KT, AK
					jetR  		= cms.double(0.8),
					jetPtmin  	= cms.double(10.0)
)
#AK4 FastJet Jets
process.ak4jets = cms.EDFilter('b2g_miniAodAnalyzer_jets',
					vertecies	= cms.InputTag("offlineSlimmedPrimaryVertices"),
					pfcands		= cms.InputTag("packedPFCandidates"),
					jetAlgo  	= cms.string("AK"), # CA, KT, AK
					jetR  		= cms.double(0.4),
					jetPtmin  	= cms.double(10.0)
)
#AK8 FastJet Jets
process.ak8jets = cms.EDFilter('b2g_miniAodAnalyzer_jets',
					vertecies	= cms.InputTag("offlineSlimmedPrimaryVertices"),
					pfcands		= cms.InputTag("packedPFCandidates"),
					jetAlgo  	= cms.string("AK"), # CA, KT, AK
					jetR  		= cms.double(0.8),
					jetPtmin  	= cms.double(100.0)
)

#######################################################
#######             PDF Weights            ############
#######################################################

process.pdfWeights = cms.EDProducer("PdfWeightProducer",
      FixPOWHEG   = cms.untracked.string(""),
      GenTag 	  = cms.untracked.InputTag("prunedGenParticles"),
      PdfInfoTag  = cms.untracked.InputTag("generator"),
      PdfSetNames = cms.untracked.vstring( 
                                             "CT10.LHgrid",   # This is for POWHEG with tune CT10 and LO, which is the tune for /TT_CT10_TuneZ2star_8TeV-powheg-tauola/
                                             "GJR08VFnloE.LHgrid", #This is for GJR
                                             "cteq66.LHgrid",  # This is for Madgraph sample with tune cteq66
                                            # "CT10nlo.LHgrid" # This is for POWHEG with tune CT10 and NLO

        )
    )

#Final Output Path
process.p = cms.Path(
	process.chs *
	process.extraJetCols *
	process.trigger *
	process.general *
	process.ca8jets *
	process.ak4jets *
	process.pdfWeights*
	process.ak8jets
)
print "---+++---+++---+++---"
