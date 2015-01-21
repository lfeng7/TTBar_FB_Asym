#makes some plots based on information in the nTuples and genParticles and stuff
#mostly for dealing with ttbar MC or whatever I guess

import ROOT
import glob
from optparse import OptionParser
from DataFormats.FWLite import Events, Handle
from math import *

# COMMAND LINE OPTIONS
parser = OptionParser()
parser.add_option('--dir', metavar='F', type='string', action='store', dest='directory', help='') ## Sets which files to run on
parser.add_option('--generator', metavar='F', type='string', action='store', default='powheg', dest='generator', help='') ## MC generator: 'powheg' (default),
																													## 'madgraph','mcatnlo', or 'pythia8'
parser.add_option('--maxEvents', metavar='F', type='int', action='store', default='100000', dest='maxEvents', help='') ## hard cut on number of events (100k default)
parser.add_option('--out', metavar='F', type='string', action='store', default='NEW_NTUPLE_PLOTS', dest='out', help='') ## name of output file that holds plots
(options, args) = parser.parse_args()

#Read in input files
FILES = options.directory + "/*.root"
files = glob.glob(FILES)
events = Events(files)

#Set up handles and labels for use
#MC GenParticle variables
GenHandle = Handle( "vector<reco::GenParticle>" )
GenLabel = ( "prunedGenParticles", "" )
#pythia8 nTuple GenParticles
genPartHandles = [];	genPartLabels = []
genPartLabels.append(('genPart','genPartPt')); 	   genPartHandles.append(Handle('vector<float>'))
genPartLabels.append(('genPart','genPartEta'));    genPartHandles.append(Handle('vector<float>'))
genPartLabels.append(('genPart','genPartPhi'));    genPartHandles.append(Handle('vector<float>'))
genPartLabels.append(('genPart','genPartMass'));   genPartHandles.append(Handle('vector<float>'))
genPartLabels.append(('genPart','genPartID')); 	   genPartHandles.append(Handle('vector<float>'))
genPartLabels.append(('genPart','genPartMomID'));  genPartHandles.append(Handle('vector<float>'))
genPartLabels.append(('genPart','genPartStatus')); genPartHandles.append(Handle('vector<float>'))

#formatted particle names based on pdgIds
ids = ( [(2212,'p'),(1,'d'),(2,'u'),(3,'s'),(4,'c'),(5,'b'),(6,'t'),(11,'e'),(12,'ve'),(13,'mu'),(14,'vmu'),(15,'tau'),(16,'vtau'),(21,'g'),
		 (22,'photon'),(23,'Z'),(24,'W'),(310,'K_S'),(513,'B*0'),(511,'B0'),(413,'D*(2010)+'),(421,'D0'),(411,'D+'),(111,'pi0'),(521,'B+')] )
def getId(pdgid) :
	name = ''
	if pdgid < 0 :
		name = name+'~'
		pdgid = -1*pdgid
	for tup in ids :
		if tup[0] == pdgid :
			name = name+tup[1]
	if name == '~' or name == '' :
		name = name + str(pdgid)
	return name

#histograms
all_histos = []
lep_b_pT_boosted_tops = ROOT.TH1F('lep_b_pT_boosted_tops','p_{T} of leptonic b-jet in boosted semileptonic t#bar{t} events; p_{T} (GeV)',50,0.,250.) 
all_histos.append(lep_b_pT_boosted_tops)
lep_b_dR_boosted_tops = ROOT.TH1F('lep_b_dR_boosted_tops','#Delta R(lepton,leptonic b-jet) in boosted semileptonic t#bar{t} events; #Delta R',60,0.,3.) 
all_histos.append(lep_b_dR_boosted_tops)

count = 0
for event in events :
	count+=1
	#hard cut on max events
	if count==options.maxEvents+1 :
		print 'ITERATION '+str(count)+' REACHED, EXITING'
		break
	#print progress every, like, 5000 events or whatever
	if (count-1)%5000 == 0 :
		print 'At iteration '+str(count)+''
	#open the GenParticles structure
	event.getByLabel( GenLabel, GenHandle )
	if not GenHandle.isValid() :
		print 'INVALID GENPARTICLE HANDLE, EXITING'
		break
	GenParticles = GenHandle.product()
	#attributes of the event
	is_semilep = False
	is_boosted = False
	#things we need
	#look for the t, tbar, and leptonic b
	t_vec = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
	tbar_vec = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
	lep_vec = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
	lep_b_vec = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
	#count the number of leptonic W decay products
	n_leps_from_Ws = 0
#	print 'New Event -----------------------------------------' #DEBUGGING
	#loop over all the genPartticles
	for ig in GenParticles :
		#Powheg algorithm
		if options.generator.lower() == 'powheg' :
			if fabs(ig.pdgId()) == 6 and ig.status() == 3 :
#				print 'found a top' #DEBUGGING
#				print ('ig = ('+str(ig.pt())+','+str(ig.eta())+','+str(ig.phi())+','+str(ig.mass())+		  #DEBUGGING
#							') = ('+str(ig.px())+','+str(ig.py())+','+str(ig.pz())+','+str(ig.energy())+')' ) #DEBUGGING
				#found a top
				if ig.pdgId() == 6 : 
					t_vec.SetPtEtaPhiM(ig.pt(),ig.eta(),ig.phi(),ig.mass())
				elif ig.pdgId() == -6 : 
					tbar_vec.SetPtEtaPhiM(ig.pt(),ig.eta(),ig.phi(),ig.mass())
				#look through its daughters to see if its W decayed leptonically
				this_W_leptonic = False
				for i in range(ig.numberOfDaughters()) :
					if fabs(ig.daughter(i).pdgId()) == 24 and ig.daughter(i).status()==3 :
#						print 'found its W' #DEBUGGING
						dauW = ig.daughter(i)
						for j in range(dauW.numberOfDaughters()) :
							if fabs(dauW.daughter(j).pdgId()) in range(11,17) and dauW.daughter(j).status()==3 :
								n_leps_from_Ws+=1
								this_W_leptonic = True
								lep = dauW.daughter(j)
								if fabs(lep.pdgId()) == 11 or fabs(lep.pdgId()) == 13 or fabs(lep.pdgId()) == 15 :
									lep_vec.SetPtEtaPhiM(lep.pt(),lep.eta(),lep.phi(),lep.mass())
#						print 'this_W_leptonic = '+str(this_W_leptonic)+'' #DEBUGGING
				#if its W decayed leptonically, the b is the leptonic b
				if this_W_leptonic :
					for i in range(ig.numberOfDaughters()) :
						if fabs(ig.daughter(i).pdgId()) == 5 and ig.daughter(i).status()==3 :
#							print 'found its leptonic b' #DEBUGGING
							daub = ig.daughter(i)
							lep_b_vec.SetPtEtaPhiM(daub.pt(),daub.eta(),daub.phi(),daub.mass())
		#pythia8 algorithm
		elif options.generator.lower() == 'pythia8' :
			if fabs(ig.pdgId()) == 6 and ig.status() == 22 and ig.numberOfMothers()==2 and ig.mother(0).pdgId()==2212 and ig.mother(1).pdgId()==2212 :
#				print 'found a top' #DEBUGGING
				#found a top
				#crank the particle down through its daughters because pythia8. . .
				thist = ig
				while thist.daughter(0).pdgId() == thist.pdgId() :
#					print '	top crank' #DEBUGGING
					thist = thist.daughter(0)
				if thist.pdgId() == 6 : 
					t_vec.SetPtEtaPhiM(thist.pt(),thist.eta(),thist.phi(),thist.mass())
				elif thist.pdgId() == -6 : 
					tbar_vec.SetPtEtaPhiM(thist.pt(),thist.eta(),thist.phi(),thist.mass())
				#look through its daughters to see if its W decayed leptonically
				this_W_leptonic = False
				for i in range(thist.numberOfDaughters()) :
					if fabs(thist.daughter(i).pdgId()) == 24 :
#						print 'found its W' #DEBUGGING
						dauW = thist.daughter(i)
						#crank the W down as well
						while dauW.daughter(0).pdgId() == dauW.pdgId() :
#							print '	W crank' #DEBUGGING
							dauW = dauW.daughter(0)
						for j in range(dauW.numberOfDaughters()) :
							if fabs(dauW.daughter(j).pdgId()) in range(11,17) :
								n_leps_from_Ws+=1
								this_W_leptonic = True
								lep = dauW.daughter(j)
								#crank down the lepton
								while lep.numberOfDaughters()>0 and lep.daughter(0).pdgId()==lep.pdgId() :
									lep = lep.daughter(0)
								if fabs(lep.pdgId()) == 11 or fabs(lep.pdgId()) == 13 or fabs(lep.pdgId()) == 15 :
									lep_vec.SetPtEtaPhiM(lep.pt(),lep.eta(),lep.phi(),lep.mass())
#						print 'this_W_leptonic = '+str(this_W_leptonic)+'' #DEBUGGING
				#if its W decayed leptonically, the b is the leptonic b
				if this_W_leptonic :
					for i in range(thist.numberOfDaughters()) :
						if fabs(thist.daughter(i).pdgId()) == 5 :
#							print 'found its leptonic b' #DEBUGGING
							daub = thist.daughter(i)
							lep_b_vec.SetPtEtaPhiM(daub.pt(),daub.eta(),daub.phi(),daub.mass())
		
	#If the ttbar mass is greater than, like, 700, then the event is pretty darn boosted.
#	print 't = ('+str(t_vec.Px())+','+str(t_vec.Py())+','+str(t_vec.Pz())+','+str(t_vec.E())+')' #DEBUGGING
#	print 'tbar = ('+str(tbar_vec.Px())+','+str(tbar_vec.Py())+','+str(tbar_vec.Pz())+','+str(tbar_vec.E())+')' #DEBUGGING
#	print 'ttbar_mass = '+str((t_vec+tbar_vec).M())+'' #DEBUGGING
	is_boosted = (t_vec+tbar_vec).M()>700.
	#If the number of leptonic W decay products is 2, then the event was semileptonic
#	print 'n_leps_from_Ws = '+str(n_leps_from_Ws)+'' #DEBUGGING
	is_semilep = n_leps_from_Ws==2

	#add to histograms
	if is_boosted and is_semilep :
		lep_b_pT_boosted_tops.Fill(lep_b_vec.Pt())
		lep_b_dR_boosted_tops.Fill(lep_vec.DeltaR(lep_b_vec))
	
#	genPartVars = []
#	for i in range(len(genPartHandles)) :
#		event.getByLabel(genPartLabels[i],genPartHandles[i])
#		if not genPartHandles[i].isValid() :
#			print 'invalid handle'
#			break
#		genPartVars.append(genPartHandles[i].product())
#	print 'New Event ------------------------------------------------------------'
#	for i in range(len(genPartVars[0])) :
#		if genPartVars[5][i] == genPartVars[4][i] :
#			continue
#		#if genPartVars[6][i] != 21 :
#		#	continue
#		print '{'+getId(genPartVars[5][i])+'} -> ('+getId(genPartVars[4][i])+')	status = '+str(genPartVars[6][i])+' (%.3f,%.3f,%.3f,%.3f)'%(genPartVars[0][i],genPartVars[1][i],genPartVars[2][i],genPartVars[3][i])

#setup the output file and write histograms
outfile = ROOT.TFile(options.out+'.root','recreate')
outfile.cd()
for histo in all_histos :
	histo.Write()
outfile.Write()
outfile.Close()


