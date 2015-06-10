#lepHelper: contains helper functions for lepton selection
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

import ROOT
from math import *

#muonCuts puts high quality muons above high quality electrons in the list
def muonCuts(mulist) :
	return mulist #DEBUG RETURN



#findNearestJet returns the jet object of the jet closest to this lepton
def findNearestJet(lepvec,jets_list) :
	closestDR = jets_list[0].vec.DeltaR(lepvec)
	closestJet = jets_list[0]
	for jet in jets_list :
		checkDR = jet.vec.DeltaR(lepvec)
		if checkDR < closestDR :
			closestDR = checkDR
			closestJet = jet
	return closestJet


