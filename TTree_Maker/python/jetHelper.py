import ROOT

#jetHelper contains helper functions for jet selection

#selectJets takes in the list of jet variables
#returns:
#	1) a list of selected jet indices for the analysis, ordered by pT OR
#	2) 1-element list of negative value of a cutflow failure point
def selectJets(jetvars) :
	return [0,1,2,3] #DEBUG RETURN