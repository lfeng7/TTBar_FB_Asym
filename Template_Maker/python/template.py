#imports
from ROOT import *
from array import array

#global variables
#histogram limits
XBINS = array('d',[-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
YBINS = array('d',[0.,0.075,0.15,0.3,0.7])
ZBINS = array('d',[500.,700.,900.,1100.,1300.,1500.,2500.])

#TDR Style
#gROOT.Macro('rootlogon.C')

##############################		   Template Class  		##############################

class template :
	#docstring
	"""template class"""
	
	#__init__function
	def __init__(self,name,formatted_name) :
		print '				Adding template with name '+name
		self.name = name
		self.formatted_name = formatted_name
		self.histo_3D = TH3D(name,	   formatted_name+'; c*; |x_{F}|; M (GeV)',len(XBINS)-1,XBINS,len(YBINS)-1,YBINS,len(ZBINS)-1,ZBINS)
		self.histo_x  = TH1D(name+'_x',formatted_name+' X Projection; c*',len(XBINS)-1,XBINS)
		self.histo_y  = TH1D(name+'_y',formatted_name+' Y Projection; |x_{F}|',len(YBINS)-1,YBINS)
		self.histo_z  = TH1D(name+'_z',formatted_name+' Z Projection; M (GeV)',len(ZBINS)-1,ZBINS)
		self.histo_3D.SetDirectory(0); self.histo_x.SetDirectory(0); self.histo_y.SetDirectory(0); self.histo_z.SetDirectory(0)

	def Fill(self,c,x,m,w) :
		self.histo_3D.Fill(c,x,m,w)
		self.histo_x.Fill(c,w)
		self.histo_y.Fill(x,w)
		self.histo_z.Fill(m,w)

	#convertTo1D takes a 3D distribution and makes it 1D for use with theta
	def convertTo1D(self) :
		nBins = self.histo_3D.GetNbinsX()*self.histo_3D.GetNbinsY()*self.histo_3D.GetNbinsZ()
		newHisto = TH1F(self.histo_3D.GetName(),self.histo_3D.GetTitle(),nBins,0.,nBins-1.)
		newHisto.SetDirectory(0)
		for k in range(nBins) :
			newHisto.SetBinContent(k,self.histo_3D.GetBinContent(k))
		return newHisto