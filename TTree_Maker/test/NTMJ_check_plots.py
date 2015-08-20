from ROOT import *
from array import array

LUMINOSITY = 19748.

inputfile = TFile('NTMJ_skim_all.root')

tree = inputfile.Get('tree')

allplots = []; allplotnames = []
x_pass  = TH1D('x_pass','NTMJ background template morphing, c* projection; c*',20,-1.0,1.0); allplots.append(x_pass); allplotnames.append('x_pass')
x_fail  = TH1D('x_fail','NTMJ background template morphing, c* projection; c*',20,-1.0,1.0); allplots.append(x_fail); allplotnames.append('x_fail')
x_morph = TH1D('x_morph','NTMJ background template morphing, c* projection; c*',20,-1.0,1.0); allplots.append(x_morph); allplotnames.append('x_morph')
y_pass  = TH1D('y_pass','NTMJ background template morphing, x_{F} projection; |x_{F}|',30,0.0,0.6); allplots.append(y_pass); allplotnames.append('y_pass')
y_fail  = TH1D('y_fail','NTMJ background template morphing, x_{F} projection; |x_{F}|',30,0.0,0.6); allplots.append(y_fail); allplotnames.append('y_fail')
y_morph = TH1D('y_morph','NTMJ background template morphing, x_{F} projection; |x_{F}|',30,0.0,0.6); allplots.append(y_morph); allplotnames.append('y_morph')
z_pass  = TH1D('z_pass','NTMJ background template morphing, M projection; M (GeV)',20,500.,2500.); allplots.append(z_pass); allplotnames.append('z_pass')
z_fail  = TH1D('z_fail','NTMJ background template morphing, M projection; M (GeV)',20,500.,2500.); allplots.append(z_fail); allplotnames.append('z_fail')
z_morph = TH1D('z_morph','NTMJ background template morphing, M projection; M (GeV)',20,500.,2500.); allplots.append(z_morph); allplotnames.append('z_morph')

for i in range(len(allplots)) :
	name = allplotnames[i]
	plot = allplots[i]
	plot.SetLineWidth(4)
	plot.SetMarkerStyle(21)
	plot.SetFillStyle(0)
	plot.SetLineStyle(1)
	plot.SetMinimum(0.0)
	if 'pass' in name :
		plot.SetLineColor(kBlack)
	elif 'fail' in name :
		plot.SetLineColor(kRed)
	elif 'morph' in name :
		plot.SetLineColor(kBlue)

leg = TLegend(0.62,0.67,0.9,0.9)
leg.AddEntry(x_pass,'passing','F')
leg.AddEntry(x_fail,'failing','F')
leg.AddEntry(x_morph,'morphed','F')

x_canv = TCanvas('x_canv','x_canv',1200,900)
y_canv = TCanvas('y_canv','y_canv',1200,900)
z_canv = TCanvas('z_canv','z_canv',1200,900)

cstar 	   = array('d',[1.0]); tree.SetBranchAddress('cstar_scaled',cstar)
x_F 	   = array('d',[1.0]); tree.SetBranchAddress('x_F_scaled',x_F)
M 		   = array('d',[1.0]); tree.SetBranchAddress('M_scaled',M)
hadt_M 	   = array('d',[1.0]); tree.SetBranchAddress('hadt_M',hadt_M)
hadt_tau32 = array('d',[1.0]); tree.SetBranchAddress('hadt_tau32',hadt_tau32)
weight 	   = array('d',[1.0]); tree.SetBranchAddress('weight',weight)
sf_pileup  = array('d',[1.0]); tree.SetBranchAddress('sf_pileup',sf_pileup)
sf_top_pT  = array('d',[1.0]); tree.SetBranchAddress('sf_top_pT',sf_top_pT)

slope = 0.00131510840887313755
intercept = -0.0954673765551862785
morphing_function = TF1('morphing_function','[0]*x+[1]',100.,500.)
morphing_function.SetParameter(0,slope)
morphing_function.SetParameter(1,intercept)

nEntries = tree.GetEntriesFast()
print '	# of entries: '+str(nEntries)
for entry in range(nEntries) :
	percent_done = 100.*entry/nEntries
	if percent_done%10 < 100./nEntries :
		print '	'+str((int)(percent_done))+'%'
	tree.GetEntry(entry)
	eventweight = LUMINOSITY*weight[0]*sf_pileup[0]*sf_top_pT[0]
	morphing_reweight = morphing_function.Eval(hadt_M[0])*eventweight
	if hadt_tau32[0] < 0.55 : #passing
		x_pass.Fill(cstar[0],eventweight)
		y_pass.Fill(x_F[0],eventweight)
		z_pass.Fill(M[0],eventweight)
	elif hadt_tau32 > 0.55 : #failing
		x_fail.Fill(cstar[0],eventweight)
		y_fail.Fill(x_F[0],eventweight)
		z_fail.Fill(M[0],eventweight)
		x_morph.Fill(cstar[0],morphing_reweight)
		y_morph.Fill(x_F[0],morphing_reweight)
		z_morph.Fill(M[0],morphing_reweight)		

x_pass.SetMaximum(1.05*max(x_pass.GetMaximum()+sqrt(x_pass.GetMaximum()),x_morph.GetMaximum()))
y_pass.SetMaximum(1.05*max(y_pass.GetMaximum()+sqrt(y_pass.GetMaximum()),y_morph.GetMaximum()))
z_pass.SetMaximum(1.05*max(z_pass.GetMaximum()+sqrt(z_pass.GetMaximum()),z_morph.GetMaximum()))

x_canv.cd()
#x_fail.Draw()
x_pass.Draw('PE1')
x_morph.Draw('SAME')
leg.Draw()
y_canv.cd()
#y_fail.Draw()
y_pass.Draw('PE1')
y_morph.Draw('SAME')
leg.Draw()
z_canv.cd()
#z_fail.Draw()
z_pass.Draw('PE1')
z_morph.Draw('SAME')
leg.Draw()

outfile = TFile('NTMJ_check_plots_morphed.root','recreate')
x_canv.Write()
y_canv.Write()
z_canv.Write()
outfile.Close()