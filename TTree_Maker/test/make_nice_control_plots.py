#make_nice_control_plots.py makes nice control plots based on cuts
#reads in a TTree File with a control_plots_folder and the plots

import ROOT
from optparse import OptionParser

#Global vars
MET_PT_MIN = 10. #GeV
MU_PT_MIN = 45 #GeV
MU_ETA_MAX = 2.1
EL_PT_MIN = 35 #GeV
EL_ETA_MAX = 2.5
OTHER_MU_PT_MIN = 35 #GeV
OTHER_MU_ETA_MAX = 2.1
OTHER_EL_PT_MIN = 35 #GeV
OTHER_EL_ETA_MAX = 2.5

parser = OptionParser()
#Run options
parser.add_option('--infile', type='string', action='store', default='file', dest='infile', help='Path to input TTree file')
parser.add_option('--outfile', type='string', action='store', default='file', dest='outfile', help='name of output file')
(options, args) = parser.parse_args()

#get the folder of plots
f = TFile(options.infile).Get('control_plots_folder')

#output file
fout = TFile(options.outfile,'recreate')
fout.cd()

#plots

#metpt 
met_pt = f.Get('met_pt')
line = ROOT.TLine(MET_PT_MIN,0.0,MET_PT_MIN,met_pt.GetMaximum())
arrow = ROOT.TArrow(MET_PT_MIN,met_pt.GetMaximum()/2.0,MET_PT_MIN+40.,met_pt.GetMaximum()/2.0)
met_pt.SetLineWidth(2); line.SetLineWidth(2); arrow.SetLineWidth(2)
line.SetLineColor(ROOT.kRed); arrow.SetLineColor(ROOT.kRed)
met_pt_canv = ROOT.TCanvas('met_pt_canv','met_pt_canv',900.,750.)
met_pt_canv.cd()
met_pt.Draw(); line.Draw('SAME'); arrow.Draw('SAME')
fout.cd(); met_pt_canv.Write()

#lepton 1 pT
lep1_pt = f.Get('lep1_pt')
line = ROOT.TLine(MU_PT_MIN,0.0,MU_PT_MIN,lep1_pt.GetMaximum())
arrow = ROOT.TArrow(MU_PT_MIN,lep1_pt.GetMaximum()/2.0,MU_PT_MIN+40.,lep1_pt.GetMaximum()/2.0)
lep1_pt.SetLineWidth(2); line.SetLineWidth(2); arrow.SetLineWidth(2)
line.SetLineColor(ROOT.kRed); arrow.SetLineColor(ROOT.kRed)
lep1_pt_canv = ROOT.TCanvas('lep1_pt_canv','lep1_pt_canv',900.,750.)
lep1_pt_canv.cd()
lep1_pt.Draw(); line.Draw('SAME'); arrow.Draw('SAME')
fout.cd(); lep1_pt_canv.Write()

#lepton 2 pT
lep2_pt = f.Get('lep2_pt')
line = ROOT.TLine(MU_PT_MIN,0.0,MU_PT_MIN,lep2_pt.GetMaximum())
arrow = ROOT.TArrow(MU_PT_MIN,lep2_pt.GetMaximum()/2.0,MU_PT_MIN-15.,lep2_pt.GetMaximum()/2.0,'<')
lep2_pt.SetLineWidth(2); line.SetLineWidth(2); arrow.SetLineWidth(2)
line.SetLineColor(ROOT.kRed); arrow.SetLineColor(ROOT.kRed)
lep2_pt_canv = ROOT.TCanvas('lep2_pt_canv','lep2_pt_canv',900.,750.)
lep2_pt_canv.cd()
lep2_pt.Draw(); line.Draw('SAME'); arrow.Draw('SAME')
fout.cd(); lep2_pt_canv.Write()

#other lepton 1 pT
other_lep1_pt = f.Get('other_lep1_pt')
line = ROOT.TLine(OTHER_EL_PT_MIN,0.0,OTHER_EL_PT_MIN,other_lep1_pt.GetMaximum())
arrow = ROOT.TArrow(OTHER_EL_PT_MIN,other_lep1_pt.GetMaximum()/2.0,OTHER_EL_PT_MIN-15.,other_lep1_pt.GetMaximum()/2.0,'<')
other_lep1_pt.SetLineWidth(2); line.SetLineWidth(2); arrow.SetLineWidth(2)
line.SetLineColor(ROOT.kRed); arrow.SetLineColor(ROOT.kRed)
other_lep1_pt_canv = ROOT.TCanvas('other_lep1_pt_canv','other_lep1_pt_canv',900.,750.)
other_lep1_pt_canv.cd()
other_lep1_pt.Draw(); line.Draw('SAME'); arrow.Draw('SAME')
fout.cd(); other_lep1_pt_canv.Write()

#other lepton 2 pT
other_lep2_pt = f.Get('other_lep2_pt')
line = ROOT.TLine(OTHER_EL_PT_MIN,0.0,OTHER_EL_PT_MIN,other_lep2_pt.GetMaximum())
arrow = ROOT.TArrow(OTHER_EL_PT_MIN,other_lep2_pt.GetMaximum()/2.0,OTHER_EL_PT_MIN-15.,other_lep2_pt.GetMaximum()/2.0,'<')
other_lep2_pt.SetLineWidth(2); line.SetLineWidth(2); arrow.SetLineWidth(2)
line.SetLineColor(ROOT.kRed); arrow.SetLineColor(ROOT.kRed)
other_lep2_pt_canv = ROOT.TCanvas('other_lep2_pt_canv','other_lep2_pt_canv',900.,750.)
other_lep2_pt_canv.cd()
other_lep2_pt.Draw(); line.Draw('SAME'); arrow.Draw('SAME')
fout.cd(); other_lep2_pt_canv.Write()

#lepton 1 eta
lep1_eta = f.Get('lep1_eta')
line1 = ROOT.TLine(MU_ETA_MAX,0.0,MU_ETA_MAX,lep1_eta.GetMaximum())
line2 = ROOT.TLine(-1.0*MU_ETA_MAX,0.0,-1.0*MU_ETA_MAX,lep1_eta.GetMaximum())
arrow1 = ROOT.TArrow(MU_ETA_MAX,lep1_eta.GetMaximum()/2.0,MU_ETA_MAX-0.3.,lep1_eta.GetMaximum()/2.0,'<')
arrow2 = ROOT.TArrow(-1.0*MU_ETA_MAX,lep1_eta.GetMaximum()/2.0,-1.0*MU_ETA_MAX+0.3.,lep1_eta.GetMaximum()/2.0,'>')
lep1_eta.SetLineWidth(2); line1.SetLineWidth(2); arrow1.SetLineWidth(2); line2.SetLineWidth(2); arrow2.SetLineWidth(2)
line1.SetLineColor(ROOT.kRed); arrow1.SetLineColor(ROOT.kRed); line2.SetLineColor(ROOT.kRed); arrow2.SetLineColor(ROOT.kRed)
lep1_eta_canv = ROOT.TCanvas('lep1_eta_canv','lep1_eta_canv',900.,750.)
lep1_eta_canv.cd()
lep1_eta.Draw(); line1.Draw('SAME'); arrow1.Draw('SAME'); line2.Draw('SAME'); arrow2.Draw('SAME')
fout.cd(); lep1_eta_canv.Write()

#lepton 2 eta
lep2_eta = f.Get('lep2_eta')
line1 = ROOT.TLine(MU_ETA_MAX,0.0,MU_ETA_MAX,lep2_eta.GetMaximum())
line2 = ROOT.TLine(-1.0*MU_ETA_MAX,0.0,-1.0*MU_ETA_MAX,lep2_eta.GetMaximum())
arrow1 = ROOT.TArrow(MU_ETA_MAX,lep2_eta.GetMaximum()/2.0,MU_ETA_MAX-0.3.,lep2_eta.GetMaximum()/2.0,'<')
arrow2 = ROOT.TArrow(-1.0*MU_ETA_MAX,lep2_eta.GetMaximum()/2.0,-1.0*MU_ETA_MAX+0.3.,lep2_eta.GetMaximum()/2.0,'>')
lep2_eta.SetLineWidth(2); line1.SetLineWidth(2); arrow1.SetLineWidth(2); line2.SetLineWidth(2); arrow2.SetLineWidth(2)
line1.SetLineColor(ROOT.kRed); arrow1.SetLineColor(ROOT.kRed); line2.SetLineColor(ROOT.kRed); arrow2.SetLineColor(ROOT.kRed)
lep2_eta_canv = ROOT.TCanvas('lep2_eta_canv','lep2_eta_canv',900.,750.)
lep2_eta_canv.cd()
lep2_eta.Draw(); line1.Draw('SAME'); arrow1.Draw('SAME'); line2.Draw('SAME'); arrow2.Draw('SAME')
fout.cd(); lep2_eta_canv.Write()

#other lepton 1 eta
other_lep1_eta = f.Get('other_lep1_eta')
line1 = ROOT.TLine(OTHER_EL_ETA_MAX,0.0,OTHER_EL_ETA_MAX,other_lep1_eta.GetMaximum())
line2 = ROOT.TLine(-1.0*OTHER_EL_ETA_MAX,0.0,-1.0*OTHER_EL_ETA_MAX,other_lep1_eta.GetMaximum())
arrow1 = ROOT.TArrow(OTHER_EL_ETA_MAX,other_lep1_eta.GetMaximum()/2.0,OTHER_EL_ETA_MAX+0.3.,other_lep1_eta.GetMaximum()/2.0,'>')
arrow2 = ROOT.TArrow(-1.0*OTHER_EL_ETA_MAX,other_lep1_eta.GetMaximum()/2.0,-1.0*OTHER_EL_ETA_MAX-0.3.,other_lep1_eta.GetMaximum()/2.0,'<')
other_lep1_eta.SetLineWidth(2); line1.SetLineWidth(2); arrow1.SetLineWidth(2); line2.SetLineWidth(2); arrow2.SetLineWidth(2)
line1.SetLineColor(ROOT.kRed); arrow1.SetLineColor(ROOT.kRed); line2.SetLineColor(ROOT.kRed); arrow2.SetLineColor(ROOT.kRed)
other_lep1_eta_canv = ROOT.TCanvas('other_lep1_eta_canv','other_lep1_eta_canv',900.,750.)
other_lep1_eta_canv.cd()
other_lep1_eta.Draw(); line1.Draw('SAME'); arrow1.Draw('SAME'); line2.Draw('SAME'); arrow2.Draw('SAME')
fout.cd(); other_lep1_eta_canv.Write()

#other lepton 2 eta
other_lep2_eta = f.Get('other_lep2_eta')
line1 = ROOT.TLine(OTHER_EL_ETA_MAX,0.0,OTHER_EL_ETA_MAX,other_lep2_eta.GetMaximum())
line2 = ROOT.TLine(-1.0*OTHER_EL_ETA_MAX,0.0,-1.0*OTHER_EL_ETA_MAX,other_lep2_eta.GetMaximum())
arrow1 = ROOT.TArrow(OTHER_EL_ETA_MAX,other_lep2_eta.GetMaximum()/2.0,OTHER_EL_ETA_MAX+0.3.,other_lep2_eta.GetMaximum()/2.0,'>')
arrow2 = ROOT.TArrow(-1.0*OTHER_EL_ETA_MAX,other_lep2_eta.GetMaximum()/2.0,-1.0*OTHER_EL_ETA_MAX-0.3.,other_lep2_eta.GetMaximum()/2.0,'<')
other_lep2_eta.SetLineWidth(2); line1.SetLineWidth(2); arrow1.SetLineWidth(2); line2.SetLineWidth(2); arrow2.SetLineWidth(2)
line1.SetLineColor(ROOT.kRed); arrow1.SetLineColor(ROOT.kRed); line2.SetLineColor(ROOT.kRed); arrow2.SetLineColor(ROOT.kRed)
other_lep2_eta_canv = ROOT.TCanvas('other_lep2_eta_canv','other_lep2_eta_canv',900.,750.)
other_lep2_eta_canv.cd()
other_lep2_eta.Draw(); line1.Draw('SAME'); arrow1.Draw('SAME'); line2.Draw('SAME'); arrow2.Draw('SAME')
fout.cd(); other_lep2_eta_canv.Write()