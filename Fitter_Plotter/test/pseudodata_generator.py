import os
import sys
import glob
from array import array
from ROOT import *
from optparse import OptionParser

#Parse command line options
parser = OptionParser()
parser.add_option('--name', metavar='F', type='string', action='store',dest='name',default='pseudodata',help='')
parser.add_option('--nEvents', metavar='F', type='int', action='store',dest='nEvents',default='225000',help='')
#HARDER JET CUT DEFAULTS BELOW
parser.add_option('--Afb', metavar='F', type='float', action='store',dest='Afb',default='-0.0249',help='')
parser.add_option('--Rbck_4jet', metavar='F', type='float', action='store',dest='Rbck_4jet',default='0.2585',help='')
parser.add_option('--qqbar_frac_4jet', metavar='F', type='float', action='store',dest='qqbar_frac_4jet',default='0.13173',help='')
parser.add_option('--Rbck_5jet', metavar='F', type='float', action='store',dest='Rbck_5jet',default='0.1856',help='')
parser.add_option('--qqbar_frac_5jet', metavar='F', type='float', action='store',dest='qqbar_frac_5jet',default='0.08316392',help='')
parser.add_option('--xi_4jet', metavar='F', type='float', action='store',dest='xi_4jet',default='0.00',help='')
parser.add_option('--xi_5jet', metavar='F', type='float', action='store',dest='xi_5jet',default='0.00',help='')
parser.add_option('--delta_4jet', metavar='F', type='float', action='store',dest='delta_4jet',default='-0.256',help='')
parser.add_option('--delta_5jet', metavar='F', type='float', action='store',dest='delta_5jet',default='0.143',help='')
#SOFTER JET CUT DEFAULTS BELOW
#parser.add_option('--Afb', metavar='F', type='float', action='store',dest='Afb',default='-0.0210',help='')
#parser.add_option('--Rbck_4jet', metavar='F', type='float', action='store',dest='Rbck_4jet',default='0.3316',help='')
#parser.add_option('--qqbar_frac_4jet', metavar='F', type='float', action='store',dest='qqbar_frac_4jet',default='0.2276',help='')
#parser.add_option('--Rbck_5jet', metavar='F', type='float', action='store',dest='Rbck_5jet',default='0.2155',help='')
#parser.add_option('--qqbar_frac_5jet', metavar='F', type='float', action='store',dest='qqbar_frac_5jet',default='0.1405',help='')
#parser.add_option('--xi_4jet', metavar='F', type='float', action='store',dest='xi_4jet',default='0.00',help='')
#parser.add_option('--xi_5jet', metavar='F', type='float', action='store',dest='xi_5jet',default='0.00',help='')
#parser.add_option('--delta_4jet', metavar='F', type='float', action='store',dest='delta_4jet',default='0.228',help='')
#parser.add_option('--delta_5jet', metavar='F', type='float', action='store',dest='delta_5jet',default='-0.010',help='')
(options, args) = parser.parse_args()
if options.name == 'semilep_all' :
	print "Pick a different name for the output file. I'mma crash now."
	sys.exit('You made a mistake!')

#set constants
x = RooRealVar('x','x',-1.0,1.0)
y = RooRealVar('y','y',0.0,0.6)
z = RooRealVar('z','z',350.,1750.)
l = RooRealVar('l','l',-50.,100.)
#filename for total templates
filename  = 'aggregated_template_distributions_Powheg_no_rescale.root'

#load up file
file_pointer  = TFile(filename)
#load the template histograms
fqqs_plus_4jet = file_pointer.Get("f_qqs_plus_4jet")
fqqs_xi_plus_4jet = file_pointer.Get("f_qqs_xi_plus_4jet")
fqqs_delta_plus_4jet = file_pointer.Get("f_qqs_delta_plus_4jet")
fqqa_plus_4jet = file_pointer.Get("f_qqa_plus_4jet")
fqqa_xi_plus_4jet = file_pointer.Get("f_qqa_xi_plus_4jet")
fqqa_delta_plus_4jet = file_pointer.Get("f_qqa_delta_plus_4jet")
fqqs_minus_4jet = file_pointer.Get("f_qqs_minus_4jet")
fqqs_xi_minus_4jet = file_pointer.Get("f_qqs_xi_minus_4jet")
fqqs_delta_minus_4jet = file_pointer.Get("f_qqs_delta_minus_4jet")
fqqa_minus_4jet = file_pointer.Get("f_qqa_minus_4jet")
fqqa_xi_minus_4jet = file_pointer.Get("f_qqa_xi_minus_4jet")
fqqa_delta_minus_4jet = file_pointer.Get("f_qqa_delta_minus_4jet")
fqqs_plus_5jet = file_pointer.Get("f_qqs_plus_5jet")
fqqs_xi_plus_5jet = file_pointer.Get("f_qqs_xi_plus_5jet")
fqqs_delta_plus_5jet = file_pointer.Get("f_qqs_delta_plus_5jet")
fqqa_plus_5jet = file_pointer.Get("f_qqa_plus_5jet")
fqqa_xi_plus_5jet = file_pointer.Get("f_qqa_xi_plus_5jet")
fqqa_delta_plus_5jet = file_pointer.Get("f_qqa_delta_plus_5jet")
fqqs_minus_5jet = file_pointer.Get("f_qqs_minus_5jet")
fqqs_xi_minus_5jet = file_pointer.Get("f_qqs_xi_minus_5jet")
fqqs_delta_minus_5jet = file_pointer.Get("f_qqs_delta_minus_5jet")
fqqa_minus_5jet = file_pointer.Get("f_qqa_minus_5jet")
fqqa_xi_minus_5jet = file_pointer.Get("f_qqa_xi_minus_5jet")
fqqa_delta_minus_5jet = file_pointer.Get("f_qqa_delta_minus_5jet")
fgg_plus_4jet = file_pointer.Get("f_gg_plus_4jet")
fbck_plus_4jet = file_pointer.Get("f_bck_plus_4jet")
fgg_minus_4jet = file_pointer.Get("f_gg_minus_4jet")
fbck_minus_4jet = file_pointer.Get("f_bck_minus_4jet")
fgg_plus_5jet = file_pointer.Get("f_gg_plus_5jet")
fbck_plus_5jet = file_pointer.Get("f_bck_plus_5jet")
fgg_minus_5jet = file_pointer.Get("f_gg_minus_5jet")
fbck_minus_5jet = file_pointer.Get("f_bck_minus_5jet")
gtt_4jet = file_pointer.Get("gtt_4jet_local")
gbk_4jet = file_pointer.Get("gbk_4jet_local")
gtt_5jet = file_pointer.Get("gtt_5jet_local")
gbk_5jet = file_pointer.Get("gbk_5jet_local")
#renormalize
qq_rescale = 1.0/(fqqs_plus_4jet.Integral()+fqqs_minus_4jet.Integral()+fqqs_plus_5jet.Integral()+fqqs_minus_5jet.Integral())
fqqs_plus_4jet.Scale(qq_rescale)
fqqs_xi_plus_4jet.Scale(qq_rescale)
fqqs_delta_plus_4jet.Scale(qq_rescale)
fqqs_minus_4jet.Scale(qq_rescale)
fqqs_xi_minus_4jet.Scale(qq_rescale)
fqqs_delta_minus_4jet.Scale(qq_rescale)
fqqa_plus_4jet.Scale(qq_rescale)
fqqa_xi_plus_4jet.Scale(qq_rescale)
fqqa_delta_plus_4jet.Scale(qq_rescale)
fqqa_minus_4jet.Scale(qq_rescale)
fqqa_xi_minus_4jet.Scale(qq_rescale)
fqqa_delta_minus_4jet.Scale(qq_rescale)
fqqs_plus_5jet.Scale(qq_rescale)
fqqs_xi_plus_5jet.Scale(qq_rescale)
fqqs_delta_plus_5jet.Scale(qq_rescale)
fqqs_minus_5jet.Scale(qq_rescale)
fqqs_xi_minus_5jet.Scale(qq_rescale)
fqqs_delta_minus_5jet.Scale(qq_rescale)
fqqa_plus_5jet.Scale(qq_rescale)
fqqa_xi_plus_5jet.Scale(qq_rescale)
fqqa_delta_plus_5jet.Scale(qq_rescale)
fqqa_minus_5jet.Scale(qq_rescale)
fqqa_xi_minus_5jet.Scale(qq_rescale)
fqqa_delta_minus_5jet.Scale(qq_rescale)
gg_rescale = 1.0/(fgg_plus_4jet.Integral()+fgg_minus_4jet.Integral()+fgg_plus_5jet.Integral()+fgg_minus_5jet.Integral())
fgg_plus_4jet.Scale(gg_rescale)
fgg_minus_4jet.Scale(gg_rescale)
fgg_plus_5jet.Scale(gg_rescale)
fgg_minus_5jet.Scale(gg_rescale)
bck_rescale = 1.0/(fbck_plus_4jet.Integral()+fbck_minus_4jet.Integral()+fbck_plus_5jet.Integral()+fbck_minus_5jet.Integral())
fbck_plus_4jet.Scale(bck_rescale)
fbck_minus_4jet.Scale(bck_rescale)
fbck_plus_5jet.Scale(bck_rescale)
fbck_minus_5jet.Scale(bck_rescale)
fbck_plus_5jet.Scale(bck_rescale)
fbck_minus_5jet.Scale(bck_rescale)
gtt_rescale = 1.0/(gtt_4jet.Integral() + gtt_5jet.Integral())
gbk_rescale = 1.0/(gbk_4jet.Integral() + gbk_5jet.Integral())
gtt_4jet.Scale(gtt_rescale)
gbk_4jet.Scale(gbk_rescale)
gtt_5jet.Scale(gtt_rescale)
gbk_5jet.Scale(gbk_rescale)
F_xi_4jet = fqqs_xi_plus_4jet.Integral()+fqqs_xi_minus_4jet.Integral()
F_delta_4jet = fqqs_delta_plus_4jet.Integral()+fqqs_delta_minus_4jet.Integral()
F_xi_5jet = fqqs_xi_plus_5jet.Integral()+fqqs_xi_minus_5jet.Integral()
F_delta_5jet = fqqs_delta_plus_5jet.Integral()+fqqs_delta_minus_5jet.Integral()
#Add the symmetric and asymmetric histograms with xi and delta built in to get the total qqbar distributions
qqs_plus_4jet  = fqqs_plus_4jet.Clone('qqs_plus_4jet')
qqs_plus_4jet.Add(fqqs_xi_plus_4jet,options.xi_4jet)
qqs_plus_4jet.Add(fqqs_delta_plus_4jet,options.delta_4jet)
qqa_plus_4jet  = fqqa_plus_4jet.Clone('qqa_plus_4jet')
qqa_plus_4jet.Add(fqqa_xi_plus_4jet,options.xi_4jet)
qqa_plus_4jet.Add(fqqa_delta_plus_4jet,options.delta_4jet)
qqs_minus_4jet  = fqqs_minus_4jet.Clone('qqs_minus_4jet')
qqs_minus_4jet.Add(fqqs_xi_minus_4jet,options.xi_4jet)
qqs_minus_4jet.Add(fqqs_delta_minus_4jet,options.delta_4jet)
qqa_minus_4jet  = fqqa_minus_4jet.Clone('qqa_minus_4jet')
qqa_minus_4jet.Add(fqqa_xi_minus_4jet,options.xi_4jet)
qqa_minus_4jet.Add(fqqa_delta_minus_4jet,options.delta_4jet)
qqs_plus_5jet  = fqqs_plus_5jet.Clone('qqs_plus_5jet')
qqs_plus_5jet.Add(fqqs_xi_plus_5jet,options.xi_5jet)
qqs_plus_5jet.Add(fqqs_delta_plus_5jet,options.delta_5jet)
qqa_plus_5jet  = fqqa_plus_5jet.Clone('qqa_plus_5jet')
qqa_plus_5jet.Add(fqqa_xi_plus_5jet,options.xi_5jet)
qqa_plus_5jet.Add(fqqa_delta_plus_5jet,options.delta_5jet)
qqs_minus_5jet  = fqqs_minus_5jet.Clone('qqs_minus_5jet')
qqs_minus_5jet.Add(fqqs_xi_minus_5jet,options.xi_5jet)
qqs_minus_5jet.Add(fqqs_delta_minus_5jet,options.delta_5jet)
qqa_minus_5jet  = fqqa_minus_5jet.Clone('qqa_minus_5jet')
qqa_minus_5jet.Add(fqqa_xi_minus_5jet,options.xi_5jet)
qqa_minus_5jet.Add(fqqa_delta_minus_5jet,options.delta_5jet)
qq_plus_4jet = qqs_plus_4jet.Clone('qq_plus_4jet')
qq_plus_4jet.Add(qqa_plus_4jet,options.Afb)
qq_minus_4jet = qqs_minus_4jet.Clone('qq_minus_4jet')
qq_minus_4jet.Add(qqa_minus_4jet,options.Afb)
qq_plus_5jet = qqs_plus_5jet.Clone('qq_plus_5jet')
qq_plus_5jet.Add(qqa_plus_5jet,options.Afb)
qq_minus_5jet = qqs_minus_5jet.Clone('qq_minus_5jet')
qq_minus_5jet.Add(qqa_minus_5jet,options.Afb)
new_qq_rescale = 1.0/(qq_plus_4jet.Integral()+qq_minus_4jet.Integral()+qq_plus_5jet.Integral()+qq_minus_5jet.Integral())
qq_plus_4jet.Scale(options.qqbar_frac_4jet*new_qq_rescale)
qq_minus_4jet.Scale(options.qqbar_frac_4jet*new_qq_rescale)
qq_plus_5jet.Scale(options.qqbar_frac_5jet*new_qq_rescale)
qq_minus_5jet.Scale(options.qqbar_frac_5jet*new_qq_rescale)
sig_plus_4jet = qq_plus_4jet.Clone('sig_plus_4jet')
sig_plus_4jet.Add(fgg_plus_4jet,(1.0-options.qqbar_frac_4jet))
sig_minus_4jet = qq_minus_4jet.Clone('sig_minus_4jet')
sig_minus_4jet.Add(fgg_minus_4jet,(1.0-options.qqbar_frac_4jet))
sig_plus_5jet = qq_plus_5jet.Clone('sig_plus_5jet')
sig_plus_5jet.Add(fgg_plus_5jet,(1.0-options.qqbar_frac_5jet))
sig_minus_5jet = qq_minus_5jet.Clone('sig_minus_5jet')
sig_minus_5jet.Add(fgg_minus_5jet,(1.0-options.qqbar_frac_5jet))
sig_rescale = 1.0/(sig_plus_4jet.Integral()+sig_minus_4jet.Integral()+sig_plus_5jet.Integral()+sig_minus_5jet.Integral())
sig_plus_4jet.Scale(sig_rescale)
sig_minus_4jet.Scale(sig_rescale)
sig_plus_5jet.Scale(sig_rescale)
sig_minus_5jet.Scale(sig_rescale)
#Define the number of events to generate for each distribution
tt4 = gtt_4jet.Integral()
bk4 = gbk_4jet.Integral()
Rb4 = options.Rbck_4jet
Rb5 = options.Rbck_5jet
ns4 = 1./(1.+((1.-tt4)/tt4)*(1.+(Rb5/(1.-Rb5))*(1.+(bk4/(1.-bk4)))))
ns5 = 1./(1.+(tt4/(1.-tt4))*(1.+(Rb4/(1.-Rb4))*(1.+((1.-bk4)/bk4))))
nb4 = 1./(1.+((1.-bk4)/bk4)*(1.+((1.-Rb5)/Rb5)*(1.+(tt4/(1.-tt4)))))
nb5 = 1./(1.+(bk4/(1.-bk4))*(1.+((1.-Rb4)/Rb4)*(1.+((1.-tt4)/tt4))))
nsig = ns4 + ns5
nbck = nb4 + nb5
print 'TOTAL NORMALIZATION IS:  '+str(ns4)+' + '+str(ns5)+' + '+str(nb4)+' + '+str(nb5)+' = '+str(ns4+ns5+nb4+nb5)+''
nsemilep_plus_4jet  = RooRealVar('nsemilep_plus_4jet','nsemilep_plus_4jet',nsig*sig_plus_4jet.Integral())
nsemilep_minus_4jet = RooRealVar('nsemilep_minus_4jet','nsemilep_minus_4jet',nsig*sig_minus_4jet.Integral())
nbck_plus_4jet  	= RooRealVar('nbck_plus_4jet','nbck_plus_4jet',nbck*fbck_plus_4jet.Integral())
nbck_minus_4jet 	= RooRealVar('nbck_minus_4jet','nbck_minus_4jet',nbck*fbck_minus_4jet.Integral())
nsemilep_plus_5jet  = RooRealVar('nsemilep_plus_5jet','nsemilep_plus_5jet',nsig*sig_plus_5jet.Integral())
nsemilep_minus_5jet = RooRealVar('nsemilep_minus_5jet','nsemilep_minus_5jet',nsig*sig_minus_5jet.Integral())
nbck_plus_5jet  	= RooRealVar('nbck_plus_5jet','nbck_plus_5jet',nbck*fbck_plus_5jet.Integral())
nbck_minus_5jet 	= RooRealVar('nbck_minus_5jet','nbck_minus_5jet',nbck*fbck_minus_5jet.Integral())
ntt_4jet = RooRealVar('ntt_4jet','ntt_4jet',ns4)
nbk_4jet = RooRealVar('nbk_4jet','nbk_4jet',nb4)
ntt_5jet = RooRealVar('ntt_5jet','ntt_5jet',ns5)
nbk_5jet = RooRealVar('nbk_5jet','nbk_5jet',nb5)
#build RooDataHists
semilep_plus_4jet_rdh = RooDataHist('semilep_plus_4jet_rdh','semilep_plus_4jet_rdh',RooArgList(x,y,z),sig_plus_4jet)
semilep_minus_4jet_rdh = RooDataHist('semilep_minus_4jet_rdh','semilep_minus_4jet_rdh',RooArgList(x,y,z),sig_minus_4jet)
fbckplus_4jet_rdh = RooDataHist('fbckplus_4jet_rdh','fbckplus_4jet_rdh',RooArgList(x,y,z),fbck_plus_4jet)
fbckminus_4jet_rdh = RooDataHist('fbckminus_4jet_rdh','fbckminus_4jet_rdh',RooArgList(x,y,z),fbck_minus_4jet)
semilep_plus_5jet_rdh = RooDataHist('semilep_plus_5jet_rdh','semilep_plus_5jet_rdh',RooArgList(x,y,z),sig_plus_5jet)
semilep_minus_5jet_rdh = RooDataHist('semilep_minus_5jet_rdh','semilep_minus_5jet_rdh',RooArgList(x,y,z),sig_minus_5jet)
fbckplus_5jet_rdh = RooDataHist('fbckplus_5jet_rdh','fbckplus_5jet_rdh',RooArgList(x,y,z),fbck_plus_5jet)
fbckminus_5jet_rdh = RooDataHist('fbckminus_5jet_rdh','fbckminus_5jet_rdh',RooArgList(x,y,z),fbck_minus_5jet)
gtt_4jet_rdh = RooDataHist('gtt_4jet_rdh','gtt_4jet_rdh',RooArgList(l),gtt_4jet)
gbk_4jet_rdh = RooDataHist('gbk_4jet_rdh','gbk_4jet_rdh',RooArgList(l),gbk_4jet)
gtt_5jet_rdh = RooDataHist('gtt_5jet_rdh','gtt_5jet_rdh',RooArgList(l),gtt_5jet)
gbk_5jet_rdh = RooDataHist('gbk_5jet_rdh','gbk_5jet_rdh',RooArgList(l),gbk_5jet)
#build RooHistPDFs
#build RooHistPDFs
semilep_plus_4jet_pdf = RooHistPdf('semilep_plus_4jet_pdf','semilep_plus_4jet_pdf',RooArgSet(x,y,z),semilep_plus_4jet_rdh)
semilep_minus_4jet_pdf = RooHistPdf('semilep_minus_4jet_pdf','semilep_minus_4jet_pdf',RooArgSet(x,y,z),semilep_minus_4jet_rdh)
fbckplus_4jet_pdf = RooHistPdf('fbckplus_4jet_pdf','fbckplus_4jet_pdf',RooArgSet(x,y,z),fbckplus_4jet_rdh)
fbckminus_4jet_pdf = RooHistPdf('fbckminus_4jet_pdf','fbckminus_4jet_pdf',RooArgSet(x,y,z),fbckminus_4jet_rdh)
semilep_plus_5jet_pdf = RooHistPdf('semilep_plus_5jet_pdf','semilep_plus_5jet_pdf',RooArgSet(x,y,z),semilep_plus_5jet_rdh)
semilep_minus_5jet_pdf = RooHistPdf('semilep_minus_5jet_pdf','semilep_minus_5jet_pdf',RooArgSet(x,y,z),semilep_minus_5jet_rdh)
fbckplus_5jet_pdf = RooHistPdf('fbckplus_5jet_pdf','fbckplus_5jet_pdf',RooArgSet(x,y,z),fbckplus_5jet_rdh)
fbckminus_5jet_pdf = RooHistPdf('fbckminus_5jet_pdf','fbckminus_5jet_pdf',RooArgSet(x,y,z),fbckminus_5jet_rdh)
gtt_4jet_pdf = RooHistPdf('gtt_4jet_pdf','gtt_4jet_pdf',RooArgSet(l),gtt_4jet_rdh)
gbk_4jet_pdf = RooHistPdf('gbk_4jet_pdf','gbk_4jet_pdf',RooArgSet(l),gbk_4jet_rdh)
gtt_5jet_pdf = RooHistPdf('gtt_5jet_pdf','gtt_5jet_pdf',RooArgSet(l),gtt_5jet_rdh)
gbk_5jet_pdf = RooHistPdf('gbk_5jet_pdf','gbk_5jet_pdf',RooArgSet(l),gbk_5jet_rdh)
#Extend PDFs
semilep_plus_4jet_ext = RooExtendPdf('semilep_plus_4jet_ext','Background distribution, positive leptons',semilep_plus_4jet_pdf,nsemilep_plus_4jet)
semilep_minus_4jet_ext = RooExtendPdf('semilep_minus_4jet_ext','Background distribution, negative leptons',semilep_minus_4jet_pdf,nsemilep_minus_4jet)
fbckplus_4jet_ext = RooExtendPdf('fbckplus_4jet_ext','Background distribution, positive leptons',fbckplus_4jet_pdf,nbck_plus_4jet)
fbckminus_4jet_ext = RooExtendPdf('fbckminus_4jet_ext','Background distribution, negative leptons',fbckminus_4jet_pdf,nbck_minus_4jet)
semilep_plus_5jet_ext = RooExtendPdf('semilep_plus_5jet_ext','Background distribution, positive leptons',semilep_plus_5jet_pdf,nsemilep_plus_5jet)
semilep_minus_5jet_ext = RooExtendPdf('semilep_minus_5jet_ext','Background distribution, negative leptons',semilep_minus_5jet_pdf,nsemilep_minus_5jet)
fbckplus_5jet_ext = RooExtendPdf('fbckplus_5jet_ext','Background distribution, positive leptons',fbckplus_5jet_pdf,nbck_plus_5jet)
fbckminus_5jet_ext = RooExtendPdf('fbckminus_5jet_ext','Background distribution, negative leptons',fbckminus_5jet_pdf,nbck_minus_5jet)
gtt_4jet_ext = RooExtendPdf('gtt_4jet_ext','Likelihood distribution, signal',gtt_4jet_pdf,ntt_4jet)
gbk_4jet_ext = RooExtendPdf('gbk_4jet_ext','Likelihood distribution, signal',gbk_4jet_pdf,nbk_4jet)
gtt_5jet_ext = RooExtendPdf('gtt_5jet_ext','Likelihood distribution, signal',gtt_5jet_pdf,ntt_5jet)
gbk_5jet_ext = RooExtendPdf('gbk_5jet_ext','Likelihood distribution, signal',gbk_5jet_pdf,nbk_5jet)
#make total PDFs
total_PDF_plus_4jet_signal  = RooAddPdf('total_PDF_plus_4jet_signal','total_PDF_plus_4jet_signal',RooArgList(semilep_plus_4jet_ext))
total_PDF_minus_4jet_signal = RooAddPdf('total_PDF_minus_4jet_signal','total_PDF_minus_4jet_signal',RooArgList(semilep_minus_4jet_ext))
total_PDF_plus_5jet_signal  = RooAddPdf('total_PDF_plus_5jet_signal','total_PDF_plus_5jet_signal',RooArgList(semilep_plus_5jet_ext))
total_PDF_minus_5jet_signal = RooAddPdf('total_PDF_minus_5jet_signal','total_PDF_minus_5jet_signal',RooArgList(semilep_minus_5jet_ext))
#generate pseudodata events and save in an ASCII file
print '#signal, positive leptons, 4jets = '+str(int(options.nEvents*nsig*sig_plus_4jet.Integral()))
print '#signal, negative leptons, 4jets = '+str(int(options.nEvents*nsig*sig_minus_4jet.Integral()))
print '#signal, positive leptons, 5jets = '+str(int(options.nEvents*nsig*sig_plus_5jet.Integral()))
print '#signal, negative leptons, 5jets = '+str(int(options.nEvents*nsig*sig_minus_5jet.Integral()))
RooRandom.randomGenerator().SetSeed(0)
new_signal_plus_4jet  = total_PDF_plus_4jet_signal.generate(RooArgSet(x,y,z),int(options.nEvents*nsig*sig_plus_4jet.Integral()))
new_signal_minus_4jet = total_PDF_minus_4jet_signal.generate(RooArgSet(x,y,z),int(options.nEvents*nsig*sig_minus_4jet.Integral()))
new_signal_plus_5jet  = total_PDF_plus_5jet_signal.generate(RooArgSet(x,y,z),int(options.nEvents*nsig*sig_plus_5jet.Integral()))
new_signal_minus_5jet = total_PDF_minus_5jet_signal.generate(RooArgSet(x,y,z),int(options.nEvents*nsig*sig_minus_5jet.Integral()))
print '#background, positive leptons, 4jets = '+str(int(options.nEvents*nbck*fbck_plus_4jet.Integral()))
print '#background, negative leptons, 4jets = '+str(int(options.nEvents*nbck*fbck_minus_4jet.Integral()))
print '#background, positive leptons, 5jets = '+str(int(options.nEvents*nbck*fbck_plus_5jet.Integral()))
print '#background, negative leptons, 5jets = '+str(int(options.nEvents*nbck*fbck_minus_5jet.Integral()))
RooRandom.randomGenerator().SetSeed(0)
new_background_plus_4jet  = fbckplus_4jet_ext.generate(RooArgSet(x,y,z),int(options.nEvents*nbck*fbck_plus_4jet.Integral()))
new_background_minus_4jet = fbckminus_4jet_ext.generate(RooArgSet(x,y,z),int(options.nEvents*nbck*fbck_minus_4jet.Integral()))
new_background_plus_5jet  = fbckplus_5jet_ext.generate(RooArgSet(x,y,z),int(options.nEvents*nbck*fbck_plus_5jet.Integral()))
new_background_minus_5jet = fbckminus_5jet_ext.generate(RooArgSet(x,y,z),int(options.nEvents*nbck*fbck_minus_5jet.Integral()))
print '# signal events in 4 jets likelihood file = '+str(int(options.nEvents*nsig*(sig_plus_4jet.Integral()+sig_minus_4jet.Integral())))
print '# signal events in 5 jets likelihood file = '+str(int(options.nEvents*nsig*(sig_plus_5jet.Integral()+sig_minus_5jet.Integral())))
print '# background events in 4 jets likelihood file = '+str(int(options.nEvents*nbck*(fbck_plus_4jet.Integral()+fbck_minus_4jet.Integral())))
print '# background events in 5 jets likelihood file = '+str(int(options.nEvents*nbck*(fbck_plus_5jet.Integral()+fbck_minus_5jet.Integral())))
RooRandom.randomGenerator().SetSeed(0) 
new_signal_4jet_likelihood = gtt_4jet_ext.generate(RooArgSet(l),int(options.nEvents*nsig*(sig_plus_4jet.Integral()+sig_minus_4jet.Integral())))
new_signal_5jet_likelihood = gtt_5jet_ext.generate(RooArgSet(l),int(options.nEvents*nsig*(sig_plus_5jet.Integral()+sig_minus_5jet.Integral())))
new_background_4jet_likelihood = gbk_4jet_ext.generate(RooArgSet(l),int(options.nEvents*nbck*(fbck_plus_4jet.Integral()+fbck_minus_4jet.Integral())))
new_background_5jet_likelihood = gbk_5jet_ext.generate(RooArgSet(l),int(options.nEvents*nbck*(fbck_plus_5jet.Integral()+fbck_minus_5jet.Integral())))
new_signal_plus_4jet.write(options.name+'signal_plus_4jet.txt')
new_signal_minus_4jet.write(options.name+'signal_minus_4jet.txt')
new_signal_plus_5jet.write(options.name+'signal_plus_5jet.txt')
new_signal_minus_5jet.write(options.name+'signal_minus_5jet.txt')
#new_background_plus_4jet.write(options.name+'background_plus_4jet.txt')
#new_background_minus_4jet.write(options.name+'background_minus_4jet.txt')
#new_background_plus_5jet.write(options.name+'background_plus_5jet.txt')
#new_background_minus_5jet.write(options.name+'background_minus_5jet.txt')
new_signal_4jet_likelihood.write(options.name+'signal_4jet_likelihood.txt')
#new_background_4jet_likelihood.write(options.name+'background_4jet_likelihood.txt')
new_signal_5jet_likelihood.write(options.name+'signal_5jet_likelihood.txt')
#new_background_5jet_likelihood.write(options.name+'background_5jet_likelihood.txt')
#build list of likelihood values because I am stupid and I don't know how else to do this well
lines = [line.strip() for line in open(options.name+'signal_4jet_likelihood.txt','r')]
signal_4jet_likelihood_list = []
for line in lines :
	values = line.split()
	signal_4jet_likelihood_list.append(float(values[0]))
lines = [line.strip() for line in open(options.name+'signal_5jet_likelihood.txt','r')]
signal_5jet_likelihood_list = []
for line in lines :
	values = line.split()
	signal_5jet_likelihood_list.append(float(values[0]))
#lines = [line.strip() for line in open(options.name+'background_4jet_likelihood.txt','r')]
#background_4jet_likelihood_list =[]
#for line in lines :
#	values = line.split()
#	background_4jet_likelihood_list.append(float(values[0]))
#lines = [line.strip() for line in open(options.name+'background_5jet_likelihood.txt','r')]
#background_5jet_likelihood_list =[]
#for line in lines :
#	values = line.split()
#	background_5jet_likelihood_list.append(float(values[0]))
#build Ttree to write out
filename = 'angles_data_' + options.name + '.root'
f = TFile(filename, 'Recreate' )
output = TTree('angles_data', 'angles_data')
output.SetDirectory(0)
ttbar_mass = array('d',[0.])
Qt = array('d',[0.])
cos_theta_cs = array('d',[0.])
Feynman_x = array('d',[0.])
w_a = array('d',[0.])
w_s_xi = array('d',[0.])
w_s_delta = array('d',[0.])
w_a_xi = array('d',[0.])
w_a_delta = array('d',[0.])
w_a_opp = array('d',[0.])
w_s_xi_opp = array('d',[0.])
w_s_delta_opp = array('d',[0.])
w_a_xi_opp = array('d',[0.])
w_a_delta_opp = array('d',[0.])
Q_l = array('i',[0])
cos_theta_mc = array('d',[0.])
Feynman_x_mc = array('d',[0.])
ttbar_mass_mc = array('d',[0.])
lnL = array('d',[0.])
n_valid_jets = array('i',[0])
n_bTags = array('i',[0])
raw_pts = array('d',4*[0.])
event_weight = array('f',[0.0])
output.Branch('ttbar_mass',ttbar_mass,'ttbar_mass/D')
output.Branch('Qt',Qt,'Qt/D')
output.Branch('cos_theta_cs',cos_theta_cs,'cos_theta_cs/D')
output.Branch('Feynman_x',Feynman_x,'Feynman_x/D')
output.Branch('Q_l',Q_l,'Q_l/I')
output.Branch('cos_theta_mc',cos_theta_mc,'cos_theta_mc/D')
output.Branch('Feynman_x_mc',Feynman_x_mc,'Feynman_x_mc/D')
output.Branch('ttbar_mass_mc',ttbar_mass_mc,'ttbar_mass_mc/D')
output.Branch('lnL',lnL,'lnL/D')
output.Branch('n_valid_jets',n_valid_jets,'n_valid_jets/I')
output.Branch('n_bTags',n_bTags,'n_bTags/I')
output.Branch('raw_pts',raw_pts,'raw_pts[4]/D')
output.Branch('event_weight',event_weight,'event_weight/f')
output.Branch('w_a',w_a,'w_a/D')
output.Branch('w_s_xi',w_s_xi,'w_s_xi/D')
output.Branch('w_s_delta',w_s_delta,'w_s_delta/D')
output.Branch('w_a_xi',w_a_xi,'w_a_xi/D')
output.Branch('w_a_delta',w_a_delta,'w_a_delta/D')
output.Branch('w_a_opp',w_a_opp,'w_a_opp/D')
output.Branch('w_s_xi_opp',w_s_xi_opp,'w_s_xi_opp/D')
output.Branch('w_s_delta_opp',w_s_delta_opp,'w_s_delta_opp/D')
output.Branch('w_a_xi_opp',w_a_xi_opp,'w_a_xi_opp/D')
output.Branch('w_a_delta_opp',w_a_delta_opp,'w_a_delta_opp/D')
#open ASCII files and write contents to new Ttree
signal_4jet_ascii_file_names = [options.name+'signal_plus_4jet.txt',options.name+'signal_minus_4jet.txt']
iterator = 0
for ascii_file in signal_4jet_ascii_file_names :
	lines = [line.strip() for line in open(ascii_file,'r')]
	for line in lines :
		values = line.split()
		costheta = float(values[0])
		xF = float(values[1])
		mttbar = float(values[2])
		cos_theta_cs[0] = costheta
		Feynman_x[0] = xF
		ttbar_mass[0] = mttbar
		lnL[0] = signal_4jet_likelihood_list[iterator]
		Qt[0] = 0.
		cos_theta_mc[0] = 0.
		Feynman_x_mc[0] = 0.
		ttbar_mass_mc[0] = 0.
		n_bTags[0] = 0
		w_a[0] = 0.
		w_s_xi[0] = 0.
		w_s_delta[0] = 0.
		w_a_xi[0] = 0.
		w_a_delta[0] = 0.
		w_a_opp[0] = 0.
		w_s_xi_opp[0] = 0.
		w_s_delta_opp[0] = 0.
		w_a_xi_opp[0] = 0.
		w_a_delta_opp[0] = 0.
		n_valid_jets[0] = 4
		event_weight[0] = 1.0
		if ascii_file == options.name+'signal_plus_4jet.txt' :
			Q_l[0] = 1
		if ascii_file == options.name+'signal_minus_4jet.txt' :
			Q_l[0] = -1
		for i in range(4) :
			raw_pts[i] = 100.
		output.Fill()
		iterator =iterator+1
signal_5jet_ascii_file_names = [options.name+'signal_plus_5jet.txt',options.name+'signal_minus_5jet.txt']
iterator = 0
for ascii_file in signal_5jet_ascii_file_names :
	lines = [line.strip() for line in open(ascii_file,'r')]
	for line in lines :
		values = line.split()
		costheta = float(values[0])
		xF = float(values[1])
		mttbar = float(values[2])
		cos_theta_cs[0] = costheta
		Feynman_x[0] = xF
		ttbar_mass[0] = mttbar
		lnL[0] = signal_5jet_likelihood_list[iterator]
		Qt[0] = 0.
		w_a[0] = 0.
		w_s_xi[0] = 0.
		w_s_delta[0] = 0.
		w_a_xi[0] = 0.
		w_a_delta[0] = 0.
		w_a_opp[0] = 0.
		w_s_xi_opp[0] = 0.
		w_s_delta_opp[0] = 0.
		w_a_xi_opp[0] = 0.
		w_a_delta_opp[0] = 0.
		cos_theta_mc[0] = 0.
		Feynman_x_mc[0] = 0.
		ttbar_mass_mc[0] = 0.
		n_bTags[0] = 0
		n_valid_jets[0] = 5
		event_weight[0] = 1.0
		if ascii_file == options.name+'signal_plus_5jet.txt' :
			Q_l[0] = 1
		if ascii_file == options.name+'signal_minus_5jet.txt' :
			Q_l[0] = -1
		for i in range(4) :
			raw_pts[i] = 100.
		output.Fill()
		iterator =iterator+1
#background_4jet_ascii_file_names = [options.name+'background_plus_4jet.txt',options.name+'background_minus_4jet.txt']
#iterator = 0
#for ascii_file in background_4jet_ascii_file_names :
#	lines = [line.strip() for line in open(ascii_file,'r')]
#	for line in lines :
#		values = line.split()
#		costheta = float(values[0])
#		xF = float(values[1])
#		mttbar = float(values[2])
#		cos_theta_cs[0] = costheta
#		Feynman_x[0] = xF
#		ttbar_mass[0] = mttbar
#		lnL[0] = background_4jet_likelihood_list[iterator]
#		Qt[0] = 0.
#		w_a[0] = 0.
#		w_s_xi[0] = 0.
#		w_s_delta[0] = 0.
#		w_a_xi[0] = 0.
#		w_a_delta[0] = 0.
#		w_a_opp[0] = 0.
#		w_s_xi_opp[0] = 0.
#		w_s_delta_opp[0] = 0.
#		w_a_xi_opp[0] = 0.
#		w_a_delta_opp[0] = 0.
#		cos_theta_mc[0] = 0.
#		Feynman_x_mc[0] = 0.
#		ttbar_mass_mc[0] = 0.
#		n_bTags[0] = 0
#		n_valid_jets[0] = 4
#		event_weight[0] = 1.0
#		if ascii_file == options.name+'background_plus_4jet.txt' :
#			Q_l[0] = 1
#		if ascii_file == options.name+'background_minus_4jet.txt' :
#			Q_l[0] = -1
#		for i in range(4) :
#			raw_pts[i] = 100.
#		output.Fill()
#		iterator =iterator+1
#background_5jet_ascii_file_names = [options.name+'background_plus_5jet.txt',options.name+'background_minus_5jet.txt']
#iterator = 0
#for ascii_file in background_5jet_ascii_file_names :
#	lines = [line.strip() for line in open(ascii_file,'r')]
#	for line in lines :
#		values = line.split()
#		costheta = float(values[0])
#		xF = float(values[1])
#		mttbar = float(values[2])
#		cos_theta_cs[0] = costheta
#		Feynman_x[0] = xF
#		ttbar_mass[0] = mttbar
#		lnL[0] = background_5jet_likelihood_list[iterator]
#		Qt[0] = 0.
#		w_a[0] = 0.
#		w_s_xi[0] = 0.
#		w_s_delta[0] = 0.
#		w_a_xi[0] = 0.
#		w_a_delta[0] = 0.
#		w_a_opp[0] = 0.
#		w_s_xi_opp[0] = 0.
#		w_s_delta_opp[0] = 0.
#		w_a_xi_opp[0] = 0.
#		w_a_delta_opp[0] = 0.
#		cos_theta_mc[0] = 0.
#		Feynman_x_mc[0] = 0.
#		ttbar_mass_mc[0] = 0.
#		n_bTags[0] = 0
#		n_valid_jets[0] = 5
#		event_weight[0] = 1.0
#		if ascii_file == options.name+'background_plus_5jet.txt' :
#			Q_l[0] = 1
#		if ascii_file == options.name+'background_minus_5jet.txt' :
#			Q_l[0] = -1
#		for i in range(4) :
#			raw_pts[i] = 100.
#		output.Fill()
#		iterator =iterator+1

f.cd()
output.Write()
f.Close()
for name in signal_4jet_ascii_file_names :
	os.system('rm '+name)
#for name in background_4jet_ascii_file_names :
#	os.system('rm '+name)
os.system('rm '+options.name+'signal_4jet_likelihood.txt')
#os.system('rm '+options.name+'background_4jet_likelihood.txt')
for name in signal_5jet_ascii_file_names :
	os.system('rm '+name)
#for name in background_5jet_ascii_file_names :
#	os.system('rm '+name)
os.system('rm '+options.name+'signal_5jet_likelihood.txt')
#os.system('rm '+options.name+'background_5jet_likelihood.txt')
