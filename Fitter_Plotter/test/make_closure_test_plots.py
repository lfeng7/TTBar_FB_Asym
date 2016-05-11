from ROOT import *
from array import array

#settings
par_of_int = 'AFB'
#par_of_int = 'mu'
#par_of_int = 'd'
#data_fit_central_value = 0.105512937777 #AFB WITH SYS
data_fit_central_value = 0.114700063076 #AFB NO SYS
#data_fit_central_value = 0.0706306681281 #MU WITH SYS
#data_fit_central_value = 0.0831232619589 #MU NO SYS
#data_fit_central_value = 9.58443635568e-09 #D WITH SYS
#data_fit_central_value = 2.65831801016e-12 #D NO SYS
input_values = [-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,-0.025,0.0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.3,0.35,0.4,0.45,0.5]
if par_of_int!='AFB' :
	input_values = [-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
#append = '_everything'
#append = '_fix_mu_and_d'
#append = '_fix_AFB_and_d'
#append = '_fix_AFB_and_mu'
append = '_only_'+par_of_int

par_values_filename = 'par_Afb_values'+append+'.txt'
if par_of_int != 'AFB' :
	par_values_filename = 'par_'+par_of_int+'_values'+append+'.txt'
output_filename = 'par_'+par_of_int+'_closure_test_plots'+append+'.root'

#output file
outfile = TFile(output_filename,'recreate')

#TGraphErrors setup
n = len(input_values)
x_values = array('d',n*[1.0])
y_values = array('d',n*[1.0])
x_errs   = array('d',n*[0.1])
y_1_errs = array('d',n*[1.0])
y_2_errs = array('d',n*[1.0])
onesigma_graph = TGraphAsymmErrors(n,x_values,y_values,x_errs,x_errs,y_1_errs,y_1_errs)
twosigma_graph = TGraphAsymmErrors(n,x_values,y_values,x_errs,x_errs,y_2_errs,y_2_errs)

#closure test plot setup
binwidth = input_values[1]-input_values[0]
lowvalue = input_values[0]-binwidth/2.
hivalue  = input_values[n-1]+binwidth/2.
title = 'Fitting Closure Tests'
if append.find('_only_')!=-1 :
	title += ' (no systematics)'
if par_of_int == 'AFB' :
	title += '; Input A_{FB}; Fitted A_{FB}'
elif par_of_int == 'mu' :
	title += '; Input #mu; Fitted #mu'
elif par_of_int == 'd' :
	title += '; Input d; Fitted d'
closure_test_histo = TH2D('closure_test_histo',title,n,lowvalue,hivalue,4*n,2*lowvalue,2*hivalue)

#Set graph attributes
onesigma_graph.SetFillColor(kGreen)
onesigma_graph.SetLineStyle(2)
onesigma_graph.SetLineWidth(3)
onesigma_graph.SetLineColor(kBlack)
onesigma_graph.SetMarkerStyle(20)
onesigma_graph.SetMarkerColor(kBlack)
twosigma_graph.SetFillColor(kYellow)
if par_of_int=='AFB' :
	twosigma_graph.SetTitle('Neyman Construction for A_{FB}')
	twosigma_graph.GetXaxis().SetTitle('Input A_{FB}')
	twosigma_graph.GetYaxis().SetTitle('Fitted A_{FB}')
elif par_of_int=='mu' :
	twosigma_graph.SetTitle('Neyman Construction for #mu')
	twosigma_graph.GetXaxis().SetTitle('Input #mu')
	twosigma_graph.GetYaxis().SetTitle('Fitted #mu')
elif par_of_int=='d' :
	twosigma_graph.SetTitle('Neyman Construction for d')
	twosigma_graph.GetXaxis().SetTitle('Input d')
	twosigma_graph.GetYaxis().SetTitle('Fitted d')

#for each of the data points tested
for i in range(n) :
	#build the filepath
	filepath = './CLOSURE_TEST_'+par_of_int.upper()+'_'
	if input_values[i] < 0 :
		filepath+='-'
	filepath+=str(abs(input_values[i])).split('.')[0]+str(abs(input_values[i])).split('.')[1]
	filepath+='/'+par_values_filename
	print 'Getting '+par_of_int+' values from file at path '+filepath
	#open the afb values file
	par_values_file = open(filepath)
	#put the values in a list
	lines = par_values_file.readlines()
	#read the values into the closure test histogram
	realvalues = []
	for value in lines :
		realvalue = eval(value)#-0.1+input_values[i]
		realvalues.append(realvalue)
		closure_test_histo.Fill(input_values[i],realvalue,1./len(lines))
	print '	Got '+str(len(realvalues))+' values'
	#find the relevant values for the Neyman construction
	twosigma_low = realvalues[int(round(0.025*len(realvalues)))-1]
	onesigma_low = realvalues[int(round(0.160*len(realvalues)))-1]
	mean 		 = realvalues[int(round(0.500*len(realvalues)))-1]
	onesigma_hi  = realvalues[int(round(0.840*len(realvalues)))-1]
	twosigma_hi  = realvalues[int(round(0.975*len(realvalues)))-1]
	#Set graphs' point values
	onesigma_graph.SetPoint(i,input_values[i],mean)
	twosigma_graph.SetPoint(i,input_values[i],mean)
	onesigma_graph.SetPointError(i,x_errs[i],x_errs[i],abs(onesigma_low-mean),abs(onesigma_hi-mean))
	twosigma_graph.SetPointError(i,x_errs[i],x_errs[i],abs(twosigma_low-mean),abs(twosigma_hi-mean))

#Interpolate given the central value
data_fit_one_sigma_down = 2*lowvalue
data_fit_one_sigma_up   = 2*hivalue
data_fit_mean = data_fit_central_value
for i in range(n-1) :
	thisx = array('d',[0.0]); thisymean = array('d',[0.0])
	nextx = array('d',[0.0]); nextymean = array('d',[0.0])
	onesigma_graph.GetPoint(i,thisx,thisymean); onesigma_graph.GetPoint(i+1,nextx,nextymean)
	thisyhi  = thisymean[0]+onesigma_graph.GetErrorYhigh(i)
	thisylow = thisymean[0]-onesigma_graph.GetErrorYlow(i)
	nextyhi  = nextymean[0]+onesigma_graph.GetErrorYhigh(i+1)
	nextylow = nextymean[0]-onesigma_graph.GetErrorYlow(i+1)
	if thisyhi <= data_fit_central_value and nextyhi > data_fit_central_value :
		slope = (nextyhi-thisyhi)/(input_values[i+1]-input_values[i])
		data_fit_one_sigma_down = (data_fit_central_value-thisyhi)/slope+input_values[i]
	if thisylow <= data_fit_central_value and nextylow > data_fit_central_value :
		slope = (nextylow-thisylow)/(input_values[i+1]-input_values[i])
		data_fit_one_sigma_up = (data_fit_central_value-thisylow)/slope+input_values[i]
	if thisymean[0] <= data_fit_central_value and nextymean[0] > data_fit_central_value :
		slope = (nextymean[0]-thisymean[0])/(input_values[i+1]-input_values[i])
		data_fit_mean = (data_fit_central_value-thisymean[0])/slope+input_values[i]

print 'Central value = %.5f, plus one sigma = %.5f, minus one sigma = %.5f'%(data_fit_mean,data_fit_one_sigma_up,data_fit_one_sigma_down)
print 'FINAL RESULT: Parameter %s = %.5f + %.5f - %.5f'%(par_of_int,data_fit_mean,data_fit_one_sigma_up-data_fit_mean,abs(data_fit_one_sigma_down-data_fit_mean))

#Build the lines to indicate based on the data fit value
lines = []
x_low = twosigma_graph.GetXaxis().GetBinLowEdge(twosigma_graph.GetXaxis().GetFirst())
x_high = twosigma_graph.GetXaxis().GetBinUpEdge(twosigma_graph.GetXaxis().GetLast())
y_low = twosigma_graph.GetYaxis().GetBinLowEdge(twosigma_graph.GetYaxis().GetFirst())
lines.append(TLine(x_low,data_fit_central_value,data_fit_one_sigma_up,data_fit_central_value))
lines.append(TLine(data_fit_one_sigma_down,data_fit_central_value,data_fit_one_sigma_down,y_low))
lines.append(TLine(data_fit_one_sigma_up,data_fit_central_value,data_fit_one_sigma_up,y_low))
lines.append(TLine(data_fit_mean,data_fit_central_value,data_fit_mean,y_low))
for line in lines :
	line.SetLineWidth(3); line.SetLineColor(kRed); line.SetLineStyle(2)

#Build a legend
leg = TLegend(0.62,0.67,0.9,0.9)
leg.AddEntry(onesigma_graph,'Mean values','PL')
leg.AddEntry(onesigma_graph,'#pm 1 #sigma','F')
leg.AddEntry(twosigma_graph,'#pm 2 #sigma','F')
leg.AddEntry(lines[0],'data fit','L')

#The line to go on the closure test plots
otherline = TLine(lowvalue,lowvalue,hivalue,hivalue)
otherline.SetLineColor(kBlack)
otherline.SetLineWidth(4)

#Plot the neyman construction histogram
neyman_canv = TCanvas('neyman_canv','neyman_canv',900,900)
neyman_canv.cd()
twosigma_graph.Draw('A E3')
onesigma_graph.Draw('SAME E3')
onesigma_graph.Draw('SAME PLX')
for line in lines :
	line.Draw()
leg.Draw()

#plot the closure test histogram
closure_test_canv = TCanvas('closure_test_canv','closure_test_canv',900,900)
closure_test_canv.cd()
closure_test_histo.Draw('COL')
otherline.Draw()

neyman_canv.Write()
closure_test_canv.Write()
outfile.Close()
