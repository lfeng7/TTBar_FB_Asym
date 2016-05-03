from ROOT import *
from array import array

#settings
input_values = [-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5]
filename = 'afb_values.txt'

#output file
outfile = TFile('closure_test_plot.root','recreate')

#TGraphErrors setup
n = len(input_values)
x_values = array('d',n*[1.0])
y_values = array('d',n*[1.0])
x_errs   = array('d',n*[0.1])
y_1_errs = array('d',n*[1.0])
y_2_errs = array('d',n*[1.0])
onesigma_graph = TGraphAsymmErrors(n,x_values,y_values,x_errs,x_errs,y_1_errs,y_1_errs)
twosigma_graph = TGraphAsymmErrors(n,x_values,y_values,x_errs,x_errs,y_2_errs,y_2_errs)

#Set graph attributes
onesigma_graph.SetFillColor(kGreen)
onesigma_graph.SetLineStyle(2)
onesigma_graph.SetLineWidth(3)
onesigma_graph.SetLineColor(kBlack)
onesigma_graph.SetMarkerStyle(20)
onesigma_graph.SetMarkerColor(kBlack)
twosigma_graph.SetFillColor(kYellow)
twosigma_graph.SetTitle('Closure test on A_{FB} (~900 toy datasets each)')

#for each of the data points tested
for i in range(len(input_values)) :
	#build the filepath
	filepath = './CLOSURE_TEST_AFB_'
	if input_values[i] < 0 :
		filepath+='-'
	filepath+=str(abs(input_values[i])).split('.')[0]+str(abs(input_values[i])).split('.')[1]
	filepath+='/'+filename
	print 'Getting A_FB values from file at path '+filepath
	#open the afb values file
	afb_values_file = open(filepath)
	#put the values in a list
	lines = afb_values_file.readlines()
	#find the relevant values
	twosigma_low = eval(lines[int(round(0.025*len(lines)))-1])
	onesigma_low = eval(lines[int(round(0.160*len(lines)))-1])
	mean 		 = eval(lines[int(round(0.500*len(lines)))-1])
	onesigma_hi  = eval(lines[int(round(0.840*len(lines)))-1])
	twosigma_hi  = eval(lines[int(round(0.975*len(lines)))-1])
	#Set graphs' point values
	onesigma_graph.SetPoint(i,input_values[i],mean)
	twosigma_graph.SetPoint(i,input_values[i],mean)
	onesigma_graph.SetPointError(i,x_errs[i],x_errs[i],abs(onesigma_low-mean),abs(onesigma_hi-mean))
	twosigma_graph.SetPointError(i,x_errs[i],x_errs[i],abs(twosigma_low-mean),abs(twosigma_hi-mean))

#Reset Axis ranges

#Build a legend
leg = TLegend(0.62,0.67,0.9,0.9)
leg.AddEntry(onesigma_graph,'Mean values','PL')
leg.AddEntry(onesigma_graph,'#pm 1 #sigma','F')
leg.AddEntry(twosigma_graph,'#pm 2 #sigma','F')

#Plot the graphs
canv = TCanvas('canv','canv',900,900)
canv.cd()
twosigma_graph.Draw('A E4')
onesigma_graph.Draw('SAME E4')
onesigma_graph.Draw('SAME PLX')
leg.Draw()

canv.Write()
outfile.Close()
