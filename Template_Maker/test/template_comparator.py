from ROOT import *

filename = 'templates_both.root'

x_canv = TCanvas('x_canv','x_canv',1200,900)
y_canv = TCanvas('y_canv','y_canv',1200,900)
z_canv = TCanvas('z_canv','z_canv',1200,900)

tfile = TFile(filename)

all_histo_names = []; all_histo_colors = []
all_histo_names.append('fqqs_minus_x'); all_histo_names.append('fqqs_minus_y'); all_histo_names.append('fqqs_minus_z')
all_histo_colors.append(kBlack); all_histo_colors.append(kBlack); all_histo_colors.append(kBlack)
all_histo_names.append('fqqa_minus_x'); all_histo_names.append('fqqa_minus_y'); all_histo_names.append('fqqa_minus_z')
all_histo_colors.append(kBlue); all_histo_colors.append(kBlue); all_histo_colors.append(kBlue)
all_histo_names.append('fgg_minus_x');  all_histo_names.append('fgg_minus_y');  all_histo_names.append('fgg_minus_z')
all_histo_colors.append(kRed); all_histo_colors.append(kRed); all_histo_colors.append(kRed)
all_histo_names.append('fbck_minus_x'); all_histo_names.append('fbck_minus_y'); all_histo_names.append('fbck_minus_z')
all_histo_colors.append(kGreen); all_histo_colors.append(kGreen); all_histo_colors.append(kGreen)

all_histos = []
for histo_name in all_histo_names :
	all_histos.append(tfile.Get(histo_name))

for i in range(len(all_histos)) :
	all_histos[i].SetLineWidth(4)
	all_histos[i].SetLineColor(all_histo_colors[i])
	all_histos[i].SetFillColor(all_histo_colors[i])
	all_histos[i].SetMarkerStyle(21)
	all_histos[i].SetFillStyle(0)
	all_histos[i].SetLineStyle(1)

x_min = 0.; x_max = 0.
y_min = 0.; y_max = 0.
z_min = 0.; z_max = 0.

for i in range(len(all_histos)/3) :
	if all_histos[3*i].GetMaximum()>x_max :
		x_max = all_histos[3*i].GetMaximum()
	if all_histos[3*i].GetMinimum()<x_min :
		x_min = all_histos[3*i].GetMinimum()
	if all_histos[3*i+1].GetMaximum()>y_max :
		y_max = all_histos[3*i+1].GetMaximum()
	if all_histos[3*i+1].GetMinimum()<y_min :
		y_min = all_histos[3*i+1].GetMinimum()
	if all_histos[3*i+2].GetMaximum()>z_max :
		z_max = all_histos[3*i+2].GetMaximum()
	if all_histos[3*i+2].GetMinimum()<z_min :
		z_min = all_histos[3*i+2].GetMinimum()

x_canv.cd()
all_histos[0].SetMaximum(1.02*x_max)
all_histos[0].SetMinimum(1.02*x_min)
all_histos[0].Draw()
y_canv.cd()
all_histos[1].SetMaximum(1.02*y_max)
all_histos[1].SetMinimum(1.02*y_min)
all_histos[1].Draw()
z_canv.cd()
all_histos[2].SetMaximum(1.02*z_max)
all_histos[2].SetMinimum(1.02*z_min)
all_histos[2].Draw()

for i in range(3,len(all_histos)) :
	if '_x' in all_histo_names[i] :
		x_canv.cd()
		all_histos[i].Draw('SAME')
	if '_y' in all_histo_names[i] :
		y_canv.cd()
		all_histos[i].Draw('SAME')
	if '_z' in all_histo_names[i] :
		z_canv.cd()
		all_histos[i].Draw('SAME')

x_leg = TLegend(0.62,0.67,0.9,0.9)
y_leg = TLegend(0.62,0.67,0.9,0.9)
z_leg = TLegend(0.62,0.67,0.9,0.9)
for i in range(len(all_histo_colors)) :
	if '_x' in all_histo_names[i] :
		x_leg.AddEntry(all_histos[i],all_histo_names[i],"F")
	if '_y' in all_histo_names[i] :
		y_leg.AddEntry(all_histos[i],all_histo_names[i],"F")
	if '_z' in all_histo_names[i] :
		z_leg.AddEntry(all_histos[i],all_histo_names[i],"F")

x_canv.cd(); x_leg.Draw()
y_canv.cd(); y_leg.Draw()
z_canv.cd(); z_leg.Draw()

outfile = TFile('template_comparator_plots_minus.root','recreate')
outfile.cd()
x_canv.Write()
y_canv.Write()
z_canv.Write()
