from ROOT import *
import os

#Build templates with the initial parameter value predictions
#All the input options
TEMPLATE_FILE_NAME = 'templates'
SUM_CHARGES = 'yes'
INITIAL_PARAMETERS_FILE = 'initial_parameters.txt'
REFINED_PARAMETERS_FILE = 'refined_parameters.txt'
FINAL_PARAMETERS_FILE   = 'final_parameters.txt'
ON_GRID = 'no'
#Global parameters list
parameters = []
parameters.append(['Rbck',  1.00, 0.80,  1.20])
parameters.append(['Rntmj', 1.00, 0.80,  1.20])
parameters.append(['Rqqbar',1.00, 0.80,  1.20])
parameters.append(['Afb',   0.10, -0.40, 0.60])
parameters.append(['scale', 1.00, 0.80,  1.20])

def build_parameter_file(name) :
	newfile = open(name,'w')
	if name == INITIAL_PARAMETERS_FILE :
		newfile.write('#Initial parameter values for the first iteration of the fit\n')
	elif name == REFINED_PARAMETERS_FILE :
		newfile.write('#Refined parameter values for the second run of fitting\n')
	elif name == FINAL_PARAMETERS_FILE :
		newfile.write('#FINAL parameter values')
	newfile.write('#\n')
	if name == FINAL_PARAMETERS_FILE :
		newfile.write('#Parameter     Initial Value     "Down" Value     "Up" Value     Sigma     \n')
	else :	
		newfile.write('#Parameter     Initial Value     "Down" Value     "Up" Value     \n')
	newfile.write('#\n')
	for par in parameters :
		newline = par[0].ljust(15)
		newline+= str(par[1]).ljust(18)
		newline+= str(par[2]).ljust(17)
		newline+= str(par[3]).ljust(15)
		if name == FINAL_PARAMETERS_FILE :
			newline+= str(par[4]).ljust(10)
		newline+='\n'
		newfile.write(newline)
	newfile.write('#\n')
	newfile.close()
	print 'New parameter file: '
	os.system('cat '+name)

def run_fit(templatefilename) :
	thetafeedname = 'theta_feed'+templatefilename.split('templates')[1]
	print 'running fit with template file '+thetafeedname+'. . .'
	#build the model
	model = build_model_from_rootfile(thetafeedname, include_mc_uncertainties=True)
	#take care of the zero bins
	model.fill_histogram_zerobins()	
	#expand the available range of the parameters to 5 sigma
	for p in model.distribution.get_parameters():
	    model.distribution.set_distribution_parameters(p, range = [-5.0, 5.0])
	#run the fit
	signal_process_groups = {'': []}
	parVals = mle(model, input = 'data', n=1, signal_process_groups = signal_process_groups)
	#get back the results
	newpars = []
	print 'results of fit: '
	for p in model.get_parameters([]):
		newpars.append([])
		print '-------------------------------------------------------'
		print 'name = '+p
		newpars[len(newpars)-1].append(p)
		print 'true value from theta = '+str(parVals[''][p][0][0])
		newpars[len(newpars)-1].append(parVals[''][p][0][0])
		print 'true sigma from theta = '+str(parVals[''][p][0][1])
		newpars[len(newpars)-1].append(parVals[''][p][0][1])
	print '-------------------------------------------------------'
	#Reset the global list of parameters
	for i in range(len(parameters)) :
		for j in range(len(newpars)) :
			if newpars[j][0].find(parameters[i][0])!=-1 :
				sigma = parameters[i][1]-parameters[i][2]
				truesigma = sigma*newpars[j][2]
				#if it's the last fit also reset the central values and put the sigma in the list too for good measure
				if templatefilename == refined_templates_filename :
					parameters[i][1] = parameters[i][1]+sigma*newpars[j][1]
					parameters[i].append(truesigma)
				parameters[i][2] = parameters[i][1]-truesigma
				parameters[i][3] = parameters[i][1]+truesigma
				break
	#Write the details to the html file
	#model_summary(model)
	#report.write_html(templatefilename.rstrip('.root')+'_htmlout')

def make_comparison_plots() :
	#build the command and run the final template plot file
	#Build the command, run the initial templates, and figure out the name of the file they're in
	cmd = 'python '
	if ON_GRID=='yes':
		cmd+='./tardir/'
	else :
		cmd+='../../Template_Maker/python/'
	cmd  += 'run_templates.py --parameters '+FINAL_PARAMETERS_FILE +' '
	cmd += '--out_name '+TEMPLATE_FILE_NAME+'_plots --sum_charges '+SUM_CHARGES+''
	if ON_GRID=='yes' :
		cmd+=' --on_grid yes'
	os.system(cmd)

#Build the initial parameters file
print 'Building initial parameters file. . .'
build_parameter_file(INITIAL_PARAMETERS_FILE)
#Build the command, run the initial templates, and figure out the name of the file they're in
cmd = 'python '
if ON_GRID=='yes':
	cmd+='./tardir/'
else :
	cmd+='../../Template_Maker/python/'
cmd  += 'run_templates.py --parameters '+INITIAL_PARAMETERS_FILE +' '
cmd += '--out_name '+TEMPLATE_FILE_NAME+'_initial --sum_charges '+SUM_CHARGES+''
if ON_GRID=='yes' :
	cmd+=' --on_grid yes'
os.system(cmd)
initial_templates_filename = TEMPLATE_FILE_NAME+'_initial'
refined_templates_filename = TEMPLATE_FILE_NAME+'_refined'
if SUM_CHARGES == 'yes' :
	initial_templates_filename+= '_charge_summed'
initial_templates_filename += '.root'
#Run the fit the first time with the initial parameter guesses
run_fit(initial_templates_filename)
#Build another input parameter file with the new parameter values
build_parameter_file(REFINED_PARAMETERS_FILE)
#Build new templates
cmd = 'python '
if ON_GRID=='yes':
	cmd+='./tardir/'
else :
	cmd+='../../Template_Maker/python/'
cmd  += 'run_templates.py --parameters '+REFINED_PARAMETERS_FILE +' '
cmd += '--out_name '+TEMPLATE_FILE_NAME+'_refined --sum_charges '+SUM_CHARGES+''
if ON_GRID=='yes' :
	cmd+=' --on_grid yes'
os.system(cmd)
if SUM_CHARGES == 'yes' :
	refined_templates_filename+= '_charge_summed'
refined_templates_filename += '.root'
#Run the fit again
run_fit(refined_templates_filename)
#Calculate the final parameter values and put them in a file
build_parameter_file(FINAL_PARAMETERS_FILE)
#Make comparison plots
make_comparison_plots()