from ROOT import *
import os
import copy
#Filtering functions
#keeps all the generated templates
def everything(name) :
    return True
#tosses anything that isn't a fit parameter
def no_systematics(name) :
    namesplit = name.split('__')
    if len(namesplit)==2 or (len(namesplit)>2 and namesplit[2].split('_')[0]=='par') :
        return True
    else :
    	print '	REMOVING HISTOGRAM '+name
        return False
#return false for the 'mu' and 'd' parameter histograms
def fix_mu_and_d(name) :
	namesplit = name.split('__')
	if len(namesplit)>2 and (namesplit[2]=='par_mu' or namesplit[2]=='par_d') :
		print '	REMOVING HISTOGRAM '+name
		return False
	else :
		return True
#return false for the 'Afb' and 'd' parameter histograms
def fix_AFB_and_d(name) :
	namesplit = name.split('__')
	if len(namesplit)>2 and (namesplit[2]=='par_Afb' or namesplit[2]=='par_d') :
		print '	REMOVING HISTOGRAM '+name
		return False
	else :
		return True
#return false for the 'Afb' and 'mu' parameter histograms
def fix_AFB_and_mu(name) :
	namesplit = name.split('__')
	if len(namesplit)>2 and (namesplit[2]=='par_Afb' or namesplit[2]=='par_mu') :
		print '	REMOVING HISTOGRAM '+name
		return False
	else :
		return True
#only allows A_FB and nuisance parameters to float
def only_AFB(name) :
	return (no_systematics(name) and fix_mu_and_d(name))
#only allows mu and nuisance parameters to float
def only_mu(name) :
	return (no_systematics(name) and fix_AFB_and_d(name))
#only allows d and nuisance parameters to float
def only_d(name) :
	return (no_systematics(name) and fix_AFB_and_mu(name))

#Build templates with the initial parameter value predictions
#All the input options
TEMPLATE_FILE_NAME = 'templates'
SUM_CHARGES = False
INCLUDE_PDF = True
INCLUDE_JEC = True
INITIAL_PARAMETERS_FILE = 'initial_parameters.txt'
REFINED_PARAMETERS_FILE = 'refined_parameters.txt'
FINAL_PARAMETERS_FILE   = 'final_parameters.txt'
#ON_GRID = False
ON_GRID = True
FILTER_FUNCTION = everything
#FILTER_FUNCTION = fix_mu_and_d
#FILTER_FUNCTION = fix_AFB_and_d
#FILTER_FUNCTION = fix_AFB_and_mu
#FILTER_FUNCTION = only_AFB
#FILTER_FUNCTION = only_mu
#FILTER_FUNCTION = only_d
USE_TOYS = False
#USE_TOYS = True
N_TOYS   = 1700
AFB_CENTRAL_VALUE = 0.1
MU_CENTRAL_VALUE  = 0.0
D_CENTRAL_VALUE   = 0.0
#STEP = 'all'
STEP = 'initial_templates'
#STEP = 'initial_fit'
#STEP = 'refined_templates'
#STEP = 'refined_fit'
#STEP = 'final_plots'
#Global parameters list
AFB_sigma = min(min(abs(1.-AFB_CENTRAL_VALUE)/2,abs(AFB_CENTRAL_VALUE+1.)/2),0.2)
parameters = []
parameters.append(['Rbck',  1.00, 0.90,  1.10])
parameters.append(['Rntmj', 1.00, 0.90,  1.10])
parameters.append(['Rqqbar',1.00, 0.90,  1.10])
parameters.append(['Afb',   AFB_CENTRAL_VALUE, AFB_CENTRAL_VALUE-AFB_sigma, AFB_CENTRAL_VALUE+AFB_sigma])
parameters.append(['mu',	MU_CENTRAL_VALUE, MU_CENTRAL_VALUE-1.00, MU_CENTRAL_VALUE+1.00])
parameters.append(['d',		D_CENTRAL_VALUE, D_CENTRAL_VALUE-1.00, D_CENTRAL_VALUE+1.00])
parameters.append(['scale', 1.00, 0.90,  1.10])

sum_charges = 'no'
if SUM_CHARGES :
	sum_charges = 'yes'
include_PDF = 'no'
if INCLUDE_PDF :
	include_PDF = 'yes'
include_JEC = 'no'
if INCLUDE_JEC :
	include_JEC = 'yes'

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
		newline = par[0].ljust(14)+' '
		newline+= str(par[1]).ljust(17)+' '
		newline+= str(par[2]).ljust(16)+' '
		newline+= str(par[3]).ljust(14)+' '
		if name == FINAL_PARAMETERS_FILE and len(par)>4 :
			newline+= str(par[4]).ljust(10)
		newline+='\n'
		newfile.write(newline)
	newfile.write('#\n')
	newfile.close()
	print 'New parameter file: '
	os.system('cat '+name)

def run_fit(templatefilename) :
	thetafeedname = templatefilename+'.root'
	print 'running fit with template file '+thetafeedname+'. . .'
	#build the model
	model = build_model_from_rootfile(thetafeedname, histogram_filter=FILTER_FUNCTION, include_mc_uncertainties=True)
	#take care of the zero bins
	model.fill_histogram_zerobins(0.00001)
	#copy the model's distribution to use for toy generation, setting the parameters of interest precisely
	dist_for_toys = copy.deepcopy(model.distribution)
	for p in dist_for_toys.get_parameters() :
		if p=='par_Afb' or p=='par_mu' or p=='par_d' :
			dist_for_toys.set_distribution_parameters(p,width=0.000000001,range=[-1.0,1.0])
	#alter the parameter priors and ranges for fitting
	for p in model.distribution.get_parameters() :
		if p=='par_Afb' :
			low_afb = (-1.0-AFB_CENTRAL_VALUE)/AFB_sigma
			hi_afb  = (1.0-AFB_CENTRAL_VALUE)/AFB_sigma
			model.distribution.set_distribution_parameters(p,typ='flat_distribution',range=[low_afb,hi_afb])
		elif p=='par_mu' or p=='par_Afb' :
			model.distribution.set_distribution_parameters(p,typ='flat_distribution',range=['-inf','inf'])
		else :
			d = model.distribution.get_distribution(p)
			if d['typ'] == 'gauss' :
				model.distribution.set_distribution_parameters(p, range = [-5.0, 5.0])
	#Set some options
	options = Options()
	options.set('global','debug','True')
	options.set('minimizer','strategy','robust')
	options.set('minimizer','minuit_tolerance_factor','10')
	if USE_TOYS and N_TOYS>1 :
		nthreads = str(N_TOYS/10)
		options.set('main','n_threads',nthreads)
	#run the fit
	signal_process_groups = {'': []}
	if USE_TOYS :
		parVals = mle(model, input = 'toys:0.0', n=N_TOYS, with_covariance = False, nuisance_prior_toys=dist_for_toys, signal_process_groups = signal_process_groups)
	else :
		parVals = mle(model, input = 'data', n=1, with_covariance = True, signal_process_groups = signal_process_groups)
	print 'PARVALS returned from mle = '+str(parVals)
	#if there was more than one run, save the parameter of interest values
	toy_parameter_values = []
	parameter_of_interest = 'par_Afb'
	if FILTER_FUNCTION == only_mu or FILTER_FUNCTION == fix_AFB_and_d :
		parameter_of_interest = 'par_mu'
	elif FILTER_FUNCTION == only_d or FILTER_FUNCTION == fix_AFB_and_mu :
		parameter_of_interest = 'par_d'
	for i in range(len(parVals[''][parameter_of_interest])) :
		sigma = parameters[3][1]-parameters[3][2]
		v     = parameters[3][1]+parVals[''][parameter_of_interest][i][0]*sigma
		print '	Found parameter value '+parameter_of_interest+' = '+str(v)
		toy_parameter_values.append(v)
	if len(toy_parameter_values)>1 :
		toy_parameter_values.sort()
		par_value_filename = parameter_of_interest+'_values_'
		if FILTER_FUNCTION == everything :
			par_value_filename+='everything'
		elif FILTER_FUNCTION == fix_mu_and_d :
			par_value_filename+='fix_mu_and_d'
		elif FILTER_FUNCTION == only_AFB :
			par_value_filename+='only_AFB'
		elif FILTER_FUNCTION == fix_AFB_and_d :
			par_value_filename+='fix_AFB_and_d'
		elif FILTER_FUNCTION == only_mu :
			par_value_filename+='only_mu'
		elif FILTER_FUNCTION == fix_AFB_and_mu :
			par_value_filename+='fix_AFB_and_mu'
		elif FILTER_FUNCTION == only_d :
			par_value_filename+='only_d'
		par_value_filename+='.txt'
		toy_parameter_values_file = open(par_value_filename,'w')
		for v in toy_parameter_values :
			toy_parameter_values_file.write(str(v)+'\n')
	parameter_values = {}
	#get back the results
	newpars = []
	print 'results of fit: '
	print 'NLL = '+str(parVals['']['__nll'][0])
	for p in model.get_parameters([]):
		parameter_values[p] = parVals[''][p][0][0]
		newpars.append([])
		print '-------------------------------------------------------'
		print 'name = '+p
		newpars[len(newpars)-1].append(p)
		print 'true value from theta = '+str(parVals[''][p][0][0])
		newpars[len(newpars)-1].append(parVals[''][p][0][0])
		print 'true sigma from theta = '+str(parVals[''][p][0][1])
		newpars[len(newpars)-1].append(parVals[''][p][0][1])
	print '-------------------------------------------------------'
	print 'length of parameter array = '+str(len(parVals[''][parameter_of_interest]))
	#Save fit histograms
	histos = evaluate_prediction(model,parameter_values,include_signal = False)
	write_histograms_to_rootfile(histos,'postfit_histos_'+thetafeedname)
	#Reset the global list of parameters
	for i in range(len(parameters)) :
		for j in range(len(newpars)) :
			if newpars[j][0] == parameters[i][0] or newpars[j][0].lstrip('par_') == parameters[i][0] :
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
	if ON_GRID :
		cmd+='./tardir/'
	else :
		cmd+='../../../Template_Maker/python/'
	cmd  += 'run_templates.py --parameters '+FINAL_PARAMETERS_FILE +' '
	cmd += '--out_name '+TEMPLATE_FILE_NAME+'_plots --sum_charges '+sum_charges+' --include_PDF '+include_PDF+' --include_JEC '+include_JEC+''
	if ON_GRID :
		cmd+=' --on_grid yes'
	os.system(cmd)

initial_templates_filename = TEMPLATE_FILE_NAME+'_initial'
refined_templates_filename = TEMPLATE_FILE_NAME+'_refined'

if STEP == 'initial_templates' or STEP == 'all' :
	#Build the initial parameters file
	print 'Building initial parameters file. . .'
	build_parameter_file(INITIAL_PARAMETERS_FILE)
	#Build the command, run the initial templates, and figure out the name of the file they're in
	cmd = 'python '
	if ON_GRID :
		cmd+='./tardir/'
	else :
		cmd+='../../../Template_Maker/python/'
	cmd  += 'run_templates.py --parameters '+INITIAL_PARAMETERS_FILE +' '
	cmd += '--out_name '+initial_templates_filename+' --sum_charges '+sum_charges+' --include_PDF '+include_PDF+' --include_JEC '+include_JEC+''
	if ON_GRID :
		cmd+=' --on_grid yes'
	os.system(cmd)

if STEP == 'initial_fit' or STEP == 'all' :
	#Run the fit the first time with the initial parameter guesses
	run_fit(initial_templates_filename)
	#Build another input parameter file with the new parameter values
	build_parameter_file(REFINED_PARAMETERS_FILE)

if STEP == 'refined_templates' or STEP == 'all' :
	#Build new templates
	cmd = 'python '
	if ON_GRID :
		cmd+='./tardir/'
	else :
		cmd+='../../../Template_Maker/python/'
	cmd  += 'run_templates.py --parameters '+REFINED_PARAMETERS_FILE +' '
	cmd += '--out_name '+refined_templates_filename+' --sum_charges '+sum_charges+' --include_PDF '+include_PDF+' --include_JEC '+include_JEC+''
	if ON_GRID :
		cmd+=' --on_grid yes'
	os.system(cmd)

if STEP == 'refined_fit' or STEP == 'all' :
	#Run the fit again
	run_fit(refined_templates_filename)
	#Calculate the final parameter values and put them in a file
	build_parameter_file(FINAL_PARAMETERS_FILE)

if STEP == 'final_plots' or STEP == 'all' :
	#Make comparison plots
	make_comparison_plots()