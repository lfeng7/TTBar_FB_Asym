from fit import *

#build the model
model = build_model_from_rootfile(THIS_TEMPLATE_FILE_NAME, include_mc_uncertainties=True)
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
		if parameters[i][0] == newpars[j][0] :
			sigma = parameters[i][1]-parameters[i][2]
			truesigma = sigma*newpars[j][2]
			#if it's the last fit also reset the central values and put the sigma in the list too for good measure
			if THIS_TEMPLATE_FILE_NAME == refined_templates_filename :
				parameters[i][1] = parameters[i][1]+sigma*newpars[j][1]
				parameters[i].append(truesigma)
			parameters[i][2] = parameters[i][1]-truesigma
			parameters[i][3] = parameters[i][1]+truesigma
			break
#Write the details to the html file
model_summary(model)
report.write_html(THIS_TEMPLATE_FILE_NAME.rstrip('.root')+'_htmlout')