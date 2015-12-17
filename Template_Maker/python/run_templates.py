#run_template.py is the helper code that reads the input file for a set of templates and runs the
#template_maker workhorse code to produce the templates needed
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

from optparse import OptionParser
import os
from template_maker import template_group

##########								Parser Options								##########

parser = OptionParser()
parser.add_option('--input', 	  type='string', action='store', default='input',	 dest='input',		help='Path to input file holding list of files to run on')
parser.add_option('--parameters', type='string', action='store', default='initial_parameters',	 		dest='parameters',		help='Path to input file holding list of fitting parameters and their initial values')
parser.add_option('--out_name',	  type='string', action='store', default='templates',dest='out_name',   help='Name of output file that will have all the templates in it')
parser.add_option('--sum_charges',type='string', action='store', default='no',		 dest='sum_charges',help='Whether or not to integrate over the lepton charge in building templates')
(options, args) = parser.parse_args()

#Start up the output file
output_name = options.out_name
if options.sum_charges.lower() == 'yes' :
	output_name+= '_charge_summed_'
else :
	output_name+='_'
parfilename = options.parameters
if '.txt' not in parfilename :
	parfilename+='.txt'
#Start up the group of templates
print 'Creating template group'
templates = template_group(output_name,parfilename)
print 'Done'
#Open the input file
input_file_path ='./'
input_file_path+=options.input
if not '.txt' in options.input :
	input_file_path+='.txt'
input_file = open(input_file_path,'r')
#Make appropriate templates from each line in the file
print 'Adding all files to distributions'
for line in input_file :
	if not line.startswith('#') :
		ttree_dir_path = line.rstrip()
		print '	Adding files from '+ttree_dir_path
		templates.add_file_to_distributions(ttree_dir_path)
		print '	Done'
print 'Done'
#Build the templates
print 'Building templates from distributions'
templates.build_templates()
print 'Done'
#Build the NTMJ template
print 'Building NTMJ templates'
templates.build_NTMJ_templates()
print 'Done'
#Save all the templates for this group
print 'Writing templates to files'
templates.write_to_files()
print 'Done'
#Make and save the theta feed file
'Making theta feed file from separate template files'
step = options.parameters.split('_')[0]
if step == 'initial' or step == 'final' :
	separate_files = ['nominal_'+step,'simple_systematics_'+step,'fit_parameters_'+step,'JEC_'+step,'PDF_systematics_'+step,'NTMJ_'+step]
else : 
	separate_files = ['nominal_initial','simple_systematics_initial','fit_parameters_'+step,'JEC_initial','PDF_systematics_initial','NTMJ_'+step]
cmd = 'hadd -f theta_feed'+output_name.lstrip('templates')+step+'.root'
for separate_file in separate_files :
	cmd += ' '+output_name+separate_file+'.root'
print '	'+cmd
os.system(cmd)
print 'Done'
#Make and save comparison plots if necessary
if options.parameters.find('final') != -1 :
	print 'Making comparison plots'
	templates.make_plots()
	print 'Done'