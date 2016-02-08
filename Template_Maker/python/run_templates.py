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
parser.add_option('--on_grid', 	  type='string', action='store', default='no',	 dest='on_grid',		help='Path to on_grid file holding list of files to run on')
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
input_file_path =''
if options.on_grid == 'yes' :
	input_file_path+='./tardir/'
else:
	input_file_path+='./'
input_file_path+=options.input
if not '.txt' in options.input :
	input_file_path+='.txt'
input_file = open(input_file_path,'r')
#Make appropriate templates from each line in the file
print 'Adding all files to distributions'
for line in input_file :
	if not line.startswith('#') :
		ttree_file_path = line.rstrip()
		print '	Adding files from '+ttree_file_path
		templates.add_file_to_distributions(ttree_file_path)
		print '	Done'
print 'Done'
#Build the NTMJ template
print 'Building NTMJ templates'
templates.build_NTMJ_templates()
print 'Done'
#Build the templates
print 'Building templates from distributions'
templates.build_templates()
print 'Done'
#Save all the templates for this group
print 'Writing templates to files'
templates.write_to_files()
print 'Done'
#Make and save the theta feed file
print 'Making theta feed file from separate template files'
step = options.parameters.split('_')[0]
if step == 'initial' or step == 'final' :
	separate_files = ['nominal','simple_systematics','fit_parameters','JEC','PDF_systematics','NTMJ']
	for i in range(len(separate_files)) :
		tmp = output_name+separate_files[i]
		separate_files[i] = tmp
else : 
	iname = output_name.split('_refined')[0]+'_initial'+output_name.split('_refined')[1]
	separate_files = [iname+'nominal',iname+'simple_systematics',output_name+'fit_parameters',iname+'JEC',iname+'PDF_systematics',output_name+'NTMJ']
cmd = 'hadd -f theta_feed'+output_name.split('templates')[1].rstrip('_')+'.root'
for separate_file in separate_files :
	cmd += ' '+separate_file+'.root'
print '	'+cmd
os.system(cmd)
print 'Done'
#Make and save comparison plots if necessary
if options.parameters.find('final') != -1 :
	print 'Making comparison plots'
	templates.make_plots()
	print 'Done'