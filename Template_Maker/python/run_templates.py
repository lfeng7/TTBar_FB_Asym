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
parser.add_option('--include_PDF',type='string', action='store', default='yes',		 dest='include_PDF',help='Whether or not to include PDF uncertainty in building templates')
parser.add_option('--include_JEC',type='string', action='store', default='yes',		 dest='include_JEC',help='Whether or not to include JEC systematics in building templates')
(options, args) = parser.parse_args()

SUM_CHARGES = False
if options.sum_charges.lower() == 'yes' :
	SUM_CHARGES = True
INCLUDE_PDF = False
if options.include_PDF.lower() == 'yes' :
	INCLUDE_PDF = True
INCLUDE_JEC = False
if options.include_JEC.lower() == 'yes' :
	INCLUDE_JEC = True

#Start up the output file
output_name = options.out_name
parfilename = options.parameters
if '.txt' not in parfilename :
	parfilename+='.txt'
#Start up the group of templates
print 'Creating template group'
templates = template_group(output_name,parfilename,SUM_CHARGES,INCLUDE_PDF,INCLUDE_JEC)
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
if options.parameters.find('final') == -1 :
	#Save all the templates for this group
	print 'Writing templates to file'
	templates.write_to_files()
	print 'Done'
#Make and save comparison plots if necessary
if options.parameters.find('final') != -1 or options.parameters.find('refined') != -1 :
	print 'Making comparison plots'
	templates.make_plots()
	print 'Done'