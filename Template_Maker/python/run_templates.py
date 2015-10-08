#run_template.py is the helper code that reads the input file for a set of templates and runs the
#template_maker workhorse code to produce the templates needed
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

from optparse import OptionParser
from template_maker import template_file

##########								Parser Options								##########

parser = OptionParser()
parser.add_option('--input', 	  type='string', action='store', default='input',	 dest='input',		help='Path to input file holding list of files to run on')
parser.add_option('--parameters', type='string', action='store', default='initial_parameters',	 		dest='parameters',		help='Path to input file holding list of fitting parameters and their initial values')
parser.add_option('--out_name',	  type='string', action='store', default='templates',dest='out_name',   help='Name of output file that will have all the templates in it')
parser.add_option('--sum_charges',type='string', action='store', default='no',		 dest='sum_charges',help='Whether or not to integrate over the lepton charge in building templates')
parser.add_option('--plots', 	  type='string', action='store', default='no',		 dest='plots', 		help='Whether or not to make comparison plots to data')
(options, args) = parser.parse_args()

#Start up the output file
output_name = options.out_name
if options.sum_charges.lower() == 'yes' :
	output_name+= '_charge_summed'
output_name += '.root'
parfilename = options.parameters
if '.txt' not in parfilename :
	parfilename+='.txt'
output_file = template_file(output_name,parfilename,options.sum_charges.lower())
#Open the input file
input_file_path ='./'
input_file_path+=options.input
if not '.txt' in options.input :
	input_file_path+='.txt'
input_file = open(input_file_path,'r')
#Make appropriate templates from each line in the file
for line in input_file :
	if line.startswith('#') :
		continue
	[ttree_dir_path,name,ifd] = line.rstrip().split()
	ifd = ifd.lower()
	output_file.addToDistributions(ttree_dir_path,name,ifd)
#Build the templates
output_file.build_templates()
#Build the NTMJ template
output_file.build_NTMJ_templates()
#Make and save comparison plots if necessary
if options.plots.lower() == 'yes' :
	output_file.make_plots()
#clean up after yourself
del output_file