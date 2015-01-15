#run_template.py is the helper code that reads the input file for a set of templates and runs the
#template_maker workhorse code to produce the templates needed
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

from optparse import OptionParser
from template_maker import background_template, signal_template, template_file

##########								Parser Options								##########

parser = OptionParser()
parser.add_option('--input', 	  type='string', action='store', default='input',	 dest='input',		help='Path to input file holding list of files to run on')
parser.add_option('--on_grid', 	  type='string', action='store', default='no',		 dest='on_grid',	help='Changes everything to relative paths if running on the grid, default is "no"')
parser.add_option('--out_name',	  type='string', action='store', default='templates',dest='out_name',   help='Name of output file that will have all the templates in it')
parser.add_option('--sum_charges',type='string', action='store', default='no',		 dest='sum_charges',help='Whether or not to integrate over the lepton charge in building templates')
(options, args) = parser.parse_args()

#Start up the output file
output_name = options.out_name
if '.root' not in output_name :
	output_name += '.root'
output_file = template_file(output_name)
#Open the input file
input_file_path = ''
if options.on_grid=='yes' :
	input_file_path+='tardir/'
else :
	input_file_path+='.'
input_file_path+=options.input
if not '.' in options.input :
	input_file_path+='.txt'
input_file = open(input_file_path,'r')
#Make appropriate templates from each line in the file
for line in input_file :
	if line.startswith('#') :
		continue
	[ttree_filepath,name,ifd] = line.rstrip().split()
	output_file.addToTemplate(ttree_filepath,name,ifd.lower())
#Normalize the total distributions in preparation for fitting
output_file.normalizeDistributions()
#clean up after yourself
del output_file