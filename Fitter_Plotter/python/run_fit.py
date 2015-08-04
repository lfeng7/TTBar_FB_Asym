#run_fitter.py is the helper code that reads the input files of MC templates and data TTrees and runs the 
#fitter workhorse code to find the qqbar and background fractions and A_FB, also makes fit comparison plots
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

from optparse import OptionParser
from time import strftime
from fitter import fitter
from ROOT import *

##########								Parser Options								##########

parser = OptionParser()
parser.add_option('--run_name', 	  type='string',action='store',default='new_run',   dest='run_name', 	   help='Distinguishing name for this run')
parser.add_option('--templates_file', type='string',action='store',default='templates', dest='templates_file', help='Path to input file holding MC templates to fit to data')
parser.add_option('--on_grid',		  type='string',action='store',default='no',		dest='on_grid',		   help='Changes everything to relative paths if running on the grid, default is "no"')
parser.add_option('--out_name',		  type='string',action='store',default='fit_plots', dest='out_name',	   help='Name of output file that will have all the plots and stuff in it')
(options, args) = parser.parse_args()

#set up the input and output files and run name
input_file_path = ''
if options.on_grid=='yes' :
	input_file_path+='tardir/'
else :
	input_file_path+='./'
templates_filename = input_file_path+options.templates_file
if not templates_filename.endswith('.root') :
	templates_filename+='.root'
output_name = options.out_name
if not output_name.endswith('.root') :
	output_name += '.root'
runname = options.run_name
if runname == 'new_run' :
	runname+='_'+strftime('%Y-%m-%d_%X')
#make a new fitter
fit_obj = fitter(runname,options.on_grid,templates_filename,output_name)
#fit the MC templates to the data
fit_obj.fit()
#make plots
save_pdfs = True
fit_obj.makeComparisonPlots(save_pdfs)
	#more plots go here or whatever
#clean up after yourself
del fit_obj